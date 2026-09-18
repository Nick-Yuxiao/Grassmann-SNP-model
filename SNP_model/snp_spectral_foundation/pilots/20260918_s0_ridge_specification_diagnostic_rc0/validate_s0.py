"""Independent re-derivation of every S0 headline number.

Point ``--results`` at ``results_smoke/`` to check a synthetic rehearsal; with no
argument it checks the real ``results/``.

Tolerances are stated per quantity and sized to what the artefact actually
stores. TQ-G1 rc0 compared a rank statistic recomputed from a table written at
six decimals against a hard-coded 1e-9 and failed itself; nothing here is
checked more tightly than the file it is read from can support.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

PILOT_DIR = Path(__file__).resolve().parent
if str(PILOT_DIR) not in sys.path:
    sys.path.insert(0, str(PILOT_DIR))

from s0_core import SPECIFICATIONS  # noqa: E402

# RESULTS_S0.json stores full double precision, so re-derived sums agree to
# floating-point accumulation order, not to the last bit.
TOLERANCE = 1e-9


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", default=str(PILOT_DIR / "results"))
    args = parser.parse_args()
    result_dir = Path(args.results).resolve()
    smoke = result_dir.name != "results"

    results_path = result_dir / "RESULTS_S0.json"
    status_path = result_dir / "FINAL_STATUS_S0.json"
    binding_path = result_dir / "RUN_BINDING_S0.json"
    required = [results_path, status_path] if smoke else [results_path, status_path, binding_path]
    for path in required:
        if not path.exists():
            print(f"FAIL: missing artifact {path.name}")
            return 1

    results = json.loads(results_path.read_text(encoding="utf-8"))
    status = json.loads(status_path.read_text(encoding="utf-8"))
    config = json.loads((PILOT_DIR / "CONFIG_S0.json").read_text(encoding="utf-8"))

    checks: list[dict[str, object]] = []

    def check(name: str, ok: bool, detail: object = "") -> None:
        checks.append({"check": name, "pass": bool(ok), "detail": detail})

    macro = results["macro_r2"]
    counts = results["snp_counts"]

    # --- the contrasts are re-derived from the stored macro R2 --------------
    for specification in SPECIFICATIONS:
        for count in counts:
            entry = results["contrasts"][f"B-A@{specification}:{count}"]
            expected = macro[f"B@{specification}:{count}"] - macro[f"A@{specification}"]
            check(
                f"B-A@{specification}:{count} equals macro B minus macro A",
                abs(float(entry["macro_delta_r2"]) - expected) < TOLERANCE,
                {"reported": entry["macro_delta_r2"], "recomputed": expected},
            )
            check(
                f"B-A@{specification}:{count} agrees with its own bootstrap macro_r2",
                abs(float(entry["macro_r2"]["A"]) - macro[f"A@{specification}"]) < TOLERANCE
                and abs(float(entry["macro_r2"]["B"]) - macro[f"B@{specification}:{count}"])
                < TOLERANCE,
            )
            low, high = entry["paired_bootstrap_ci95"]
            check(f"B-A@{specification}:{count} interval is ordered", low <= high, [low, high])
            check(
                f"B-A@{specification}:{count} pass rule matches the frozen rule",
                bool(entry["pass"]) == bool(float(entry["macro_delta_r2"]) > 0 and low > 0),
            )
            check(
                f"B-A@{specification}:{count} resamples individuals and genes",
                entry["resampling_units"] == ["evaluation_individual", "gene"],
            )

    # --- S2 and S3 differ only in how the dosage block is fitted ------------
    check(
        "S2 and S3 agree on arm A",
        abs(macro["A@S2"] - macro["A@S3"]) < 1e-8,
        {"S2": macro["A@S2"], "S3": macro["A@S3"]},
    )

    # --- the verdict is the one the numbers imply ---------------------------
    primary = counts[-1] if smoke else int(config["snps_per_gene"])
    rows = {}
    for specification in SPECIFICATIONS:
        entry = results["contrasts"][f"B-A@{specification}:{primary}"]
        rows[specification] = bool(entry["pass"])
        reported = status["per_specification"][specification]
        check(
            f"FINAL_STATUS repeats B-A@{specification} faithfully",
            abs(float(reported["delta"]) - float(entry["macro_delta_r2"])) < TOLERANCE
            and bool(reported["clears_zero"]) == bool(entry["pass"]),
            reported,
        )
    expected_status = (
        "SPECIFICATION-ARTEFACT-SUPPORTED"
        if (rows["S2"] and rows["S3"] and not rows["S1"])
        else "SPECIFICATION-ARTEFACT-NOT-SUPPORTED"
    )
    check("status matches the frozen decision rule", status["status"] == expected_status, status["status"])
    check("the diagnostic does not claim to overturn TQ-G1", status["overturns_tqg1"] is False)

    # --- held-out and sealed-role hygiene -----------------------------------
    development = results["development"]
    check(
        "sealed roles were never read",
        set(development["sealed_roles_never_read"]) == set(config["forbidden_roles"])
        and set(status["sealed_roles_never_read"]) == set(config["forbidden_roles"]),
        development["sealed_roles_never_read"],
    )
    check(
        "only development roles were used",
        set(development["roles_used"]) == set(config["roles_used"]),
        development["roles_used"],
    )
    check(
        "every individual lands in exactly one fold",
        sum(development["fold_sizes"]) == development["n_individuals"]
        and len(development["fold_sizes"]) == development["cv_folds"],
        development["fold_sizes"],
    )
    check(
        "families outnumber folds",
        development["n_families"] >= development["cv_folds"],
        {"families": development["n_families"], "folds": development["cv_folds"]},
    )
    check("QC scope is train-only", results["qc"]["scope"] == "dev_train_only")
    check(
        "gene selection did not use phenotype signal",
        results["genes"]["selection_uses_phenotype_signal"] is False
        and config["gene_selection_uses_phenotype_signal"] is False
        and config["cis_snp_selection_uses_phenotype"] is False,
    )
    minimum = 1 if smoke else int(config["min_genes_required"])
    check(
        "gene count meets the minimum",
        results["genes"]["analysed"] >= minimum
        and results["genes"]["analysed"] == len(results["genes"]["gene_ids"]),
        {"analysed": results["genes"]["analysed"], "minimum": minimum},
    )
    check(
        "every specification and SNP count was actually run",
        set(results["specifications"]) == set(SPECIFICATIONS)
        and all(f"A@{s}" in macro for s in SPECIFICATIONS)
        and all(f"B@{s}:{c}" in macro for s in SPECIFICATIONS for c in counts),
    )
    check(
        "interpretation limits are carried through",
        set(config["interpretation_limits"]).issubset(set(results["interpretation_limits"])),
    )

    # --- artifact integrity -------------------------------------------------
    if smoke:
        print("(smoke mode: RUN_BINDING integrity checks are skipped)")
    else:
        binding = json.loads(binding_path.read_text(encoding="utf-8"))
        for name, key in (("s0_core.py", "core_sha256"), ("CONFIG_S0.json", "config_sha256")):
            recorded = binding.get(key)
            check(
                f"{name} is unchanged since the run",
                bool(recorded) and recorded == sha256_file(PILOT_DIR / name),
                {"recorded": recorded},
            )

    failures = [c for c in checks if not c["pass"]]
    report = {
        "validator": "s0-rc0",
        "n_checks": len(checks),
        "n_failures": len(failures),
        "status": "PASS" if not failures else "FAIL",
        "checks": checks,
    }
    (result_dir / "VALIDATION_S0.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    for entry in checks:
        print(f"[{'PASS' if entry['pass'] else 'FAIL'}] {entry['check']}")
    print(f"\n{report['status']}: {len(checks) - len(failures)}/{len(checks)} checks passed")
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
