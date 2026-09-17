"""Independent validator for a completed TQ-G1 run.

It re-derives every headline number from GENE_TABLE.tsv instead of trusting
RESULTS.json, and it re-checks the design invariants that make the run
interpretable. It never touches the model and never opens an outcome.

    python validate_tqg1.py
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

PILOT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(PILOT_DIR))

from tqg1_core import spearman  # noqa: E402

RESULT_DIR = PILOT_DIR / "results"
TOLERANCE = 1e-6


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--results",
        default=str(RESULT_DIR),
        help="Directory to validate. Point it at results_smoke/ to check a synthetic run.",
    )
    args = parser.parse_args()
    result_dir = Path(args.results).resolve()
    smoke = result_dir.name != "results"

    checks: list[dict[str, object]] = []

    def check(name: str, ok: bool, detail: object = "") -> None:
        checks.append({"check": name, "pass": bool(ok), "detail": detail})

    results_path = result_dir / "RESULTS.json"
    table_path = result_dir / "GENE_TABLE.tsv"
    binding_path = result_dir / "RUN_BINDING.json"
    status_path = result_dir / "FINAL_STATUS.json"
    required = [results_path, table_path] if smoke else [results_path, table_path, binding_path, status_path]
    for path in required:
        if not path.exists():
            print(f"FAIL: missing artifact {path.name}")
            return 1

    results = json.loads(results_path.read_text(encoding="utf-8-sig"))
    binding = json.loads(binding_path.read_text(encoding="utf-8-sig")) if binding_path.exists() else {}
    status = (
        json.loads(status_path.read_text(encoding="utf-8-sig"))
        if status_path.exists()
        else {
            "status": "SUPPORTED" if results["primary"]["pass"] else "NOT-SUPPORTED",
            "is_burden_test": False,
            "sealed_roles_still_sealed": sorted(
                set(results["interpretation_limits"]) & set()
            )
            or sorted(set(json.loads((PILOT_DIR / "CONFIG_FROZEN.json").read_text(encoding="utf-8-sig"))["sealed_roles"])
                      - set(results.get("sealed_roles_opened") or [])),
        }
    )
    config = json.loads((PILOT_DIR / "CONFIG_FROZEN.json").read_text(encoding="utf-8-sig"))

    with table_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    # --- the headline numbers are recomputed from the per-gene table ---------
    arms = ["A", "B", "C_gene", "C_full"]
    macro = {arm: float(np.mean([float(r[f"r2_{arm}"]) for r in rows])) for arm in arms}
    for arm in arms:
        reported = float(results["primary"]["macro_r2"][arm])
        check(
            f"macro_r2[{arm}] recomputed from GENE_TABLE",
            abs(reported - macro[arm]) < TOLERANCE,
            {"reported": reported, "recomputed": macro[arm]},
        )

    delta = macro["C_gene"] - macro["B"]
    check(
        "macro_delta_r2 equals C_gene minus B",
        abs(float(results["primary"]["macro_delta_r2"]) - delta) < TOLERANCE,
        {"reported": results["primary"]["macro_delta_r2"], "recomputed": delta},
    )

    low, high = results["primary"]["paired_bootstrap_ci95"]
    check("bootstrap interval is ordered", low <= high, [low, high])
    check(
        "pass rule matches the frozen rule",
        bool(results["primary"]["pass"]) == bool(delta > 0 and low > 0),
        {"delta": delta, "ci_low": low, "reported_pass": results["primary"]["pass"]},
    )
    check(
        "FINAL_STATUS agrees with RESULTS",
        status["status"] == ("SUPPORTED" if results["primary"]["pass"] else "NOT-SUPPORTED"),
        status["status"],
    )

    # --- rehearsal statistics ----------------------------------------------
    target = np.asarray([float(r["target_T_g"]) for r in rows])
    model_score = np.asarray([float(r["score_model_C_gene"]) for r in rows])
    linear_score = np.asarray([float(r["score_linear_B"]) for r in rows])
    for name, score in (("model", model_score), ("linear", linear_score)):
        reported = float(results["gene_effect_rehearsal"][name]["observed_spearman"])
        check(
            f"rehearsal spearman[{name}] recomputed",
            abs(reported - spearman(score, target)) < 1e-9,
            {"reported": reported, "recomputed": spearman(score, target)},
        )
    check(
        "perturbation scores are non-negative",
        bool(np.all(model_score >= 0) and np.all(linear_score >= 0)),
        {"min_model": float(model_score.min()), "min_linear": float(linear_score.min())},
    )
    check(
        "the run does not claim to be a burden test",
        results["gene_effect_rehearsal"]["is_burden_test"] is False
        and status["is_burden_test"] is False,
    )

    # --- design invariants --------------------------------------------------
    minimum = 1 if smoke else int(config["min_genes_required"])
    check(
        "gene count meets the frozen minimum",
        len(rows) >= minimum and len(rows) == results["genes"]["analysed"],
        {"rows": len(rows), "minimum": minimum},
    )
    check(
        "gene selection did not use phenotype signal",
        results["genes"]["selection_uses_phenotype_signal"] is False
        and config["gene_selection_uses_phenotype_signal"] is False
        and config["cis_snp_selection_uses_phenotype"] is False,
    )
    check(
        "encoder is shared and carries no locus lookup",
        results["encoder"]["shared_across_genes"] is True
        and results["encoder"]["config"]["use_local_block_embedding"] is False
        and results["encoder"]["config"]["local_position_mode"] == "relative_continuous",
        results["encoder"]["config"]["local_position_mode"],
    )
    check(
        "encoder never saw phenotype",
        results["encoder"]["phenotype_seen_during_pretraining"] is False,
    )
    check(
        "development and evaluation are family disjoint",
        results["split"]["family_disjoint"] is True,
    )
    opened = set(results.get("sealed_roles_opened") or [])
    still_sealed = set(status.get("sealed_roles_still_sealed") or [])
    declared = set(config["sealed_roles"])
    check(
        "sealed-role bookkeeping is consistent",
        opened | still_sealed == declared and not (opened & still_sealed),
        {"opened": sorted(opened), "still_sealed": sorted(still_sealed)},
    )
    check(
        "QC scope is train-only",
        results["qc"]["scope"] == "dev_train_only",
    )
    check(
        "interpretation limits are carried through",
        set(config["interpretation_limits"]).issubset(set(results["interpretation_limits"])),
    )

    # --- artifact integrity -------------------------------------------------
    if smoke:
        print("(smoke mode: RUN_BINDING integrity checks are skipped)")
    else:
        for name, key in (("tqg1_core.py", "core_sha256"), ("CONFIG_FROZEN.json", "config_sha256")):
            path = PILOT_DIR / name
            recorded = binding.get(key)
            check(
                f"{name} is unchanged since the run",
                bool(recorded) and recorded == sha256_file(path),
                {"recorded": recorded},
            )
    checkpoint = result_dir / "ENCODER_FROZEN.pt"
    check(
        "frozen encoder checkpoint matches its recorded hash",
        checkpoint.exists() and results["encoder"]["checkpoint_sha256"] == sha256_file(checkpoint),
    )

    failures = [c for c in checks if not c["pass"]]
    report = {
        "validator": "tqg1-rc0",
        "n_checks": len(checks),
        "n_failures": len(failures),
        "status": "PASS" if not failures else "FAIL",
        "checks": checks,
    }
    (result_dir / "VALIDATION.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    for c in checks:
        print(f"[{'PASS' if c['pass'] else 'FAIL'}] {c['check']}")
    print(f"\n{report['status']}: {len(checks) - len(failures)}/{len(checks)} checks passed")
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
