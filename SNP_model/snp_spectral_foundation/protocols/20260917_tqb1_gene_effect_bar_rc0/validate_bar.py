"""Independent validator for a completed TQ-B1 run.

Re-derives the headline numbers from GENE_TABLE.tsv rather than trusting
RESULTS.json, and re-checks the invariants that make the run interpretable.

    python validate_bar.py --results results
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from tqb1_io import dump_json, load_json, sha256_file  # noqa: E402
from tqb1_core import spearman  # noqa: E402

# GENE_TABLE carries repr-precision floats, so agreement should be near exact.
TOLERANCE = 1e-9


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", default="results")
    parser.add_argument("--config", default=None,
                        help="Config the run used. Defaults to the frozen CONFIG.json.")
    parser.add_argument("--burden", default=None, help="Burden TSV, to recompute the bar.")
    args = parser.parse_args()
    result_dir = Path(args.results)

    for name in ("RESULTS.json", "GENE_TABLE.tsv", "FINAL_STATUS.json"):
        if not (result_dir / name).exists():
            print(f"FAIL: missing artifact {name}")
            return 1
    results = load_json(result_dir / "RESULTS.json")
    status = load_json(result_dir / "FINAL_STATUS.json")
    config_path = Path(args.config) if args.config else HERE / "CONFIG.json"
    config = load_json(config_path)
    with (result_dir / "GENE_TABLE.tsv").open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    checks: list[dict[str, object]] = []

    def check(name: str, ok: bool, detail: object = "") -> None:
        checks.append({"check": name, "pass": bool(ok), "detail": detail})

    passing = list(results["passing_traits"])
    check("at least one trait passed the detectability gate", bool(passing), passing)

    for trait in passing:
        recomputed = float(np.mean([float(r[f"r2_B_eval_{trait}"]) for r in rows]))
        reported = float(results["held_out_macro"][trait]["macro_r2"]["B"])
        check(f"macro R2(B) for {trait} recomputed from GENE_TABLE",
              abs(reported - recomputed) < TOLERANCE,
              {"reported": reported, "recomputed": recomputed})
        recomputed_a = float(np.mean([float(r[f"r2_A_eval_{trait}"]) for r in rows]))
        reported_a = float(results["held_out_macro"][trait]["macro_r2"]["A"])
        check(f"macro R2(A) for {trait} recomputed from GENE_TABLE",
              abs(reported_a - recomputed_a) < TOLERANCE,
              {"reported": reported_a, "recomputed": recomputed_a})
        delta = recomputed - recomputed_a
        check(f"macro delta for {trait} equals B minus A",
              abs(float(results["held_out_macro"][trait]["macro_delta_r2"]) - delta) < TOLERANCE,
              {"recomputed": delta})
        gate = results["detectability_gate"][trait]
        check(f"gate rule for {trait} matches the frozen rule",
              bool(gate["pass"]) == bool(gate["macro_delta_r2"] > 0 and gate["paired_bootstrap_ci95"][0] > 0),
              {"delta": gate["macro_delta_r2"], "ci": gate["paired_bootstrap_ci95"]})
        scores = np.asarray([float(r[f"score_{trait}"]) for r in rows])
        check(f"perturbation scores for {trait} are non-negative", bool(np.all(scores >= 0)),
              {"min": float(scores.min())})
        check(f"perturbation scores for {trait} are not all identical",
              len(np.unique(scores)) > 1, {"distinct": int(len(np.unique(scores)))})
        check(f"GENE_TABLE keeps full precision for {trait}",
              len({r[f"score_{trait}"] for r in rows}) == len(rows),
              {"distinct_strings": len({r[f"score_{trait}"] for r in rows}), "rows": len(rows)})

    for trait in results["detectability_gate"]:
        if not results["detectability_gate"][trait]["pass"]:
            check(f"failing trait {trait} was never evaluated", trait not in passing)

    check("evaluation was opened only if a trait passed",
          bool(status["evaluation_split_opened"]) == bool(passing))
    check("gene count meets the frozen minimum",
          len(rows) >= int(config["min_genes_required"]) and len(rows) == results["genes"]["analysed"],
          {"rows": len(rows), "minimum": config["min_genes_required"]})
    check("gene and SNP selection never used phenotype",
          results["genes"]["selection_uses_phenotype"] is False
          and config["cis_snp_selection_uses_phenotype"] is False)
    check("splits are group disjoint", results["splits"]["group_disjoint"] is True)
    check("the run does not claim to be a foundation-model test",
          "this_is_an_additive_baseline_measurement_not_a_foundation_model_test"
          in results["interpretation_limits"])
    check("no GPU was used", (load_json(result_dir / "RUN_BINDING.json")["runtime"]["gpu_used"] is False)
          if (result_dir / "RUN_BINDING.json").exists() else True)

    if (result_dir / "RUN_BINDING.json").exists():
        binding = load_json(result_dir / "RUN_BINDING.json")
        for path, key in ((HERE / "tqb1_core.py", "core_sha256"), (config_path, "config_sha256")):
            recorded = binding.get(key)
            check(f"{path.name} is unchanged since the run",
                  bool(recorded) and recorded == sha256_file(path), {"recorded": recorded})

    if args.burden:
        from tqb1_io import read_burden

        burden = read_burden(Path(args.burden))
        for trait in passing:
            if trait not in results["bar"] or "matched_null" not in results["bar"][trait]:
                continue
            have = [r for r in rows if (r["gene_id"], trait) in burden]
            s = np.asarray([float(r[f"score_{trait}"]) for r in have])
            beta = np.asarray([abs(burden[(r["gene_id"], trait)][0]) for r in have])
            reported = float(results["bar"][trait]["matched_null"]["observed_spearman"])
            check(f"bar Spearman for {trait} recomputed from GENE_TABLE and burden",
                  abs(reported - spearman(s, beta)) < TOLERANCE,
                  {"reported": reported, "recomputed": spearman(s, beta), "n": len(have)})

    failures = [c for c in checks if not c["pass"]]
    dump_json(result_dir / "VALIDATION.json",
              {"validator": "tqb1-rc0", "n_checks": len(checks), "n_failures": len(failures),
               "status": "PASS" if not failures else "FAIL", "checks": checks})
    for c in checks:
        print(f"[{'PASS' if c['pass'] else 'FAIL'}] {c['check']}")
    print(f"\n{'PASS' if not failures else 'FAIL'}: {len(checks) - len(failures)}/{len(checks)} checks passed")
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
