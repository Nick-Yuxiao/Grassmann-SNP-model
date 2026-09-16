from __future__ import annotations

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
RC1 = HERE.parent / "20260915_tq2_ab_empirical_power_calibration_rc1"
SOURCE = HERE.parent / "20260915_geuvadis_chr18_stage2b_rc1"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    result = json.loads((RESULTS / "CALIBRATION_RESULTS.json").read_text(encoding="utf-8"))
    final = json.loads((RESULTS / "FINAL_STATUS.json").read_text(encoding="utf-8"))
    rc1 = json.loads((RC1 / "results" / "FINAL_STATUS.json").read_text(encoding="utf-8"))
    by_delta = {x["injected_increment"]: x for x in result["grid"]}
    script = (HERE / "run_calibration.py").read_text(encoding="utf-8")
    checks = {
        "config_hash_matches": result["binding"]["config_sha256"] == sha256(HERE / "CONFIG_FROZEN.json"),
        "protocol_hash_matches": result["binding"]["protocol_sha256"] == sha256(HERE / "FROZEN_PROTOCOL.zh-CN.md"),
        "rc1_hash_matches": result["binding"]["rc1_final_sha256"] == sha256(RC1 / "results" / "FINAL_STATUS.json"),
        "source_panel_hash_matches": result["binding"]["source_trait_panel_sha256"] == sha256(SOURCE / "results" / "TRAIT_PANEL.json"),
        "rc1_failure_preserved": rc1["status"] == "AB_DESIGN_UNDERPOWERED",
        "five_frozen_grid_points": set(by_delta) == {0.0, 0.01, 0.02, 0.05, 0.1},
        "one_hundred_seeds_each": all(len(x["seed_results"]) == 100 for x in result["grid"]),
        "five_family_disjoint_folds": len(result["folds"]) == 5 and sum(x["test_n"] for x in result["folds"]) == 445,
        "null_fpr_recomputes": result["null_false_positive_rate"] == by_delta[0.0]["empirical_pass_rate"],
        "minimum_power_recomputes": result["power_at_minimum_increment"] == by_delta[0.02]["empirical_pass_rate"],
        "status_recomputes": final["status"] == "AB_DESIGN_UNDERPOWERED" and result["power_at_minimum_increment"] < result["required_power"],
        "no_real_phenotype_access": result["real_phenotype_accessed"] is False and "load_expression_subset" not in script,
        "no_encoder_or_grassmann": result["encoder_or_grassmann_used"] is False and "LocalGenotypeEncoder" not in script and "Grassmann" not in script,
    }
    output = {"status": "PASS" if all(checks.values()) else "FAIL", "checks": checks}
    (RESULTS / "VALIDATION.json").write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    if output["status"] != "PASS":
        raise SystemExit(json.dumps(output, indent=2))
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()

