from __future__ import annotations

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
SOURCE = HERE.parent / "20260915_geuvadis_chr18_stage2b_rc1"


def load(name: str):
    return json.loads((RESULTS / name).read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    h = hashlib.sha256(path.read_bytes())
    return h.hexdigest()


def main() -> None:
    binding = load("RUN_BINDING.json")
    tq1 = load("TQ1_RESULTS.json")
    e0 = load("E0_RESULTS.json")
    final = load("FINAL_STATUS.json")
    frozen_models = json.loads((SOURCE / "results" / "FROZEN_MODELS.json").read_text(encoding="utf-8"))
    expected_checkpoints = {x["gene_id"]: x["checkpoint_sha256"] for x in frozen_models}

    checks = {
        "config_hash_matches": binding["config_sha256"] == sha256(HERE / "CONFIG_FROZEN.json"),
        "protocol_hash_matches": binding["protocol_sha256"] == sha256(HERE / "FROZEN_PROTOCOL.zh-CN.md"),
        "source_final_hash_matches": binding["source_final_sha256"] == sha256(SOURCE / "results" / "FINAL_STATUS.json"),
        "source_panel_hash_matches": binding["source_trait_panel_sha256"] == sha256(SOURCE / "results" / "TRAIT_PANEL.json"),
        "tq1_is_ab_only": set(tq1["macro_r2"]) == {"A", "B"} and all(set(x) == {"trait_index", "r2_A", "r2_B"} for x in tq1["trait_metrics"]),
        "tq1_n_is_230": tq1["n_individuals"] == 230,
        "tq1_decision_recomputes": tq1["pass"] == (tq1["macro_delta_r2"] > 0 and tq1["paired_bootstrap_ci95"][0] > 0),
        "tq1_pool_consumed": tq1["pool_status_after_run"] == "CONSUMED_BY_TQ1_PROHIBITED_FOR_CDE_BRIDGE",
        "e0_has_eight_panels": len(e0["panels"]) == 8,
        "e0_n_is_280": e0["n_individuals"] == 280,
        "e0_checkpoint_hashes_match": all(expected_checkpoints[x["gene_id"]] == x["checkpoint_sha256"] for x in e0["panels"]),
        "e0_decision_recomputes": e0["pass"] == all(
            e0["co_primary"][name]["point"] > 0 and e0["co_primary"][name]["family_cluster_ci95"][0] > 0
            for name in ("ce_empirical_minus_real", "ce_donor_minus_real")
        ),
        "no_bridge_result_file": not (RESULTS / "BRIDGE_RESULTS.json").exists() and not (RESULTS / "CDE_RESULTS.json").exists(),
        "grassmann_not_evaluated": final["grassmann_evaluated"] is False,
        "combined_status_recomputes": final["status"] == "NOT_READY_TASK" and not tq1["pass"] and e0["pass"],
        "new_bridge_not_authorized": final["new_bridge_authorized"] is False,
    }
    output = {"status": "PASS" if all(checks.values()) else "FAIL", "checks": checks}
    (RESULTS / "VALIDATION.json").write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    if output["status"] != "PASS":
        raise SystemExit(json.dumps(output, indent=2))
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()

