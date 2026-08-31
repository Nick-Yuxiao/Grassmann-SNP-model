from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
REQUIRED = [
    "PROTOCOL_ADDENDUM.v7.3.0.md",
    "ESTIMAND_FAMILY_SPEC.v7.3.0.yaml",
    "DECISIONS.v7.3.0.tsv",
    "validate_v7_3_0.py",
    "p0/audit_estimand_family_v7_3_0.py",
    "p0/run_estimand_family_audit_v7_3_0.sh",
    "p0/tests/test_v7_3_0.py",
    "server_ops/SERVER_STEPS.v7.3.0.md",
]


def main() -> None:
    missing = [name for name in REQUIRED if not (ROOT / name).is_file()]
    if missing:
        raise SystemExit(f"missing files: {missing}")
    protocol = (ROOT / "PROTOCOL_ADDENDUM.v7.3.0.md").read_text(encoding="utf-8")
    spec = (ROOT / "ESTIMAND_FAMILY_SPEC.v7.3.0.yaml").read_text(encoding="utf-8")
    runner = (ROOT / "p0" / "run_estimand_family_audit_v7_3_0.sh").read_text(
        encoding="utf-8"
    )
    with (ROOT / "DECISIONS.v7.3.0.tsv").open(
        encoding="utf-8-sig", newline=""
    ) as handle:
        decisions = list(csv.DictReader(handle, delimiter="\t"))
    checks = {
        "decision_count": len(decisions) == 12,
        "all_frozen": {row["status"] for row in decisions} == {"FROZEN"},
        "efficiency_double_firewall": (
            "Positive and negative results are both descriptive" in protocol
            and 'efficiency_positive_capacity_upgrade: "FORBIDDEN"' in spec
            and 'efficiency_negative_capacity_upgrade: "FORBIDDEN"' in spec
        ),
        "multi_budget": "20k, 30k and 40k" in protocol,
        "ld_persistence": "ELIGIBLE_FOR_CONFIRMATORY_30K" in protocol,
        "longrange_budget": "INCONCLUSIVE_BUDGET_NOT_ESTIMABLE" in protocol,
        "data_seed_unit": 'independent_replication_unit: "DATA_SEED"' in spec,
        "sqrt2_not_bound": 'sqrt2_within_arm_as_strict_bound: false' in spec,
        "tail_post_hoc": 'tail_slope_registered_in_v7_2_4: false' in spec,
        "decay_firewall": 'constant_lr_gate_reusable_for_decay: false' in spec,
        "gpu_block": (
            'gpu_required: false' in spec
            and 'gpu_authorized: false' in spec
            and "nvidia-smi" not in runner
            and "CUDA_VISIBLE_DEVICES" not in runner
        ),
        "holdout_block": (
            'decision_holdout_access: "FORBIDDEN"' in spec
            and 'hgdp_access: "FORBIDDEN"' in spec
        ),
        "formal_block": (
            'formal_a1r_authorized: false' in spec
            and 'n5_data_seed_pilot_authorized: false' in spec
            and 'hapnest_authorized: false' in spec
        ),
        "source_manifests": all(
            name in runner
            for name in (
                "BUDGET_BRIDGE_MANIFEST.v7.2.1.sha256",
                "BUDGET_EXTENSION_MANIFEST.v7.2.3.sha256",
                "FINAL_BUDGET_MANIFEST.v7.2.4.sha256",
                "C0_EXTENSION_MANIFEST.v7.1.13.sha256",
            )
        ),
    }
    failed = [name for name, passed in checks.items() if not passed]
    if failed:
        raise SystemExit(f"validation failed: {failed}")
    print(
        json.dumps(
            {
                "status": "PASS",
                "protocol_version": "v7.3.0",
                "manifest_entries": len(REQUIRED),
                "role": "CPU_ONLY_ESTIMAND_AND_MASK_FAMILY_AUDIT",
                "efficiency_budget_points": [20000, 30000, 40000],
                "ld_block_capacity": "ELIGIBLE_FOR_CONFIRMATORY_30K",
                "longrange_capacity": "INCONCLUSIVE_BUDGET_NOT_ESTIMABLE",
                "independent_replication_unit": "DATA_SEED",
                "tail_slope_role": "POST_HOC_STOP_AND_PLANNING_ONLY",
                "gpu_authorized": False,
                "formal_a1r_authorized": False,
                "architecture_decision_permitted": False,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
