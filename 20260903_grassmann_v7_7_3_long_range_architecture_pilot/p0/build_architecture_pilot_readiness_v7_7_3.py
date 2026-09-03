from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def validate_source(source: dict) -> None:
    """Accept only an immutable v7.7.2 CPU task-validity PASS with no GPU or Grassmann."""
    required = {
        "status": "LONG_RANGE_TASK_VALIDITY_PASS",
        "protocol_version": "v7.7.2",
        "next_authorized_stage": "DRAFT_V7_7_3_LONG_RANGE_ARCHITECTURE_PILOT_NO_GPU",
    }
    for key, value in required.items():
        if source.get(key) != value:
            raise ValueError(f"source {key} mismatch: {source.get(key)!r} != {value!r}")
    gates = source.get("gates", {})
    if gates.get("all_pass") is not True:
        raise ValueError("source did not report all task-validity gates passing")
    for gate in (
        "local_negative_pass",
        "oracle_positive_pass",
        "conventional_global_positive_pass",
        "target_shuffled_negative_pass",
    ):
        if gates.get(gate) is not True:
            raise ValueError(f"source gate {gate} not passed")
    authorization = source.get("authorization", {})
    if authorization.get("gpu_used") is not False:
        raise ValueError("source unexpectedly used a GPU")
    if authorization.get("grassmann_fitted") is not False:
        raise ValueError("source unexpectedly fitted Grassmann")
    if authorization.get("architecture_decision_permitted") is not False:
        raise ValueError("source unexpectedly permits an architecture decision")


def load_arms(path: Path) -> list[dict]:
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def static_contract_checks(arms: list[dict]) -> dict:
    """Only structural checks; nothing is trained and no GPU is queried."""
    decision_cells = [row for row in arms if row["architecture_decision_eligible"] == "TRUE"]
    router_present = [row for row in decision_cells if row["shared_global_router"] == "PRESENT"]
    checks = {
        "decision_eligible_cell_count": len(decision_cells),
        "full_factorial_complete": sorted(
            (row["grassmann_component"], row["shared_global_router"]) for row in decision_cells
        ) == [("ABSENT", "ABSENT"), ("ABSENT", "PRESENT"), ("PRESENT", "ABSENT"), ("PRESENT", "PRESENT")],
        "router_present_inputs_identical": len({row["router_inputs"] for row in router_present}) == 1,
        "router_present_not_source_positions": all(
            row["router_inputs"] == "FULL_TOKEN_SEQUENCE_NOT_SOURCE_POSITIONS" for row in router_present
        ),
        "no_arm_execution_authorized": all(row["execution_authorized"] == "FALSE" for row in arms),
        "controls_present": {row["role"] for row in arms} >= {
            "TASK_VALIDITY_POSITIVE_CONTROL",
            "TARGET_SHUFFLED_NEGATIVE_CONTROL",
        },
    }
    checks["all_pass"] = (
        checks["decision_eligible_cell_count"] == 4
        and checks["full_factorial_complete"]
        and checks["router_present_inputs_identical"]
        and checks["router_present_not_source_positions"]
        and checks["no_arm_execution_authorized"]
        and checks["controls_present"]
    )
    return checks


def write_manifest(output_dir: Path, names: list[str]) -> None:
    lines = [f"{sha256(output_dir / name)}  ./{name}" for name in names]
    (output_dir / "ARCHITECTURE_PILOT_READINESS_MANIFEST.v7.7.3.sha256").write_text(
        "\n".join(lines) + "\n", encoding="utf-8", newline="\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-execution", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    if not args.source_execution.is_file():
        raise FileNotFoundError(args.source_execution)
    source = load_json(args.source_execution)
    validate_source(source)

    arms_path = ROOT / "ARCHITECTURE_PILOT_ARM_MAP.v7.7.3.tsv"
    arms = load_arms(arms_path)
    checks = static_contract_checks(arms)
    if not checks["all_pass"]:
        raise ValueError(f"static architecture-pilot contract checks failed: {checks}")

    if args.output_dir.exists():
        raise FileExistsError(f"refuse overwrite: {args.output_dir}")
    args.output_dir.mkdir(parents=True)

    payload = {
        "schema_version": "1.0",
        "protocol_version": "v7.7.3",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "LONG_RANGE_ARCHITECTURE_PILOT_CONTRACT_SIGNED_NO_GPU",
        "audit_role": "PROSPECTIVE_ARCHITECTURE_PILOT_CONTRACT_ONLY",
        "source": {
            "file": str(args.source_execution.resolve()),
            "sha256": sha256(args.source_execution),
            "status": source["status"],
            "task_validity": "CONFIRMED_PASS_SYNTHETIC_ONLY",
        },
        "static_contract_checks": checks,
        "factorial": {
            "grassmann_factor_levels": ["ABSENT", "PRESENT"],
            "global_router_factor_levels": ["ABSENT", "PRESENT"],
            "shared_router_not_given_source_positions": True,
            "primary_contrast": "NLL_LR_A01_MINUS_NLL_LR_A11",
            "practical_margin_nats_per_target": 0.010,
            "go_rule": "PAIRED_ONE_SIDED_95_LCB_STRICTLY_GT_0.010_AND_VALIDITY_AND_SHUFFLED_AND_FAIRNESS_PASS",
            "secondary_estimands": [
                "GRASSMANN_MAIN_EFFECT_OVER_ROUTER_LEVELS",
                "GRASSMANN_BY_ROUTER_INTERACTION",
            ],
            "secondary_cannot_override_primary": True,
        },
        "fairness_audit": {
            "nonidentity_audit_required": True,
            "realized_compute_audit_required": True,
            "matched_parameter_label_permitted": False,
            "matched_compute_label_permitted": False,
        },
        "replication": {
            "repeat_unit": "SYNTHETIC_TRUTH_SEED",
            "init_seeds_are_independent_replicates": False,
            "formal_truth_seed_count": None,
            "formal_truth_seed_count_status": "PENDING_BLINDED_VARIANCE_PILOT",
        },
        "blinded_variance_pilot": {
            "device": "cpu",
            "scale": "REDUCED_CPU_PROXY",
            "arm_means_released": False,
            "releasable_outputs": [
                "PRIMARY_CONTRAST_DISPERSION",
                "FROZEN_FORMAL_TRUTH_SEED_COUNT",
            ],
            "pilot_data_barred_from_formal_analysis": True,
        },
        "authorization": {
            "gpu_authorized": False,
            "grassmann_training_authorized": False,
            "cpu_blinded_variance_pilot_authorized": False,
            "architecture_decision_permitted": False,
            "a3_scaling_authorized": False,
            "hapnest_authorized": False,
            "hgdp_access": "FORBIDDEN",
            "phenotype_claim_permitted": False,
            "biological_long_range_claim_permitted": False,
        },
        "next_authorized_stage": "IMPLEMENT_V7_7_4_CPU_BLINDED_VARIANCE_PILOT_NO_GPU",
    }

    readiness_name = "ARCHITECTURE_PILOT_READINESS.v7.7.3.json"
    arm_name = "ARCHITECTURE_PILOT_ARM_MAP.v7.7.3.tsv"
    (args.output_dir / readiness_name).write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n"
    )
    (args.output_dir / arm_name).write_text(
        arms_path.read_text(encoding="utf-8"), encoding="utf-8", newline="\n"
    )
    write_manifest(args.output_dir, [readiness_name, arm_name])
    print(json.dumps({
        "status": payload["status"],
        "output_dir": str(args.output_dir.resolve()),
        "gpu_authorized": False,
        "grassmann_training_authorized": False,
        "next_authorized_stage": payload["next_authorized_stage"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
