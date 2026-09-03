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


def validate_source(source: dict) -> None:
    """Accept only an immutable v7.7.7 re-plan contract as the source."""
    if source.get("status") != "LONG_RANGE_DIFFICULTY_REPLAN_CONTRACT_SIGNED_NO_LAUNCH":
        raise ValueError("source must be the v7.7.7 re-plan contract readiness")
    if source.get("protocol_version") != "v7.7.7":
        raise ValueError("source protocol_version mismatch")
    if source.get("next_authorized_stage") != "IMPLEMENT_V7_7_8_LONG_RANGE_DIFFICULTY_REPLAN_HARNESS_NO_LAUNCH":
        raise ValueError("source next stage mismatch")
    authorization = source.get("authorization", {})
    if authorization.get("gpu_authorized") is not False:
        raise ValueError("source unexpectedly authorizes GPU")
    if authorization.get("grassmann_fitted") is not False:
        raise ValueError("source unexpectedly fitted Grassmann")


def load_arms(path: Path) -> list[dict]:
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def static_arm_checks(arms: list[dict]) -> dict:
    families = {row["family"] for row in arms}
    conventional = [row for row in arms if row["role"] == "conventional_comparator"]
    roles = {row["role"] for row in arms}
    checks = {
        "conventional_family_count": len(conventional),
        "has_grassmann_arm": any(row["arm_id"] == "GRASSMANN" for row in arms),
        "conventional_families": sorted(row["family"] for row in conventional),
        "controls_present": roles >= {
            "task_validity_positive_control",
            "local_insufficiency_control",
            "target_shuffled_negative_control",
        },
        "no_arm_execution_authorized": all(row["execution_authorized"] == "FALSE" for row in arms),
    }
    checks["all_pass"] = (
        checks["conventional_family_count"] >= 3
        and checks["has_grassmann_arm"]
        and checks["controls_present"]
        and checks["no_arm_execution_authorized"]
    )
    return checks


def write_manifest(output_dir: Path, names: list[str]) -> None:
    lines = [f"{sha256(output_dir / name)}  ./{name}" for name in names]
    (output_dir / "GRASSMANN_SEPARATION_PREREG_READINESS_MANIFEST.v7.8.0.sha256").write_text(
        "\n".join(lines) + "\n", encoding="utf-8", newline="\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--replan-readiness", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    if not args.replan_readiness.is_file():
        raise FileNotFoundError(args.replan_readiness)
    source = json.load(args.replan_readiness.open(encoding="utf-8"))
    validate_source(source)

    arms_path = ROOT / "SEPARATION_ARM_MAP.v7.8.0.tsv"
    arms = load_arms(arms_path)
    checks = static_arm_checks(arms)
    if not checks["all_pass"]:
        raise ValueError(f"static separation arm-map checks failed: {checks}")

    if args.output_dir.exists():
        raise FileExistsError(f"refuse overwrite: {args.output_dir}")
    args.output_dir.mkdir(parents=True)

    readiness = {
        "schema_version": "1.0",
        "protocol_version": "v7.8.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "LONG_RANGE_GRASSMANN_SEPARATION_PREREGISTRATION_SIGNED_NO_LAUNCH",
        "audit_role": "PROSPECTIVE_SEPARATION_PREREGISTRATION_ONLY",
        "track": "AFFIRMATIVE_A1_SEPARATION",
        "source_replan_readiness": {
            "file": str(args.replan_readiness.resolve()),
            "sha256": sha256(args.replan_readiness),
            "status": source["status"],
        },
        "static_arm_checks": checks,
        "estimand": {
            "name": "SEPARATION",
            "definition": "NLL_conv_best - NLL_grassmann",
            "positive_favors": "GRASSMANN",
            "lr1_required": False,
            "distinct_from_v7_7_0_incremental_contrast": True,
        },
        "validity": {
            "lr0_margin_nats": 0.010,
            "oracle_min_improvement_nats": 0.050,
            "target_shuffled_margin_nats": 0.010,
        },
        "primary_separation_go": {
            "conventional_fails_rule": "NLL_conv_best-NLL_local PAIRED_2S_90CI_INSIDE_+-0.010",
            "grassmann_solves_rule": "NLL_local-NLL_grassmann PAIRED_1S_95LCB_GT_0.050",
            "grassmann_reaches_oracle_rule": "PER_SEED NLL_grassmann-NLL_oracle <= 0.050",
            "seed_stable_required": True,
        },
        "fairness_audit": {
            "matched_parameter_audit_required": True,
            "matched_compute_audit_required": True,
            "realized_compute_audit_required": True,
            "conventional_convergence_compute_sufficiency_audit_required": True,
            "nonidentity_audit_required": True,
            "matched_parameter_label_permitted": False,
            "matched_compute_label_permitted": False,
        },
        "replication": {
            "repeat_unit": "SYNTHETIC_TRUTH_SEED",
            "truth_seeds": [78001, 78002, 78003, 78004, 78005, 78006],
            "init_seeds_nested_and_averaged": True,
            "effective_n": 6,
        },
        "blinding": {
            "grassmann_arm_means_withheld_until_validity_and_fairness_locked": True,
            "pilot_data_barred_from_formal_analysis": True,
            "result_dependent_selection_permitted": False,
        },
        "gpu_gate": {
            "authorize_gpu_only_if": "CPU_FIRST_PRIMARY_SEPARATION_GO_WITH_ALL_FAIRNESS_AND_CONTROLS_PASS",
            "gpu_authorized_here": False,
        },
        "route_closure_branch": {
            "trigger": "NO_SEPARATION_UNDER_FAIR_CONVERGED_CONVENTIONAL_SUITE_AT_CPU_PROXY",
            "action": "CLOSE_V7_GRASSMANN_PRIMARY_ROUTE_NO_TASK_EXPANSION_NO_GPU",
        },
        "authorization": {
            "gpu_authorized": False,
            "grassmann_fitted": False,
            "pilot_execution_authorized": False,
            "architecture_decision_permitted": False,
            "a3_scaling_authorized": False,
            "hapnest_authorized": False,
            "hgdp_access": "FORBIDDEN",
            "phenotype_claim_permitted": False,
            "biological_long_range_claim_permitted": False,
        },
        "next_authorized_stage": "IMPLEMENT_V7_8_1_CPU_GRASSMANN_SEPARATION_HARNESS_NO_LAUNCH",
    }
    readiness_name = "GRASSMANN_SEPARATION_PREREG_READINESS.v7.8.0.json"
    arm_name = "SEPARATION_ARM_MAP.v7.8.0.tsv"
    (args.output_dir / readiness_name).write_text(
        json.dumps(readiness, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n"
    )
    (args.output_dir / arm_name).write_text(
        arms_path.read_text(encoding="utf-8"), encoding="utf-8", newline="\n"
    )
    write_manifest(args.output_dir, [readiness_name, arm_name])
    print(json.dumps({
        "status": readiness["status"],
        "output_dir": str(args.output_dir.resolve()),
        "conventional_family_count": checks["conventional_family_count"],
        "gpu_authorized": False,
        "next_authorized_stage": readiness["next_authorized_stage"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
