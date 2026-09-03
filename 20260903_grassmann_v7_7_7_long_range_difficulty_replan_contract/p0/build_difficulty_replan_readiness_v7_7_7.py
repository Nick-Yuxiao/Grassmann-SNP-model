from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate_source(source: dict) -> None:
    """Accept only an immutable v7.7.6 UNRESOLVED CPU execution as the source."""
    if source.get("status") != "LONG_RANGE_TASK_DIFFICULTY_UNRESOLVED":
        raise ValueError("v7.7.6 source must be LONG_RANGE_TASK_DIFFICULTY_UNRESOLVED")
    if source.get("selected_k") is not None:
        raise ValueError("v7.7.6 source unexpectedly selected a K")
    if source.get("protocol_version") != "v7.7.6":
        raise ValueError("source protocol_version mismatch")
    if source.get("next_authorized_stage") != "STOP_NO_AUTOMATIC_TASK_EXPANSION_OR_GPU":
        raise ValueError("source next stage mismatch")
    authorization = source.get("authorization", {})
    if authorization.get("gpu_used") is not False:
        raise ValueError("source unexpectedly used a GPU")
    if authorization.get("grassmann_fitted") is not False:
        raise ValueError("source unexpectedly fitted Grassmann")


def search_grid_rows() -> list[tuple[str, str, str, str]]:
    labels = ("parity", "majority_threshold", "weighted_threshold", "noisy_threshold")
    distractors = ("none", "moderate", "high")
    budgets = ("fixed_fair_low", "fixed_fair_high")
    return [
        (label, "random_positions_in_token_marker", distractor, budget)
        for label in labels
        for distractor in distractors
        for budget in budgets
    ]


def write_manifest(output_dir: Path, names: list[str]) -> None:
    lines = [f"{sha256(output_dir / name)}  ./{name}" for name in names]
    (output_dir / "LONG_RANGE_DIFFICULTY_REPLAN_READINESS_MANIFEST.v7.7.7.sha256").write_text(
        "\n".join(lines) + "\n", encoding="utf-8", newline="\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--task-difficulty-execution", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    if not args.task_difficulty_execution.is_file():
        raise FileNotFoundError(args.task_difficulty_execution)
    source = json.load(args.task_difficulty_execution.open(encoding="utf-8"))
    validate_source(source)
    if args.output_dir.exists():
        raise FileExistsError(f"refuse overwrite: {args.output_dir}")
    args.output_dir.mkdir(parents=True)

    grid = search_grid_rows()
    grid_name = "LONG_RANGE_DIFFICULTY_SEARCH_GRID.v7.7.7.tsv"
    header = "label_function\tretrieval_geometry\tdistractor_level\tbaseline_budget\texecution_authorized"
    rows = [header] + [f"{a}\t{b}\t{c}\t{d}\tfalse" for a, b, c, d in grid]
    (args.output_dir / grid_name).write_text("\n".join(rows) + "\n", encoding="utf-8", newline="\n")

    readiness = {
        "schema_version": "1.0",
        "protocol_version": "v7.7.7",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "LONG_RANGE_DIFFICULTY_REPLAN_CONTRACT_SIGNED_NO_LAUNCH",
        "audit_role": "PROSPECTIVE_TASK_DIFFICULTY_REPLAN_CONTRACT_ONLY",
        "source_task_difficulty_execution": {
            "file": str(args.task_difficulty_execution.resolve()),
            "sha256": sha256(args.task_difficulty_execution),
            "status": source["status"],
            "selected_k": source.get("selected_k"),
        },
        "replan_requirements": {
            "baseline_handed_source_positions": False,
            "bit_carried_in_marker_token": True,
            "baseline_fairly_and_adequately_trained": True,
            "undertrained_baseline_headroom_forbidden": True,
        },
        "corridor": {
            "lr0_local_equivalence_margin_nats": 0.010,
            "lr1_global_min_improvement_nats": 0.010,
            "headroom_per_truth_seed_min_nats": 0.020,
            "headroom_lcb95_min_nats": 0.010,
            "target_shuffled_margin_nats": 0.010,
            "seed_stable_required": True,
        },
        "predeclared_search_grid": {
            "rows": len(grid),
            "label_function_family": ["parity", "majority_threshold", "weighted_threshold", "noisy_threshold"],
            "retrieval_geometry": ["random_positions_in_token_marker"],
            "distractor_levels": ["none", "moderate", "high"],
            "baseline_budget_levels": ["fixed_fair_low", "fixed_fair_high"],
            "selection_rule": "SMALLEST_DIFFICULTY_SEED_STABLE_CORRIDOR_ELIGIBLE_CONFIG",
            "selection_inputs": ["conventional_baseline", "oracle"],
            "grassmann_consulted_in_selection": False,
        },
        "feasibility_findings_summary": {
            "positions_handed_parity": "BASELINE_SOLVES_HEADROOM_ZERO",
            "random_positions_majority": "BASELINE_SOLVES_HEADROOM_ZERO",
            "random_positions_parity_small_k": "PARTIAL_BUT_SEED_BIMODAL_SHUFFLED_FAILS",
            "random_positions_parity_large_k": "BASELINE_CANNOT_BEAT_LOCAL_LR1_FAILS",
            "structural_conclusion": "FAIR_BASELINE_IS_BIMODAL_CORRIDOR_NARROW",
        },
        "route_closure_branch": {
            "trigger": "NO_FAIR_SEED_STABLE_CORRIDOR_CONFIG_IN_PREDECLARED_GRID",
            "action": "CLOSE_V7_GRASSMANN_PRIMARY_ROUTE_NO_TASK_EXPANSION",
        },
        "authorization": {
            "gpu_authorized": False,
            "grassmann_fitted": False,
            "pilot_execution_authorized": False,
            "architecture_decision_permitted": False,
            "a3_scaling_authorized": False,
            "biological_long_range_claim_permitted": False,
        },
        "next_authorized_stage": "IMPLEMENT_V7_7_8_LONG_RANGE_DIFFICULTY_REPLAN_HARNESS_NO_LAUNCH",
    }
    readiness_name = "LONG_RANGE_DIFFICULTY_REPLAN_READINESS.v7.7.7.json"
    (args.output_dir / readiness_name).write_text(
        json.dumps(readiness, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n"
    )
    write_manifest(args.output_dir, [readiness_name, grid_name])
    print(json.dumps({
        "status": readiness["status"],
        "output_dir": str(args.output_dir.resolve()),
        "search_grid_rows": len(grid),
        "gpu_authorized": False,
        "next_authorized_stage": readiness["next_authorized_stage"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
