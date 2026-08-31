from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


MODELS = (
    "local_attn_8m_w256",
    "local_attn_gpc_8m_w256",
    "grassmann_full_8m_w256",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description="Build v7.1.0 T04 capacity contract from exact T03 cells")
    parser.add_argument("--grid", type=Path, required=True)
    parser.add_argument("--panel-manifest", type=Path, required=True)
    parser.add_argument("--profile-report", type=Path, required=True)
    parser.add_argument("--train-steps-per-run", type=int, required=True)
    parser.add_argument("--training-budget-reference", required=True)
    parser.add_argument("--signed-gpu-hours", type=float, required=True)
    parser.add_argument("--signed-storage-gib", type=float, required=True)
    parser.add_argument("--gpu-count", type=int, required=True)
    parser.add_argument("--concurrency-limit", type=int, required=True)
    parser.add_argument("--available-window-hours", type=float, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if min(args.train_steps_per_run, args.gpu_count, args.concurrency_limit) <= 0:
        raise SystemExit("steps, GPU count and concurrency must be positive")
    if args.concurrency_limit > args.gpu_count:
        raise SystemExit("concurrency cannot exceed signed GPU count")

    grid = json.loads(args.grid.read_text(encoding="utf-8"))
    panel = json.loads(args.panel_manifest.read_text(encoding="utf-8"))
    profile = json.loads(args.profile_report.read_text(encoding="utf-8"))
    if grid.get("protocol_version") != "v7.1.0" or panel.get("protocol_version") != "v7.1.0":
        raise SystemExit("grid and panel manifest must both be v7.1.0")
    if panel.get("status") != "SIGNED_INPUTS_PRESENT" or not panel.get("all_partitions_disjoint"):
        raise SystemExit("panel manifest is not signed/disjoint")
    if not profile.get("valid_t03_measurement"):
        raise SystemExit("profile is not a valid v7.1.0 100-step CUDA measurement")
    if profile.get("architecture", {}).get("attention_window") != 256:
        raise SystemExit("profile did not use frozen attention_window=256")
    visible = str(profile.get("cuda_visible_devices", "")).strip()
    if visible in {"", "0"}:
        raise SystemExit("profile used forbidden or unrecorded physical GPU 0")

    length = int(panel["chr22"]["variant_count"])
    lookup = {
        (row.get("model"), row.get("sequence_length")): row
        for row in profile.get("profiles", [])
        if row.get("status") == "PASS"
    }
    missing = [(model, length) for model in MODELS if (model, length) not in lookup]
    if missing:
        raise SystemExit(f"missing exact data-derived T03 profile cells: {missing}")

    planning = grid["T04_planning"]
    cells_per_model = sum(
        int(planning[key])
        for key in (
            "diagnostic_cells_per_model",
            "primary_confirmatory_cells_per_model",
            "stress_sensitivity_cells_per_model",
        )
    )
    runs_per_model = cells_per_model * int(planning["runs_per_cell"])
    raw_gpu_hours = 0.0
    raw_storage_bytes = 0
    cells = []
    for model in MODELS:
        row = lookup[(model, length)]
        per_run_hours = float(row["seconds_per_step"]) * args.train_steps_per_run / 3600.0
        gpu_hours = per_run_hours * runs_per_model
        raw_gpu_hours += gpu_hours
        raw_storage_bytes += int(row["checkpoint_bytes"]) * runs_per_model
        cells.append({
            "model": model,
            "sequence_length": length,
            "planned_cells": cells_per_model,
            "runs_per_cell": int(planning["runs_per_cell"]),
            "runs": runs_per_model,
            "seconds_per_step_measured": row["seconds_per_step"],
            "train_steps_per_run": args.train_steps_per_run,
            "gpu_hours": gpu_hours,
        })

    margin = float(planning["engineering_margin_multiplier"])
    usage_limit = float(planning["capacity_usage_limit"])
    required_gpu_hours = raw_gpu_hours * margin
    required_storage_gib = raw_storage_bytes / (2**30) * margin
    wall_hours = required_gpu_hours / args.concurrency_limit
    checks = {
        "gpu_hours_within_limit": required_gpu_hours <= usage_limit * args.signed_gpu_hours,
        "storage_within_limit": required_storage_gib <= usage_limit * args.signed_storage_gib,
        "wall_window_within_limit": wall_hours <= usage_limit * args.available_window_hours,
    }
    status = "P0_CAPACITY_GO" if all(checks.values()) else "P0_CAPACITY_NO_GO"
    contract = {
        "schema_version": "1.1",
        "protocol_version": "v7.1.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "scope": "P0 capacity only; not A1 scientific GO/NO-GO",
        "hashes": {"grid": sha256(args.grid), "panel_manifest": sha256(args.panel_manifest), "profile_report": sha256(args.profile_report)},
        "training_budget_reference": args.training_budget_reference,
        "signed_capacity": {
            "gpu_hours": args.signed_gpu_hours,
            "storage_gib": args.signed_storage_gib,
            "gpu_count": args.gpu_count,
            "concurrency_limit": args.concurrency_limit,
            "available_window_hours": args.available_window_hours,
        },
        "required": {
            "raw_gpu_hours": raw_gpu_hours,
            "gpu_hours_with_margin": required_gpu_hours,
            "storage_gib_with_margin": required_storage_gib,
            "scheduled_wall_hours_with_margin": wall_hours,
        },
        "checks": checks,
        "a1_cells": cells,
        "A2_funding": "UNFUNDED",
        "A3_funding": "UNFUNDED",
    }
    args.output_dir.mkdir(parents=True, exist_ok=False)
    output = args.output_dir / "COMPUTE_CONTRACT.v7.1.0.json"
    output.write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": status, "required_gpu_hours": required_gpu_hours, "required_storage_gib": required_storage_gib, "output": str(output)}, indent=2))
    if status != "P0_CAPACITY_GO":
        raise SystemExit(4)


if __name__ == "__main__":
    main()
