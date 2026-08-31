from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


MODELS = ("local_attn", "local_attn_gpc", "grassmann_full")
RUNS_PER_MODEL_AND_LENGTH = {
    131072: 6 * 10,  # two diagnostic cells + four confirmatory fairness cells, 5x2 repeats
    262144: 2 * 10,  # one cross-chrom mask, two fairness cells, 5x2 repeats
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description="Build the T04 capacity contract from exact T03 measurements.")
    parser.add_argument("--profile-report", type=Path, required=True)
    parser.add_argument("--train-steps-per-run", type=int, required=True)
    parser.add_argument("--training-budget-reference", required=True, help="Signed record that freezes the run step/epoch budget.")
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
        raise SystemExit("concurrency limit cannot exceed GPU count")

    payload: dict[str, Any] = json.loads(args.profile_report.read_text(encoding="utf-8"))
    if not payload.get("valid_t03_measurement"):
        raise SystemExit("PROFILE_REPORT is not a valid 100-step CUDA T03 measurement")
    lookup = {
        (row.get("model"), row.get("sequence_length")): row
        for row in payload.get("profiles", [])
        if row.get("status") == "PASS"
    }
    missing = [
        (model, length)
        for model in MODELS
        for length in RUNS_PER_MODEL_AND_LENGTH
        if (model, length) not in lookup
    ]
    if missing:
        raise SystemExit("Missing exact T03 profile cells: " + repr(missing))

    cells = []
    raw_gpu_hours = 0.0
    storage_bytes = 0
    for model in MODELS:
        for length, runs in RUNS_PER_MODEL_AND_LENGTH.items():
            row = lookup[(model, length)]
            per_run_hours = float(row["seconds_per_step"]) * args.train_steps_per_run / 3600.0
            cell_hours = per_run_hours * runs
            raw_gpu_hours += cell_hours
            storage_bytes += int(row["checkpoint_bytes"]) * runs
            cells.append({
                "model": model,
                "sequence_length": length,
                "runs": runs,
                "seconds_per_step_measured": row["seconds_per_step"],
                "train_steps_per_run": args.train_steps_per_run,
                "gpu_hours": cell_hours,
            })

    engineering_margin = 2.0
    required_gpu_hours = raw_gpu_hours * engineering_margin
    required_storage_gib = (storage_bytes / (2**30)) * engineering_margin
    gpu_capacity_pass = required_gpu_hours <= 0.8 * args.signed_gpu_hours
    storage_capacity_pass = required_storage_gib <= 0.8 * args.signed_storage_gib
    scheduled_wall_hours = required_gpu_hours / args.concurrency_limit
    window_pass = scheduled_wall_hours <= 0.8 * args.available_window_hours
    status = "P0_CAPACITY_GO" if gpu_capacity_pass and storage_capacity_pass and window_pass else "P0_CAPACITY_NO_GO"

    contract = {
        "schema_version": "1.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "scope": "P0 capacity only; not an A1 scientific verdict",
        "profile_report_sha256": sha256(args.profile_report),
        "training_budget_reference": args.training_budget_reference,
        "train_steps_per_run": args.train_steps_per_run,
        "engineering_margin_multiplier": engineering_margin,
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
            "scheduled_wall_hours_with_margin": scheduled_wall_hours,
        },
        "checks": {
            "gpu_hours_within_80pct": gpu_capacity_pass,
            "storage_within_80pct": storage_capacity_pass,
            "wall_window_within_80pct": window_pass,
        },
        "a1_cells": cells,
        "a2_funding": "UNFUNDED",
        "a3_funding": "UNFUNDED",
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    json_path = args.output_dir / "COMPUTE_CONTRACT.json"
    md_path = args.output_dir / "COMPUTE_CONTRACT.md"
    json_path.write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    md_path.write_text(
        "# V7 T04 compute contract\n\n"
        f"- Status: **{status}**\n"
        "- Scope: P0 capacity only; this is not an A1 scientific GO/NO-GO.\n"
        f"- Required GPU-hours (2x margin): {required_gpu_hours:.3f}\n"
        f"- Required storage GiB (2x margin): {required_storage_gib:.3f}\n"
        f"- Scheduled wall-hours at concurrency {args.concurrency_limit}: {scheduled_wall_hours:.3f}\n"
        f"- Signed GPU-hours / storage GiB / window h: {args.signed_gpu_hours:.3f} / {args.signed_storage_gib:.3f} / {args.available_window_hours:.3f}\n"
        f"- Training budget reference: `{args.training_budget_reference}`\n"
        f"- Profile SHA-256: `{contract['profile_report_sha256']}`\n"
        "- A2/A3: UNFUNDED\n",
        encoding="utf-8",
    )
    print(json.dumps({"status": status, "required_gpu_hours": required_gpu_hours, "required_storage_gib": required_storage_gib}, indent=2))
    if status != "P0_CAPACITY_GO":
        raise SystemExit(4)


if __name__ == "__main__":
    main()
