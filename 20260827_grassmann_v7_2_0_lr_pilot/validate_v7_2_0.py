from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
REQUIRED = [
    "PROTOCOL.v7.2.0.md", "LR_PILOT_SPEC.v7.2.0.yaml", "DECISIONS.v7.2.0.tsv",
    "LR_PILOT_SCHEDULE.v7.2.0.tsv", "validate_v7_2_0.py",
    "p0/freeze_lr_validation_firewall_v7_2_0.py", "p0/train_lr_pilot_v7_2_0.py",
    "p0/run_lr_pilot_gpu_worker_v7_2_0.py", "p0/select_shared_lr_v7_2_0.py",
    "p0/run_lr_pilot_v7_2_0_nonintrusive.sh", "p0/tests/test_v7_2_0.py",
    "server_ops/SERVER_STEPS.v7.2.0.md",
]


def main() -> None:
    missing = [name for name in REQUIRED if not (ROOT / name).is_file()]
    if missing:
        raise SystemExit(f"missing files: {missing}")
    with (ROOT / "LR_PILOT_SCHEDULE.v7.2.0.tsv").open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    models = {row["model"] for row in rows}
    masks = {row["mask"] for row in rows}
    lrs = {float(row["peak_lr"]) for row in rows}
    gpus = {int(row["preferred_physical_gpu"]) for row in rows}
    cells = {(float(row["peak_lr"]), row["model"], row["mask"]) for row in rows}
    gpu_blocks = {
        gpu: [row for row in rows if int(row["preferred_physical_gpu"]) == gpu]
        for gpu in gpus
    }
    checks = {
        "run_count": len(rows) == 24,
        "unique_runs": len({row["run_id"] for row in rows}) == 24,
        "complete_cells": len(cells) == 24,
        "models": len(models) == 3,
        "masks": len(masks) == 2,
        "lrs": lrs == {0.0001, 0.0002, 0.0004, 0.0008},
        "gpus": gpus == {1, 3, 4, 5, 6, 7},
        "gpu0_forbidden": 0 not in gpus,
        "gpu2_excluded": 2 not in gpus,
        "frozen_seed": {row["mask_seed"] for row in rows} == {"92001"} and {row["init_seed"] for row in rows} == {"82001"},
        "steps": {row["steps"] for row in rows} == {"4000"},
        "warmup": {row["warmup_steps"] for row in rows} == {"500"},
        "four_lrs_per_gpu": all(
            len(block) == 4 and {float(row["peak_lr"]) for row in block} == lrs
            for block in gpu_blocks.values()
        ),
        "one_model_mask_cell_per_gpu": all(
            len({(row["model"], row["mask"]) for row in block}) == 1
            for block in gpu_blocks.values()
        ),
        "gpu_order_permutation": all(
            {int(row["order_on_gpu"]) for row in block} == {1, 2, 3, 4}
            for block in gpu_blocks.values()
        ),
    }
    failed = [name for name, passed in checks.items() if not passed]
    if failed:
        raise SystemExit(f"validation failed: {failed}")
    print(json.dumps({
        "status": "PASS", "protocol_version": "v7.2.0", "manifest_entries": len(REQUIRED),
        "runs": len(rows), "run_steps": sum(int(row["steps"]) for row in rows),
        "allowed_physical_gpus": sorted(gpus), "gpu0": "FORBIDDEN", "gpu2": "EXCLUDED",
        "hgdp_access": "FORBIDDEN", "formal_a1r": "BLOCKED",
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
