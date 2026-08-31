from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
REQUIRED = [
    "PROTOCOL_ADDENDUM.v7.2.1.md", "BUDGET_BRIDGE_SPEC.v7.2.1.yaml",
    "DECISIONS.v7.2.1.tsv", "BUDGET_BRIDGE_SCHEDULE.v7.2.1.tsv",
    "validate_v7_2_1.py", "p0/train_budget_bridge_v7_2_1.py",
    "p0/run_budget_bridge_gpu_worker_v7_2_1.py", "p0/assess_budget_bridge_v7_2_1.py",
    "p0/run_budget_bridge_v7_2_1_nonintrusive.sh", "p0/tests/test_v7_2_1.py",
    "server_ops/SERVER_STEPS.v7.2.1.md",
]


def main() -> None:
    missing = [name for name in REQUIRED if not (ROOT / name).is_file()]
    if missing:
        raise SystemExit(f"missing files: {missing}")
    with (ROOT / "BUDGET_BRIDGE_SCHEDULE.v7.2.1.tsv").open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    gpus = {int(row["preferred_physical_gpu"]) for row in rows}
    gpu_blocks = {gpu: [row for row in rows if int(row["preferred_physical_gpu"]) == gpu] for gpu in gpus}
    cells = {(row["model"], row["mask"]) for row in rows}
    checks = {
        "run_count": len(rows) == 12 and len({row["run_id"] for row in rows}) == 12,
        "six_cells": len(cells) == 6,
        "gpus": gpus == {1,3,4,5,6,7},
        "gpu0_forbidden": 0 not in gpus,
        "gpu2_excluded": 2 not in gpus,
        "paired_lr_per_gpu": all(
            len(block) == 2 and {float(row["learning_rate"]) for row in block} == {0.0001,0.0004}
            and len({(row["model"],row["mask"]) for row in block}) == 1
            and {int(row["order_on_gpu"]) for row in block} == {1,2}
            for block in gpu_blocks.values()
        ),
        "seeds": {row["mask_seed"] for row in rows} == {"92001"} and {row["init_seed"] for row in rows} == {"82001"},
        "steps": {row["steps"] for row in rows} == {"20000"},
        "warmup": {row["warmup_steps"] for row in rows} == {"500"},
        "counterbalanced": sum(
            1 for block in gpu_blocks.values()
            if next(row for row in block if int(row["order_on_gpu"]) == 1)["learning_rate"] == "0.0001"
        ) == 3,
    }
    failed = [name for name, passed in checks.items() if not passed]
    if failed:
        raise SystemExit(f"validation failed: {failed}")
    print(json.dumps({
        "status":"PASS", "protocol_version":"v7.2.1", "manifest_entries":len(REQUIRED),
        "runs":12, "run_steps":240000, "allowed_physical_gpus":sorted(gpus),
        "gpu0":"FORBIDDEN", "gpu2":"EXCLUDED", "decision_holdout":"FORBIDDEN",
        "hgdp_access":"FORBIDDEN", "formal_a1r":"BLOCKED",
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

