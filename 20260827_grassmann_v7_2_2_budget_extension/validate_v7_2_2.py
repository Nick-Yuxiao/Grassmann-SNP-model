from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
REQUIRED = [
    "PROTOCOL_ADDENDUM.v7.2.2.md", "BUDGET_EXTENSION_SPEC.v7.2.2.yaml",
    "DECISIONS.v7.2.2.tsv", "BUDGET_EXTENSION_SCHEDULE.v7.2.2.tsv",
    "validate_v7_2_2.py", "p0/train_budget_extension_v7_2_2.py",
    "p0/run_budget_extension_gpu_worker_v7_2_2.py",
    "p0/assess_budget_extension_v7_2_2.py",
    "p0/run_budget_extension_v7_2_2_nonintrusive.sh",
    "p0/tests/test_v7_2_2.py", "server_ops/SERVER_STEPS.v7.2.2.md",
]


def main() -> None:
    missing = [name for name in REQUIRED if not (ROOT / name).is_file()]
    if missing:
        raise SystemExit(f"missing files: {missing}")
    with (ROOT / "BUDGET_EXTENSION_SCHEDULE.v7.2.2.tsv").open(
        encoding="utf-8-sig", newline=""
    ) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    gpus = {int(row["preferred_physical_gpu"]) for row in rows}
    blocks = {gpu: [row for row in rows if int(row["preferred_physical_gpu"]) == gpu] for gpu in gpus}
    checks = {
        "run_count": len(rows) == 12 and len({row["extension_run_id"] for row in rows}) == 12,
        "source_count": len({row["source_run_id"] for row in rows}) == 12,
        "six_cells": len({(row["model"], row["mask"]) for row in rows}) == 6,
        "gpus": gpus == {1, 3, 4, 5, 6, 7},
        "gpu0_forbidden": 0 not in gpus,
        "gpu2_forbidden": 2 not in gpus,
        "paired_lr_per_gpu": all(
            len(block) == 2
            and {float(row["learning_rate"]) for row in block} == {0.0001, 0.0004}
            and len({(row["model"], row["mask"]) for row in block}) == 1
            and {int(row["order_on_gpu"]) for row in block} == {1, 2}
            for block in blocks.values()
        ),
        "seeds": {row["mask_seed"] for row in rows} == {"92001"}
        and {row["init_seed"] for row in rows} == {"82001"},
        "source_step": {row["source_step"] for row in rows} == {"20000"},
        "target_step": {row["target_step"] for row in rows} == {"30000"},
    }
    failed = [name for name, passed in checks.items() if not passed]
    if failed:
        raise SystemExit(f"validation failed: {failed}")
    print(json.dumps({
        "status": "PASS", "protocol_version": "v7.2.2",
        "manifest_entries": len(REQUIRED), "runs": 12,
        "additional_steps": 120000, "source_step": 20000, "target_step": 30000,
        "allowed_physical_gpus": sorted(gpus), "gpu0": "FORBIDDEN",
        "gpu2": "FORBIDDEN", "selective_extension": "FORBIDDEN",
        "decision_holdout": "FORBIDDEN", "hgdp_access": "FORBIDDEN",
        "formal_a1r": "BLOCKED",
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
