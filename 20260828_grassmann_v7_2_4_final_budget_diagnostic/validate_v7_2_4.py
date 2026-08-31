from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
REQUIRED = [
    "PROTOCOL_ADDENDUM.v7.2.4.md", "FINAL_BUDGET_SPEC.v7.2.4.yaml",
    "DECISIONS.v7.2.4.tsv", "FINAL_BUDGET_SCHEDULE.v7.2.4.tsv",
    "validate_v7_2_4.py", "p0/train_final_budget_v7_2_4.py",
    "p0/run_final_budget_gpu_worker_v7_2_4.py",
    "p0/assess_final_budget_v7_2_4.py",
    "p0/run_final_budget_v7_2_4_nonintrusive.sh",
    "p0/tests/test_v7_2_4.py", "server_ops/SERVER_STEPS.v7.2.4.md",
]


def main() -> None:
    missing = [name for name in REQUIRED if not (ROOT / name).is_file()]
    if missing:
        raise SystemExit(f"missing files: {missing}")
    with (ROOT / "FINAL_BUDGET_SCHEDULE.v7.2.4.tsv").open(
        encoding="utf-8-sig", newline=""
    ) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    gpus = {int(row["preferred_physical_gpu"]) for row in rows}
    cells = {(row["model"], row["mask"]) for row in rows}
    protocol = (ROOT / "PROTOCOL_ADDENDUM.v7.2.4.md").read_text(encoding="utf-8")
    checks = {
        "run_count": len(rows) == 6 and len({row["final_run_id"] for row in rows}) == 6,
        "source_count": len({row["source_run_id"] for row in rows}) == 6,
        "six_cells": len(cells) == 6,
        "three_models": len({row["model"] for row in rows}) == 3,
        "two_masks": len({row["mask"] for row in rows}) == 2,
        "gpus": gpus == {1, 3, 4, 5, 6, 7},
        "one_cell_per_gpu": all(
            sum(int(row["preferred_physical_gpu"]) == gpu for row in rows) == 1
            for gpu in gpus
        ),
        "selected_lr_only": {row["learning_rate"] for row in rows} == {"0.0004"},
        "seeds": {row["mask_seed"] for row in rows} == {"92001"}
        and {row["init_seed"] for row in rows} == {"82001"},
        "source_step": {row["source_step"] for row in rows} == {"30000"},
        "target_step": {row["target_step"] for row in rows} == {"40000"},
        "hard_stop": "No run may continue beyond 40k" in protocol,
        "architecture_block": "architecture GO/NO-GO decision" in protocol,
        "complete_factorial": "complete `3 models x 2 masks` factorial" in protocol,
    }
    failed = [name for name, passed in checks.items() if not passed]
    if failed:
        raise SystemExit(f"validation failed: {failed}")
    print(json.dumps({
        "status": "PASS", "protocol_version": "v7.2.4",
        "manifest_entries": len(REQUIRED), "runs": 6,
        "additional_steps": 60000, "source_step": 30000, "target_step": 40000,
        "learning_rate": 0.0004, "allowed_physical_gpus": sorted(gpus),
        "gpu0": "FORBIDDEN", "gpu2": "FORBIDDEN",
        "failed_cell_selection": "FORBIDDEN",
        "extension_beyond_40k": "FORBIDDEN",
        "decision_holdout": "FORBIDDEN", "hgdp_access": "FORBIDDEN",
        "formal_a1r": "BLOCKED",
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

