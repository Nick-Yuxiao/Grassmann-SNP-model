from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def pilot_rows() -> list[tuple[int, int, str, str, str]]:
    cells = [
        ("R0_LOCAL", "absent", "absent"),
        ("R1_ROUTER", "present", "absent"),
        ("R2_GRASSMANN_ONLY", "absent", "present"),
        ("R3_ROUTER_GRASSMANN", "present", "present"),
    ]
    return [
        (truth_seed, init_seed, cell, router, grassmann)
        for truth_seed in range(77301, 77307)
        for init_seed in (87401, 87402)
        for cell, router, grassmann in cells
    ]

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--task-validity", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise FileExistsError(f"refuse overwrite: {args.output_dir}")
    validity = json.load(args.task_validity.open(encoding="utf-8"))
    if validity.get("status") != "LONG_RANGE_TASK_VALIDITY_PASS":
        raise ValueError("v7.7.2 task validity did not pass")
    if validity.get("authorization", {}).get("gpu_used") is not False:
        raise ValueError("v7.7.2 source must be CPU-only")
    gates = validity.get("gates", {})
    if gates.get("all_pass") is not True:
        raise ValueError("all v7.7.2 gates must pass")
    args.output_dir.mkdir(parents=True)
    schedule = args.output_dir / "LONG_RANGE_PILOT_DRAFT.v7.7.3.tsv"
    rows = ["truth_seed\tinit_seed\tcell\tglobal_router\tgrassmann_component\texecution_authorized"]
    for truth_seed, init_seed, cell, router, grassmann in pilot_rows():
        rows.append(f"{truth_seed}\t{init_seed}\t{cell}\t{router}\t{grassmann}\tfalse")
    schedule.write_text("\n".join(rows) + "\n", encoding="utf-8", newline="\n")
    readiness = {
        "schema_version": "1.0",
        "protocol_version": "v7.7.3",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "LONG_RANGE_ARCHITECTURE_PILOT_CONTRACT_SIGNED_IMPLEMENTATION_ONLY",
        "source_task_validity": {"file": str(args.task_validity.resolve()), "sha256": sha256(args.task_validity)},
        "effective_n": 6,
        "pilot_execution_cells": 48,
        "gpu_authorized": False,
        "pilot_execution_authorized": False,
        "architecture_decision_permitted": False,
        "pilot_arm_means_released": False,
        "next_authorized_stage": "IMPLEMENT_V7_7_4_LONG_RANGE_PILOT_HARNESS_NO_LAUNCH",
    }
    report = args.output_dir / "LONG_RANGE_ARCHITECTURE_PILOT_READINESS.v7.7.3.json"
    report.write_text(json.dumps(readiness, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    manifest = args.output_dir / "LONG_RANGE_ARCHITECTURE_PILOT_READINESS_MANIFEST.v7.7.3.sha256"
    manifest.write_text(
        f"{sha256(report)}  ./{report.name}\n{sha256(schedule)}  ./{schedule.name}\n",
        encoding="utf-8", newline="\n",
    )
    print(json.dumps({"status": readiness["status"], "output_dir": str(args.output_dir.resolve()),
                      "next_authorized_stage": readiness["next_authorized_stage"]}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
