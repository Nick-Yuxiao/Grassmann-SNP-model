from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--schedule", type=Path, required=True)
    parser.add_argument("--gpu", type=int, required=True)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--firewall-dir", type=Path, required=True)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--trainer", type=Path, required=True)
    args = parser.parse_args()
    with args.schedule.open(encoding="utf-8-sig", newline="") as handle:
        rows = [
            row for row in csv.DictReader(handle, delimiter="\t")
            if int(row["preferred_physical_gpu"]) == args.gpu
        ]
    if len(rows) != 1 or int(rows[0]["order_on_gpu"]) != 1:
        raise SystemExit(f"GPU {args.gpu} expected exactly one frozen row, found {len(rows)}")
    row = rows[0]
    if float(row["learning_rate"]) != 0.0004:
        raise SystemExit("non-selected LR row is forbidden")
    output = args.output_root / row["final_run_id"]
    log = args.output_root / f"{row['final_run_id']}.log"
    source = args.source_root / row["source_run_id"]
    command = [
        sys.executable, str(args.trainer), "--model", row["model"],
        "--mask", row["mask"], "--mask-seed", row["mask_seed"],
        "--init-seed", row["init_seed"], "--learning-rate", row["learning_rate"],
        "--source-step", row["source_step"], "--target-step", row["target_step"],
        "--data-dir", str(args.data_dir), "--firewall-dir", str(args.firewall_dir),
        "--source-run-dir", str(source), "--output-dir", str(output),
    ]
    with log.open("wb") as handle:
        code = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT).returncode
    if code != 0:
        raise SystemExit(f"final diagnostic failed: {row['final_run_id']} exit={code}")


if __name__ == "__main__":
    main()

