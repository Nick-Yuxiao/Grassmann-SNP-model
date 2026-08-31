from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def run(command: list[str], timeout: int = 15) -> dict[str, Any]:
    try:
        completed = subprocess.run(
            command, capture_output=True, text=True, check=False, timeout=timeout
        )
        return {
            "available": True,
            "exit_code": completed.returncode,
            "stdout": completed.stdout.strip(),
            "stderr": completed.stderr.strip(),
        }
    except (FileNotFoundError, subprocess.TimeoutExpired) as exc:
        return {
            "available": False,
            "exit_code": None,
            "stdout": "",
            "stderr": type(exc).__name__,
        }


def csv_rows(text: str, columns: list[str]) -> list[dict[str, str]]:
    if not text.strip():
        return []
    return [dict(zip(columns, row, strict=False)) for row in csv.reader(text.splitlines(), skipinitialspace=True)]


def parse_float(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def main() -> None:
    parser = argparse.ArgumentParser(description="Read-only server/GPU safety audit; never kills or signals processes.")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-idle-memory-mib", type=float, default=1024.0)
    parser.add_argument("--max-idle-utilization-pct", type=float, default=5.0)
    parser.add_argument("--require-idle-gpu", action="store_true")
    args = parser.parse_args()

    gpu_query = run([
        "nvidia-smi",
        "--query-gpu=index,uuid,name,memory.total,memory.used,memory.free,utilization.gpu,temperature.gpu,pstate,driver_version",
        "--format=csv,noheader,nounits",
    ])
    process_query = run([
        "nvidia-smi",
        "--query-compute-apps=gpu_uuid,pid,process_name,used_memory",
        "--format=csv,noheader,nounits",
    ])
    gpu_columns = [
        "index", "uuid", "name", "memory_total_mib", "memory_used_mib",
        "memory_free_mib", "utilization_gpu_pct", "temperature_c", "pstate", "driver_version",
    ]
    proc_columns = ["gpu_uuid", "pid", "process_name", "used_memory_mib"]
    gpus = csv_rows(gpu_query["stdout"], gpu_columns) if gpu_query["exit_code"] == 0 else []
    gpu_processes = csv_rows(process_query["stdout"], proc_columns) if process_query["exit_code"] == 0 else []
    busy_uuids = {row["gpu_uuid"] for row in gpu_processes}

    idle_indices: list[int] = []
    for gpu in gpus:
        used = parse_float(gpu["memory_used_mib"])
        util = parse_float(gpu["utilization_gpu_pct"])
        is_idle = (
            gpu["uuid"] not in busy_uuids
            and used is not None
            and used <= args.max_idle_memory_mib
            and util is not None
            and util <= args.max_idle_utilization_pct
        )
        gpu["idle_by_policy"] = str(is_idle).lower()
        if is_idle:
            idle_indices.append(int(gpu["index"]))

    user = os.environ.get("USER") or os.environ.get("USERNAME") or ""
    ps_query = run(["ps", "-u", user, "-o", "pid=,ppid=,stat=,etimes=,pcpu=,pmem=,comm="]) if user else {
        "available": False, "exit_code": None, "stdout": "", "stderr": "user_unknown"
    }
    scheduler: dict[str, Any] = {"kind": "none", "query": None}
    if shutil.which("squeue"):
        scheduler = {
            "kind": "slurm",
            "query": run(["squeue", "-u", user, "-h", "-o", "%i|%P|%j|%T|%M|%D|%R"]),
        }

    disk_query = run(["df", "-B1", str(args.output.parent.resolve())])
    status = "SAFE_IDLE_GPU_AVAILABLE" if idle_indices else "BUSY_OR_GPU_UNAVAILABLE"
    payload = {
        "schema_version": "1.0",
        "observed_at_utc": datetime.now(timezone.utc).isoformat(),
        "policy": {
            "read_only": True,
            "signals_or_kills_processes": False,
            "max_idle_memory_mib": args.max_idle_memory_mib,
            "max_idle_utilization_pct": args.max_idle_utilization_pct,
        },
        "gpus": gpus,
        "gpu_compute_processes": gpu_processes,
        "idle_gpu_indices": idle_indices,
        "current_user_processes": ps_query,
        "scheduler": scheduler,
        "disk": disk_query,
        "status": status,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))
    if args.require_idle_gpu and not idle_indices:
        raise SystemExit(3)


if __name__ == "__main__":
    main()
