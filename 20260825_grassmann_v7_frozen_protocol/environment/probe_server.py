from __future__ import annotations

import json
import platform
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def run(*command: str) -> dict[str, object]:
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    return {
        "command": list(command),
        "exit_code": completed.returncode,
        "stdout": completed.stdout.strip(),
        "stderr": completed.stderr.strip(),
    }


def main() -> None:
    root = Path(__file__).resolve().parent
    payload = {
        "schema_version": "1.0",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "hostname": platform.node(),
        "python": platform.python_version(),
        "nvidia_smi": run(
            "nvidia-smi",
            "--query-gpu=index,name,memory.total,memory.free,driver_version",
            "--format=csv,noheader,nounits",
        ),
        "disk": run("df", "-B1", str(root)),
        "slurm": {
            "sbatch": shutil.which("sbatch"),
            "sinfo": run("sinfo", "-h", "-o", "%P|%G|%D|%a") if shutil.which("sinfo") else None,
        },
    }
    output = root / "SERVER_RESOURCE.json"
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
