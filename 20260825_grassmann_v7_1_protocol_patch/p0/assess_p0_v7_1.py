from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


def load(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}


def main() -> None:
    parser = argparse.ArgumentParser(description="Issue v7.1.0 readiness verdict; never a scientific GO")
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--t00-smoke", type=Path, required=True)
    parser.add_argument("--t00-resource", type=Path, required=True)
    parser.add_argument("--t00-runtime-manifest", type=Path, required=True)
    parser.add_argument("--panel-manifest", type=Path, required=True)
    parser.add_argument("--profile-report", type=Path, required=True)
    parser.add_argument("--compute-contract", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    smoke = load(args.t00_smoke)
    resource = load(args.t00_resource)
    panel = load(args.panel_manifest)
    profile = load(args.profile_report)
    contract = load(args.compute_contract)
    decisions = root / "DECISIONS.v7.1.0.tsv"
    metrics = root / "METRIC_DEFINITIONS.v7.1.0.md"
    checks = {
        "T00": bool(
            smoke.get("status") == "PASS"
            and resource.get("nvidia_smi", {}).get("exit_code") == 0
            and args.t00_runtime_manifest.is_file()
        ),
        "T01": bool(panel.get("status") == "SIGNED_INPUTS_PRESENT" and panel.get("all_partitions_disjoint")),
        "T02": bool(decisions.is_file() and metrics.is_file()),
        "T03": bool(profile.get("valid_t03_measurement") and profile.get("architecture", {}).get("attention_window") == 256),
        "T04": bool(contract.get("status") == "P0_CAPACITY_GO"),
    }
    blockers = [task for task, passed in checks.items() if not passed]
    status = "GO_TO_A0" if not blockers else "NOT_READY"
    payload = {
        "schema_version": "1.0",
        "protocol_version": "v7.1.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "scope": "P0 readiness only; not A1 scientific GO/NO-GO",
        "checks": checks,
        "blockers": blockers,
    }
    args.output_dir.mkdir(parents=True, exist_ok=False)
    output = args.output_dir / "P0_VERDICT.v7.1.0.json"
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
