from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


EXPECTED_CONSTANTS = {
    "delta_min": 0.010,
    "delta_NI": 0.010,
    "delta_LD": 0.010,
    "overfit_thr": 0.010,
    "pc_control_thr": 0.005,
}


def load_json(path: Path) -> dict[str, Any] | None:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (FileNotFoundError, json.JSONDecodeError, OSError):
        return None


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate a machine P0 readiness verdict without converting missing evidence into NO_GO.")
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    checks: dict[str, dict[str, Any]] = {}

    smoke = load_json(root / "environment" / "ENV_SMOKE.json")
    resources = load_json(root / "environment" / "SERVER_RESOURCE.json")
    audit = load_json(root / "p0" / "REMOTE_SAFETY_AUDIT.json")
    checks["T00"] = {
        "pass": bool(
            smoke
            and smoke.get("status") == "PASS"
            and resources
            and audit
            and audit.get("status") in {"SAFE_IDLE_GPU_AVAILABLE", "BUSY_OR_GPU_UNAVAILABLE"}
            and (root / "environment" / "requirements-cu128.lock").is_file()
        ),
        "note": "Environment smoke, server resource record, non-interference audit and lockfile are required. A busy audit is valid evidence but no new GPU job may start until a fresh idle audit passes.",
    }

    panel = load_json(root / "p0" / "runtime_t01" / "PANEL_MANIFEST.json")
    checks["T01"] = {
        "pass": bool(panel and panel.get("status") == "SIGNED_INPUTS_PRESENT" and (root / "p0" / "runtime_t01" / "DATA_BRANCH_DECISION.md").is_file()),
        "note": "A signed branch decision and hashed panel inputs are required.",
    }

    decision_path = root / "DECISIONS.v7.0.1.tsv"
    metric_path = root / "METRIC_DEFINITIONS.md"
    decisions_text = decision_path.read_text(encoding="utf-8") if decision_path.is_file() else ""
    metric_text = metric_path.read_text(encoding="utf-8") if metric_path.is_file() else ""
    constants_present = all(name in decisions_text and name in metric_text and str(value) in decisions_text for name, value in EXPECTED_CONSTANTS.items())
    checks["T02"] = {
        "pass": constants_present,
        "note": "Every independently named frozen constant must appear in both the decision ledger and metric dictionary.",
    }

    profile = load_json(root / "p0" / "profile_runtime" / "PROFILE_REPORT.json")
    exact_cells = {
        (row.get("model"), row.get("sequence_length"))
        for row in (profile or {}).get("profiles", [])
        if row.get("status") == "PASS"
    }
    expected_cells = {
        (model, length)
        for model in ("local_attn", "local_attn_gpc", "grassmann_full")
        for length in (8192, 131072, 262144)
    }
    checks["T03"] = {
        "pass": bool(profile and profile.get("valid_t03_measurement") and expected_cells <= exact_cells),
        "note": "A 100-step CUDA profile at all exact P0/A1 planning lengths is required; CPU dry-runs are never accepted.",
    }

    contract = load_json(root / "p0" / "runtime_t04" / "COMPUTE_CONTRACT.json")
    checks["T04"] = {
        "pass": bool(contract and contract.get("status") == "P0_CAPACITY_GO"),
        "note": "The signed capacity must cover the measured A1 plan at <=80%, including the 2x engineering margin.",
    }
    blockers = [task for task, result in checks.items() if not result["pass"]]
    status = "GO_TO_A0" if not blockers else "NOT_READY"
    payload = {
        "schema_version": "1.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "scope": "P0 readiness only; not an A1 scientific GO/NO-GO",
        "checks": checks,
        "blockers": blockers,
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    json_path = args.output_dir / "P0_VERDICT.json"
    md_path = args.output_dir / "P0_VERDICT.md"
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    md_path.write_text(
        "# V7 P0 readiness verdict\n\n"
        f"- Status: **{status}**\n"
        "- Scope: readiness to enter A0 only; not the A1 scientific verdict.\n"
        f"- Blocking tasks: {', '.join(blockers) if blockers else 'none'}\n\n"
        + "\n".join(f"- {task}: {'PASS' if result['pass'] else 'BLOCKED'} — {result['note']}" for task, result in checks.items())
        + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"status": status, "blockers": blockers, "output": str(json_path)}, indent=2))


if __name__ == "__main__":
    main()
