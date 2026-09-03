from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
REQUIRED = (
    "PROTOCOL_ADDENDUM.v7.7.7.md",
    "LONG_RANGE_DIFFICULTY_REPLAN_SPEC.v7.7.7.yaml",
    "DECISIONS.v7.7.7.tsv",
    "validate_v7_7_7.py",
    "p0/build_difficulty_replan_readiness_v7_7_7.py",
    "p0/tests/test_v7_7_7.py",
    "server_ops/SERVER_STEPS.v7.7.7.md",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    missing = [item for item in REQUIRED if not (ROOT / item).is_file()]
    if missing:
        raise ValueError(f"missing files: {missing}")
    spec = (ROOT / "LONG_RANGE_DIFFICULTY_REPLAN_SPEC.v7.7.7.yaml").read_text(encoding="utf-8")
    protocol = (ROOT / "PROTOCOL_ADDENDUM.v7.7.7.md").read_text(encoding="utf-8")
    flat = " ".join(protocol.split()).lower()
    checks = {
        "source_unresolved": "LONG_RANGE_TASK_DIFFICULTY_UNRESOLVED" in spec,
        "no_handed_positions": "baseline_handed_source_positions: false" in spec,
        "no_gpu": "gpu_authorized: false" in spec,
        "no_grassmann_selection": "grassmann_consulted_in_selection: false" in spec,
        "route_closure": "CLOSE_V7_GRASSMANN_PRIMARY_ROUTE_NO_TASK_EXPANSION" in spec,
        "next_only": "IMPLEMENT_V7_7_8_LONG_RANGE_DIFFICULTY_REPLAN_HARNESS_NO_LAUNCH" in spec,
        "bimodal_finding": "bimodal" in flat,
        "fairness_finding": "fairly trained conventional baseline" in flat or "fairly and adequately trained" in flat,
    }
    failed = [name for name, ok in checks.items() if not ok]
    if failed:
        raise ValueError(f"validation failed: {failed}")
    decisions = list(csv.DictReader(
        (ROOT / "DECISIONS.v7.7.7.tsv").open(encoding="utf-8"), delimiter="\t"
    ))
    if len(decisions) != 12 or any(row["status"] != "PASS" for row in decisions) \
            or any(row["state"] != "FROZEN" for row in decisions):
        raise ValueError("decision ledger mismatch")
    manifest = ROOT / "MANIFEST.v7.7.7.sha256"
    if manifest.exists():
        rows = [line.split(maxsplit=1) for line in manifest.read_text(encoding="utf-8").splitlines() if line]
        if len(rows) != len(REQUIRED):
            raise ValueError("manifest count mismatch")
        for expected, name in rows:
            if sha256(ROOT / name) != expected:
                raise ValueError(f"manifest mismatch: {name}")
    print(json.dumps({
        "protocol_version": "v7.7.7",
        "status": "PASS",
        "source_required_status": "LONG_RANGE_TASK_DIFFICULTY_UNRESOLVED",
        "baseline_handed_source_positions": False,
        "gpu_authorized": False,
        "route_closure_branch": "CLOSE_V7_GRASSMANN_PRIMARY_ROUTE_NO_TASK_EXPANSION",
        "manifest_entries": len(REQUIRED),
        "next_authorized_stage": "IMPLEMENT_V7_7_8_LONG_RANGE_DIFFICULTY_REPLAN_HARNESS_NO_LAUNCH",
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
