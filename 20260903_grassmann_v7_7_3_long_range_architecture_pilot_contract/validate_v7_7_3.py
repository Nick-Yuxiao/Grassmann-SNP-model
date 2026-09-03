from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def fail(message: str) -> None:
    raise SystemExit(f"FAIL: {message}")


def main() -> None:
    spec = (ROOT / "LONG_RANGE_ARCHITECTURE_PILOT_SPEC.v7.7.3.yaml").read_text(encoding="utf-8")
    protocol = (ROOT / "PROTOCOL_ADDENDUM.v7.7.3.md").read_text(encoding="utf-8")
    manifest = (ROOT / "MANIFEST.v7.7.3.sha256").read_text(encoding="utf-8").splitlines()
    required = [
        "gpu_authorized: false", "pilot_execution_authorized: false",
        "architecture_decision_permitted: false", "effective_n: 6",
        "R1_ROUTER", "R3_ROUTER_GRASSMANN",
        "NLL(R1_ROUTER)-NLL(R3_ROUTER_GRASSMANN)",
        "release_pilot_arm_means: false",
    ]
    for item in required:
        if item not in spec:
            fail(f"missing specification: {item}")
    if "LONG_RANGE_TASK_VALIDITY_PASS" not in protocol:
        fail("source PASS requirement absent")
    if "no-op" not in protocol or "target-shuffled" not in protocol:
        fail("required readiness gates absent")
    if len(manifest) != 7:
        fail("unexpected manifest entry count")
    print(json.dumps({
        "protocol_version": "v7.7.3",
        "status": "PASS",
        "effective_n": 6,
        "pilot_execution_cells": 48,
        "gpu_authorized": False,
        "next_authorized_stage": "IMPLEMENT_V7_7_4_LONG_RANGE_PILOT_HARNESS_NO_LAUNCH",
        "manifest_entries": len(manifest),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
