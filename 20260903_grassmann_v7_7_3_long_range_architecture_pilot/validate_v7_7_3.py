from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
REQUIRED = (
    "PROTOCOL_ADDENDUM.v7.7.3.md",
    "LONG_RANGE_ARCHITECTURE_PILOT_SPEC.v7.7.3.yaml",
    "ARCHITECTURE_PILOT_ARM_MAP.v7.7.3.tsv",
    "DECISIONS.v7.7.3.tsv",
    "validate_v7_7_3.py",
    "p0/build_architecture_pilot_readiness_v7_7_3.py",
    "p0/tests/test_v7_7_3.py",
    "server_ops/SERVER_STEPS.v7.7.3.md",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    missing = [item for item in REQUIRED if not (ROOT / item).is_file()]
    if missing:
        raise ValueError(f"missing files: {missing}")
    protocol = (ROOT / "PROTOCOL_ADDENDUM.v7.7.3.md").read_text(encoding="utf-8")
    spec = (ROOT / "LONG_RANGE_ARCHITECTURE_PILOT_SPEC.v7.7.3.yaml").read_text(encoding="utf-8")
    flat = " ".join(protocol.split()).lower()
    checks = {
        "requires_task_validity_pass": "LONG_RANGE_TASK_VALIDITY_PASS" in protocol,
        "task_validity_not_grassmann": "not a grassmann" in flat,
        "router_not_source_positions": "not handed the two source positions" in flat,
        "primary_contrast": "nll(lr_a01) - nll(lr_a11)" in flat,
        "blinded_pilot": "blinded" in flat and "arm means are withheld" in flat,
        "no_gpu": "gpu_authorized: false" in spec,
        "no_grassmann_training": "grassmann_training_authorized: false" in spec,
        "no_architecture_decision": "architecture_decision_permitted: false" in spec,
        "next_only": "IMPLEMENT_V7_7_4_CPU_BLINDED_VARIANCE_PILOT_NO_GPU" in spec,
    }
    failed = [name for name, ok in checks.items() if not ok]
    if failed:
        raise ValueError(f"validation failed: {failed}")
    decisions = list(csv.DictReader(
        (ROOT / "DECISIONS.v7.7.3.tsv").open(encoding="utf-8"), delimiter="\t"
    ))
    if len(decisions) != 16 or any(row["status"] != "FROZEN" for row in decisions):
        raise ValueError("decision ledger mismatch")
    arms = list(csv.DictReader(
        (ROOT / "ARCHITECTURE_PILOT_ARM_MAP.v7.7.3.tsv").open(encoding="utf-8"), delimiter="\t"
    ))
    decision_cells = [row for row in arms if row["architecture_decision_eligible"] == "TRUE"]
    if len(decision_cells) != 4 or any(row["execution_authorized"] != "FALSE" for row in arms):
        raise ValueError("arm map mismatch")
    manifest = ROOT / "MANIFEST.v7.7.3.sha256"
    if manifest.exists():
        rows = [line.split(maxsplit=1) for line in manifest.read_text(encoding="utf-8").splitlines() if line]
        if len(rows) != len(REQUIRED):
            raise ValueError("manifest count mismatch")
        for expected, name in rows:
            if sha256(ROOT / name) != expected:
                raise ValueError(f"manifest mismatch: {name}")
    print(json.dumps({
        "protocol_version": "v7.7.3",
        "status": "PASS",
        "primary_contrast": "NLL_LR_A01_MINUS_NLL_LR_A11",
        "practical_margin_nats_per_target": 0.010,
        "gpu_authorized": False,
        "grassmann_training_authorized": False,
        "architecture_decision_permitted": False,
        "decision_eligible_cells": len(decision_cells),
        "manifest_entries": len(REQUIRED),
        "next_authorized_stage": "IMPLEMENT_V7_7_4_CPU_BLINDED_VARIANCE_PILOT_NO_GPU",
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
