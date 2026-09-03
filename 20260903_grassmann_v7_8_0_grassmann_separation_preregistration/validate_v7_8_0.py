from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
REQUIRED = (
    "PROTOCOL_ADDENDUM.v7.8.0.md",
    "GRASSMANN_SEPARATION_PREREG_SPEC.v7.8.0.yaml",
    "SEPARATION_ARM_MAP.v7.8.0.tsv",
    "DECISIONS.v7.8.0.tsv",
    "validate_v7_8_0.py",
    "p0/build_separation_prereg_readiness_v7_8_0.py",
    "p0/tests/test_v7_8_0.py",
    "server_ops/SERVER_STEPS.v7.8.0.md",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    missing = [item for item in REQUIRED if not (ROOT / item).is_file()]
    if missing:
        raise ValueError(f"missing files: {missing}")
    spec = (ROOT / "GRASSMANN_SEPARATION_PREREG_SPEC.v7.8.0.yaml").read_text(encoding="utf-8")
    protocol = (ROOT / "PROTOCOL_ADDENDUM.v7.8.0.md").read_text(encoding="utf-8")
    flat = " ".join(protocol.split()).lower()
    checks = {
        "separation_estimand": "SEPARATION" in spec and "NLL_conv_best_minus_NLL_grassmann" in spec,
        "lr1_not_required": "lr1_required: false" in spec,
        "gpu_gated_not_authorized": "gpu_authorized: false" in spec and "gpu_gated: true" in spec,
        "fair_suite": "not_a_strawman: true" in spec,
        "convergence_audit": "conventional_convergence_compute_sufficiency_audit_required: true" in spec,
        "no_taskshopping": "reverse_engineered_from_grassmann: false" in spec,
        "route_closure": "CLOSE_V7_GRASSMANN_PRIMARY_ROUTE_NO_TASK_EXPANSION_NO_GPU" in spec,
        "next_only": "IMPLEMENT_V7_8_1_CPU_GRASSMANN_SEPARATION_HARNESS_NO_LAUNCH" in spec,
        "a1r_trap_named": "a1-r trap" in flat,
        "affirmative_a1": "affirmative a1" in flat,
    }
    failed = [name for name, ok in checks.items() if not ok]
    if failed:
        raise ValueError(f"validation failed: {failed}")
    decisions = list(csv.DictReader(
        (ROOT / "DECISIONS.v7.8.0.tsv").open(encoding="utf-8"), delimiter="\t"
    ))
    if len(decisions) != 14 or any(row["status"] != "PASS" for row in decisions) \
            or any(row["state"] != "FROZEN" for row in decisions):
        raise ValueError("decision ledger mismatch")
    arms = list(csv.DictReader(
        (ROOT / "SEPARATION_ARM_MAP.v7.8.0.tsv").open(encoding="utf-8"), delimiter="\t"
    ))
    conventional = [r for r in arms if r["role"] == "conventional_comparator"]
    if len(conventional) < 3 or not any(r["arm_id"] == "GRASSMANN" for r in arms) \
            or any(r["execution_authorized"] != "FALSE" for r in arms):
        raise ValueError("arm map mismatch")
    manifest = ROOT / "MANIFEST.v7.8.0.sha256"
    if manifest.exists():
        rows = [line.split(maxsplit=1) for line in manifest.read_text(encoding="utf-8").splitlines() if line]
        if len(rows) != len(REQUIRED):
            raise ValueError("manifest count mismatch")
        for expected, name in rows:
            if sha256(ROOT / name) != expected:
                raise ValueError(f"manifest mismatch: {name}")
    print(json.dumps({
        "protocol_version": "v7.8.0",
        "status": "PASS",
        "estimand": "SEPARATION",
        "lr1_required": False,
        "gpu_authorized": False,
        "gpu_gated": True,
        "conventional_comparator_families": len(conventional),
        "route_closure_branch": "CLOSE_V7_GRASSMANN_PRIMARY_ROUTE_NO_TASK_EXPANSION_NO_GPU",
        "manifest_entries": len(REQUIRED),
        "next_authorized_stage": "IMPLEMENT_V7_8_1_CPU_GRASSMANN_SEPARATION_HARNESS_NO_LAUNCH",
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
