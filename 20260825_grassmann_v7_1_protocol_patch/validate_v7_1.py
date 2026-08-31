from __future__ import annotations

import argparse
import csv
import hashlib
import json
import py_compile
from pathlib import Path


ROOT = Path(__file__).resolve().parent
ROOT_FILES = (
    "PROTOCOL.v7.1.0.md",
    "DATA_BRANCH_DECISION.v7.1.0.md",
    "GRID_SPEC.v7.1.0.yaml",
    "PANEL_SPEC.v7.1.0.yaml",
    "DECISIONS.v7.1.0.tsv",
    "METRIC_DEFINITIONS.v7.1.0.md",
    "validate_v7_1.py",
    "write_manifest_v7_1.py",
)
P0_FILES = (
    "p0/profile_models_v7_1.py",
    "p0/run_t03_gpu1_v7_1_nonintrusive.sh",
    "p0/build_panel_manifest_v7_1.py",
    "p0/build_compute_contract_v7_1.py",
    "p0/assess_p0_v7_1.py",
    "p0/tests/test_v7_1.py",
)
SERVER_FILES = ("server_ops/SERVER_STEPS.v7.1.0.md",)
MANIFEST_FILES = ROOT_FILES + P0_FILES + SERVER_FILES


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_manifest() -> None:
    lines = [f"{sha256(ROOT / name)}  {name}" for name in MANIFEST_FILES]
    (ROOT / "MANIFEST.v7.1.0.sha256").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-manifest", action="store_true")
    args = parser.parse_args()
    missing = [name for name in MANIFEST_FILES if not (ROOT / name).is_file()]
    if missing:
        raise SystemExit(f"missing v7.1.0 files: {missing}")

    grid = json.loads((ROOT / "GRID_SPEC.v7.1.0.yaml").read_text(encoding="utf-8"))
    panel = json.loads((ROOT / "PANEL_SPEC.v7.1.0.yaml").read_text(encoding="utf-8"))
    decisions = {row["name"]: row["value"] for row in read_tsv(ROOT / "DECISIONS.v7.1.0.tsv")}
    assert grid["protocol_version"] == panel["protocol_version"] == "v7.1.0"
    assert grid["models"]["attention_contract"]["max_keys_per_query"] == 256
    assert float(decisions["local_attention_window"]) == 256
    assert grid["A1_chr22"]["diagnostic"][0]["rates"] == [0.5, 0.9, 0.99]
    assert panel["site_selection"]["force_sequence_length"] is False
    assert panel["site_selection"]["snpbag_compatibility_expected_L"] == 81920
    assert panel["split"]["split_before_synthesis"] is True
    assert panel["generator"]["hgdp_may_be_used_as_donor"] is False
    assert grid["A1_chr22"]["primary_gate"]["stress_cells_can_rescue_primary"] is False
    assert grid["T04_planning"]["runs_per_cell"] == 10

    for name in P0_FILES:
        path = ROOT / name
        if path.suffix == ".py":
            py_compile.compile(str(path), doraise=True)
    if args.write_manifest:
        write_manifest()
    print(json.dumps({
        "status": "PASS",
        "protocol_version": "v7.1.0",
        "manifest_entries": len(MANIFEST_FILES),
        "attention_window": 256,
        "chr22_L_policy": "data_derived",
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
