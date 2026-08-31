from __future__ import annotations

import argparse
import csv
import hashlib
import json
import py_compile
from pathlib import Path


ROOT = Path(__file__).resolve().parent
FILES = (
    "PROTOCOL_ADDENDUM.v7.1.5.md",
    "PILOT_GRID.v7.1.5.yaml",
    "DECISIONS.v7.1.5.tsv",
    "validate_v7_1_5.py",
    "p0/build_t04_contract_v7_1_5.py",
    "p0/tests/test_v7_1_5.py",
    "server_ops/SERVER_STEPS.v7.1.5.md",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-manifest", action="store_true")
    args = parser.parse_args()
    missing = [name for name in FILES if not (ROOT / name).is_file()]
    if missing:
        raise SystemExit(f"missing v7.1.5 files: {missing}")
    grid = json.loads((ROOT / "PILOT_GRID.v7.1.5.yaml").read_text(encoding="utf-8"))
    with (ROOT / "DECISIONS.v7.1.5.tsv").open(encoding="utf-8", newline="") as handle:
        decisions = {row["name"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}
    assert grid["protocol_version"] == "v7.1.5"
    assert grid["sequence_length"] == 154850
    assert grid["train_steps_per_run"] == 2000
    assert len(grid["models"]) == 3
    assert len(grid["mask_cells"]) == 4
    assert grid["hgdp_access"] == "FORBIDDEN"
    assert grid["pilot_capacity_envelope"]["physical_gpu"] == 1
    assert decisions["full_A1_status"] == "UNFUNDED_UNRESERVED"
    for name in FILES:
        path = ROOT / name
        if path.suffix == ".py":
            py_compile.compile(str(path), doraise=True)
    if args.write_manifest:
        lines = [f"{sha256(ROOT / name)}  {name}" for name in FILES]
        with (ROOT / "MANIFEST.v7.1.5.sha256").open("w", encoding="utf-8", newline="\n") as handle:
            handle.write("\n".join(lines) + "\n")
    print(json.dumps({
        "status": "PASS",
        "protocol_version": "v7.1.5",
        "pilot_runs": len(grid["models"]) * len(grid["mask_cells"]) * len(grid["pilot_data_seeds"]) * len(grid["pilot_init_seeds"]),
        "pilot_steps": grid["train_steps_per_run"],
        "hgdp_access": grid["hgdp_access"],
        "full_A1_status": decisions["full_A1_status"],
        "manifest_entries": len(FILES),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
