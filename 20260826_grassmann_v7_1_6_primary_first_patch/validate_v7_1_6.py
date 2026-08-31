from __future__ import annotations

import argparse
import csv
import hashlib
import json
import py_compile
from pathlib import Path


ROOT = Path(__file__).resolve().parent
FILES = (
    "PROTOCOL_ADDENDUM.v7.1.6.md",
    "PRIMARY_FIRST_GRID.v7.1.6.yaml",
    "DECISIONS.v7.1.6.tsv",
    "validate_v7_1_6.py",
    "p0/build_t04_primary_first_v7_1_6.py",
    "p0/tests/test_v7_1_6.py",
    "server_ops/SERVER_STEPS.v7.1.6.md",
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
        raise SystemExit(f"missing v7.1.6 files: {missing}")
    grid = json.loads((ROOT / "PRIMARY_FIRST_GRID.v7.1.6.yaml").read_text(encoding="utf-8"))
    with (ROOT / "DECISIONS.v7.1.6.tsv").open(encoding="utf-8", newline="") as handle:
        decisions = {row["name"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}
    primary = grid["primary_A1"]
    pilot = grid["convergence_pilot"]
    assert grid["protocol_version"] == "v7.1.6"
    assert grid["sequence_length"] == 154850
    assert grid["allowed_physical_gpus"] == [1, 2, 3, 4, 5, 6]
    assert grid["forbidden_physical_gpus"] == [0]
    assert len(primary["masks"]) * len(primary["fairness"]) * len(grid["models"]) * len(primary["data_seeds"]) * len(primary["init_seeds"]) == 120
    assert len(grid["models"]) * len(pilot["masks"]) * len(pilot["data_seeds"]) * len(pilot["init_seeds"]) == 12
    assert pilot["hgdp_access"] == "FORBIDDEN"
    assert decisions["per_run_early_stopping"] == "FORBIDDEN"
    for name in FILES:
        path = ROOT / name
        if path.suffix == ".py":
            py_compile.compile(str(path), doraise=True)
    if args.write_manifest:
        lines = [f"{sha256(ROOT / name)}  {name}" for name in FILES]
        with (ROOT / "MANIFEST.v7.1.6.sha256").open("w", encoding="utf-8", newline="\n") as handle:
            handle.write("\n".join(lines) + "\n")
    print(json.dumps({
        "status": "PASS",
        "protocol_version": "v7.1.6",
        "convergence_pilot_runs": 12,
        "primary_runs": 120,
        "allowed_physical_gpus": grid["allowed_physical_gpus"],
        "gpu0": "FORBIDDEN",
        "per_run_early_stopping": decisions["per_run_early_stopping"],
        "manifest_entries": len(FILES),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
