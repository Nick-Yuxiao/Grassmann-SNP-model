from __future__ import annotations

import argparse
import csv
import hashlib
import json
import py_compile
from pathlib import Path


ROOT = Path(__file__).resolve().parent
FILES = (
    "PROTOCOL_ADDENDUM.v7.1.1.md",
    "PANEL_SPEC.v7.1.1.yaml",
    "DECISIONS.v7.1.1.tsv",
    "validate_v7_1_1.py",
    "p0/freeze_samples_v7_1_1.py",
    "p0/build_panel_manifest_v7_1_1.py",
    "p0/tests/test_v7_1_1.py",
    "server_ops/SERVER_STEPS.v7.1.1.md",
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
        raise SystemExit(f"missing v7.1.1 patch files: {missing}")

    panel = json.loads((ROOT / "PANEL_SPEC.v7.1.1.yaml").read_text(encoding="utf-8"))
    with (ROOT / "DECISIONS.v7.1.1.tsv").open(encoding="utf-8", newline="") as handle:
        decisions = {row["name"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}

    assert panel["protocol_version"] == "v7.1.1"
    assert panel["sample_release_rule"]["expected_1KGP_count"] == 2496
    assert panel["sample_release_rule"]["expected_HGDP_count"] == 768
    assert panel["sample_release_rule"]["all_samples_related_is_additional_exclusion"] is False
    assert panel["split"]["validation_fraction"] == 0.10
    assert panel["split"]["seed"] == 71001
    assert panel["partitions"]["hgdp_snpbag_calibration"]["current_individual_count"] == 0
    assert panel["partitions"]["hgdp_snpbag_calibration"]["surrogate_allowed"] is False
    assert decisions["snpbag_calibration_status"] == "NOT_COMPARABLE_SOURCE_MISMATCH"
    assert decisions["snpbag_calibration_surrogate"] == "FORBIDDEN"

    for name in FILES:
        path = ROOT / name
        if path.suffix == ".py":
            py_compile.compile(str(path), doraise=True)

    if args.write_manifest:
        lines = [f"{sha256(ROOT / name)}  {name}" for name in FILES]
        with (ROOT / "MANIFEST.v7.1.1.sha256").open("w", encoding="utf-8", newline="\n") as handle:
            handle.write("\n".join(lines) + "\n")

    print(json.dumps({
        "status": "PASS",
        "protocol_version": "v7.1.1",
        "patch_manifest_entries": len(FILES),
        "release_counts": {"1KGP": 2496, "HGDP": 768},
        "donor_validation_fraction": 0.10,
        "calibration_status": "NOT_COMPARABLE_SOURCE_MISMATCH",
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
