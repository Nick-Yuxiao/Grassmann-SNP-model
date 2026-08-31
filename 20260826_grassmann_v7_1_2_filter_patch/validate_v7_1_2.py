from __future__ import annotations

import argparse
import csv
import hashlib
import json
import py_compile
from pathlib import Path


ROOT = Path(__file__).resolve().parent
FILES = (
    "PROTOCOL_ADDENDUM.v7.1.2.md",
    "PANEL_SPEC.v7.1.2.yaml",
    "DECISIONS.v7.1.2.tsv",
    "validate_v7_1_2.py",
    "p0/audit_site_stage1_v7_1_2.py",
    "p0/run_site_stage1_v7_1_2.sh",
    "p0/build_panel_manifest_v7_1_2.py",
    "p0/tests/test_v7_1_2.py",
    "server_ops/SERVER_STEPS.v7.1.2.md",
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
        raise SystemExit(f"missing v7.1.2 files: {missing}")

    panel = json.loads((ROOT / "PANEL_SPEC.v7.1.2.yaml").read_text(encoding="utf-8"))
    with (ROOT / "DECISIONS.v7.1.2.tsv").open(encoding="utf-8", newline="") as handle:
        decisions = {row["name"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}
    site = panel["site_selection"]
    assert panel["protocol_version"] == "v7.1.2"
    assert panel["source_bcf_sha256"] == "09e50c698522d19bbc2c43a2ea2d29b290d63cf8c1d6983c5198dcb4289ff3fa"
    assert site["accepted_vcf_filter_value"] == "."
    assert site["literal_PASS_available"] is False
    assert site["maf_operator"] == ">" and site["maf_threshold"] == 0.01
    assert site["force_sequence_length"] is False
    assert decisions["record_filter_rule"] == "SOURCE_PREFILTERED_FILTER_UNSET"
    assert decisions["duplicate_position_rule"] == "reject_all_records_at_duplicate_positions"

    for name in FILES:
        path = ROOT / name
        if path.suffix == ".py":
            py_compile.compile(str(path), doraise=True)
    shell_bytes = (ROOT / "p0/run_site_stage1_v7_1_2.sh").read_bytes()
    assert shell_bytes.startswith(b"#!/usr/bin/env bash\n")
    assert b"\r\n" not in shell_bytes

    if args.write_manifest:
        lines = [f"{sha256(ROOT / name)}  {name}" for name in FILES]
        with (ROOT / "MANIFEST.v7.1.2.sha256").open("w", encoding="utf-8", newline="\n") as handle:
            handle.write("\n".join(lines) + "\n")

    print(json.dumps({
        "status": "PASS",
        "protocol_version": "v7.1.2",
        "filter_semantics": "SOURCE_PREFILTERED_FILTER_UNSET",
        "maf_rule": "GT_RECOMPUTED_STRICT_GT_0.01",
        "manifest_entries": len(FILES),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
