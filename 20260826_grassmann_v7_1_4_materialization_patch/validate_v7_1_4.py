from __future__ import annotations

import argparse
import csv
import hashlib
import json
import py_compile
from pathlib import Path


ROOT = Path(__file__).resolve().parent
FILES = (
    "PROTOCOL_ADDENDUM.v7.1.4.md",
    "PANEL_SPEC.v7.1.4.yaml",
    "DECISIONS.v7.1.4.tsv",
    "validate_v7_1_4.py",
    "p0/build_materialization_audit_v7_1_4.py",
    "p0/materialize_panel_v7_1_4.sh",
    "p0/tests/test_v7_1_4.py",
    "server_ops/SERVER_STEPS.v7.1.4.md",
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
        raise SystemExit(f"missing v7.1.4 files: {missing}")

    panel = json.loads((ROOT / "PANEL_SPEC.v7.1.4.yaml").read_text(encoding="utf-8"))
    with (ROOT / "DECISIONS.v7.1.4.tsv").open(encoding="utf-8", newline="") as handle:
        decisions = {row["name"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}
    assert panel["protocol_version"] == "v7.1.4"
    assert panel["expected_L"] == 154850
    assert panel["variant_match"] == "CHROM_POS_REF_ALT_EXACT"
    assert panel["format_policy"]["retain"] == ["GT"]
    assert panel["format_policy"]["remove"] == ["PP"]
    assert panel["info_policy"]["retain"] == []
    assert panel["artifacts"]["joint_release"] == 3264
    assert panel["gpu_required"] is False
    assert decisions["source_index_policy"] == "READ_ONLY"

    for name in FILES:
        path = ROOT / name
        if path.suffix == ".py":
            py_compile.compile(str(path), doraise=True)

    if args.write_manifest:
        lines = [f"{sha256(ROOT / name)}  {name}" for name in FILES]
        (ROOT / "MANIFEST.v7.1.4.sha256").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(json.dumps({
        "status": "PASS",
        "protocol_version": "v7.1.4",
        "variant_match": panel["variant_match"],
        "expected_L": panel["expected_L"],
        "artifact_count": len(panel["artifacts"]),
        "gpu_required": panel["gpu_required"],
        "manifest_entries": len(FILES),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
