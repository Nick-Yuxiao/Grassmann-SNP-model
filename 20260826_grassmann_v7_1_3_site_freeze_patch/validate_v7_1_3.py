from __future__ import annotations

import argparse
import csv
import hashlib
import json
import py_compile
from pathlib import Path


ROOT = Path(__file__).resolve().parent
FILES = (
    "PROTOCOL_ADDENDUM.v7.1.3.md",
    "PANEL_SPEC.v7.1.3.yaml",
    "DECISIONS.v7.1.3.tsv",
    "validate_v7_1_3.py",
    "p0/finalize_sites_v7_1_3.py",
    "p0/build_panel_manifest_v7_1_3.py",
    "p0/tests/test_v7_1_3.py",
    "server_ops/SERVER_STEPS.v7.1.3.md",
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
        raise SystemExit(f"missing v7.1.3 files: {missing}")

    panel = json.loads((ROOT / "PANEL_SPEC.v7.1.3.yaml").read_text(encoding="utf-8"))
    with (ROOT / "DECISIONS.v7.1.3.tsv").open(encoding="utf-8", newline="") as handle:
        decisions = {row["name"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}
    site = panel["site_selection"]
    reporting = panel["allele_state_reporting"]
    assert panel["protocol_version"] == "v7.1.3"
    assert site["final_sequence_length"] == 154850
    assert site["hgdp_expected_AN"] == 1536
    assert site["hgdp_required_F_MISSING"] == 0.0
    assert site["hgdp_AC_used_for_selection"] is False
    assert site["force_sequence_length"] is False
    assert reporting["hgdp_monomorphic_reference_sites"] == 104
    assert reporting["hgdp_polymorphic_sites"] == 154746
    assert int(decisions["final_chr22_L"]) == 154850
    assert decisions["hgdp_AC_selection"] == "FORBIDDEN"

    for name in FILES:
        path = ROOT / name
        if path.suffix == ".py":
            py_compile.compile(str(path), doraise=True)

    if args.write_manifest:
        lines = [f"{sha256(ROOT / name)}  {name}" for name in FILES]
        with (ROOT / "MANIFEST.v7.1.3.sha256").open("w", encoding="utf-8", newline="\n") as handle:
            handle.write("\n".join(lines) + "\n")

    print(json.dumps({
        "status": "PASS",
        "protocol_version": "v7.1.3",
        "final_L": 154850,
        "hgdp_callability": "AN_1536_AND_F_MISSING_0",
        "hgdp_monomorphic_sites_retained": 104,
        "manifest_entries": len(FILES),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
