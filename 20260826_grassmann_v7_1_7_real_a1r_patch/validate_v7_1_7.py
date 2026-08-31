from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    root = args.root.resolve()
    grid = json.loads((root / "A1R_GRID.v7.1.7.yaml").read_text(encoding="utf-8"))
    protocol = (root / "PROTOCOL_ADDENDUM.v7.1.7.md").read_text(encoding="utf-8")

    assert grid["protocol_version"] == "v7.1.7"
    assert grid["sequence_length"] == 154850
    assert grid["allowed_physical_gpus"] == [1, 2, 3, 4, 5, 6]
    assert 0 in grid["forbidden_physical_gpus"]
    assert grid["data_contract"]["hgdp_access"] == "FORBIDDEN"
    assert grid["convergence_pilot"]["runs"] == 12
    assert grid["primary_A1R_100pct"]["runs"] == 120
    assert grid["sample_size_diagnostic"]["added_runs"] == 24
    assert grid["sample_size_diagnostic"]["sample_counts"] == [562, 1124]
    assert grid["inference"]["delta_min"] == 0.010
    assert grid["inference"]["mask_seeds_are_independent_biological_datasets"] is False
    for phrase in (
        "NO_GO_EQUIVALENT_OR_WORSE",
        "INCONCLUSIVE_SAMPLE_LIMITED_HAPNEST",
        "every one of the 120 primary runs",
        "HGDP access is forbidden",
    ):
        assert phrase in protocol, phrase

    manifest = root / "MANIFEST.v7.1.7.sha256"
    entries = 0
    for line in manifest.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        expected, relative = line.split("  ", 1)
        path = root / relative
        assert path.is_file(), relative
        assert sha256(path) == expected, relative
        entries += 1
    print(json.dumps({
        "status": "PASS",
        "protocol_version": "v7.1.7",
        "pilot_runs": 12,
        "primary_runs": 120,
        "diagnostic_added_runs": 24,
        "post_pilot_total_runs": 144,
        "hgdp_access": "FORBIDDEN",
        "negative_default": "INCONCLUSIVE_SAMPLE_LIMITED_HAPNEST",
        "manifest_entries": entries,
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

