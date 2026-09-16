from __future__ import annotations

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def main() -> None:
    contract = json.loads((HERE / "CONTRACT.json").read_text(encoding="utf-8"))
    assert contract["status"] == "FROZEN_EXECUTION_BLOCKED"
    assert contract["run_authorized"] is False
    assert contract["route"] == "0526_stage2_only"
    assert contract["prohibited_detour"] == "new_foundation_model_attribution_audit"
    assert set(contract["arms"]) == {"A", "B", "C", "D", "E"}
    assert [item["id"] for item in contract["primary_family"]] == ["P1", "P2", "P3"]
    assert contract["multiplicity"]["family"] == ["P1", "P2", "P3"]
    assert contract["uncertainty"]["replicates"] == 5000
    assert abs(sum(contract["split"][key] for key in ("train_fraction", "validation_fraction", "test_fraction")) - 1.0) < 1e-12
    assert all(value is False for key, value in contract["gates"].items() if key != "test_outcome_unopened")
    assert contract["gates"]["test_outcome_unopened"] is True

    manifest = (HERE / "SOURCE_MANIFEST.tsv").read_text(encoding="utf-8").splitlines()
    assert manifest[0] == "evidence_id\tpath\tsha256\trole"
    for row in (line for line in manifest[1:] if line.strip()):
        evidence_id, raw_path, expected, _ = row.split("\t", 3)
        path = Path(raw_path)
        if not path.is_absolute():
            path = ROOT / path
        assert path.is_file(), f"missing {evidence_id}: {path}"
        assert sha256(path) == expected, f"hash mismatch for {evidence_id}: {path}"

    frozen_manifest = HERE / "FROZEN_MANIFEST.sha256"
    if frozen_manifest.is_file():
        for row in frozen_manifest.read_text(encoding="utf-8").splitlines():
            if not row.strip():
                continue
            expected, relative = row.split(maxsplit=1)
            target = HERE / relative.lstrip("*")
            assert target.is_file(), f"missing frozen protocol file: {target}"
            assert sha256(target) == expected.upper(), f"frozen protocol hash mismatch: {target}"

    print("PASS_FROZEN_EXECUTION_BLOCKED: contract is intact; phenotype execution remains NO-GO")


if __name__ == "__main__":
    main()
