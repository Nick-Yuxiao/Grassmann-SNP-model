from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


EXPECTED = {
    "all_source": (154850, 4091),
    "joint_release": (154850, 3264),
    "donor_train": (154850, 2247),
    "donor_validation": (154850, 249),
    "hgdp_primary": (154850, 768),
}
SOURCE_SHA256 = "09e50c698522d19bbc2c43a2ea2d29b290d63cf8c1d6983c5198dcb4289ff3fa"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit v7.1.4 materialized panel inventory")
    parser.add_argument("--inventory", type=Path, required=True)
    parser.add_argument("--source-bcf", type=Path, required=True)
    parser.add_argument("--frozen-ids", type=Path, required=True)
    parser.add_argument("--actual-ids", type=Path, required=True)
    parser.add_argument("--format-ids", required=True)
    parser.add_argument("--info-ids", required=True)
    parser.add_argument("--missing-sites", type=int, required=True)
    parser.add_argument("--unphased-sites", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    with args.inventory.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    by_role = {row["role"]: row for row in rows}
    if set(by_role) != set(EXPECTED) or len(rows) != len(EXPECTED):
        raise SystemExit(f"artifact roles mismatch: {sorted(by_role)}")

    for role, (expected_variants, expected_samples) in EXPECTED.items():
        row = by_role[role]
        if int(row["actual_variants"]) != expected_variants:
            raise SystemExit(f"{role}: variant count mismatch")
        if int(row["actual_samples"]) != expected_samples:
            raise SystemExit(f"{role}: sample count mismatch")
        if row["sample_set_match"] != "PASS":
            raise SystemExit(f"{role}: sample-set mismatch")
        artifact = Path(row["artifact"])
        artifact_path = artifact if artifact.is_absolute() else args.inventory.parent / artifact
        if not artifact_path.is_file():
            raise SystemExit(f"{role}: artifact missing: {artifact_path}")

    source_sha = sha256(args.source_bcf)
    if source_sha != SOURCE_SHA256:
        raise SystemExit(f"source SHA-256 mismatch: {source_sha}")
    frozen_ids_sha = sha256(args.frozen_ids)
    actual_ids_sha = sha256(args.actual_ids)
    if frozen_ids_sha != actual_ids_sha:
        raise SystemExit("materialized variant order/key list differs from frozen IDs")
    if args.format_ids != "GT":
        raise SystemExit(f"unexpected FORMAT IDs after cleanup: {args.format_ids!r}")
    if args.info_ids != "":
        raise SystemExit(f"unexpected INFO IDs after cleanup: {args.info_ids!r}")
    if args.missing_sites != 0:
        raise SystemExit(f"joint release has missing-GT sites: {args.missing_sites}")
    if args.unphased_sites != 0:
        raise SystemExit(f"joint release has unphased-GT sites: {args.unphased_sites}")

    audit = {
        "schema_version": "1.0",
        "protocol_version": "v7.1.4",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PASS",
        "variant_match": "CHROM_POS_REF_ALT_EXACT",
        "source_bcf_sha256": source_sha,
        "final_variant_ids_sha256": frozen_ids_sha,
        "materialized_variant_ids_sha256": actual_ids_sha,
        "format_ids": ["GT"],
        "info_ids": [],
        "joint_release_qc": {
            "missing_gt_sites": args.missing_sites,
            "unphased_gt_sites": args.unphased_sites,
        },
        "artifacts": {
            role: {
                "file": by_role[role]["artifact"],
                "variants": int(by_role[role]["actual_variants"]),
                "samples": int(by_role[role]["actual_samples"]),
                "sample_set_match": True,
            }
            for role in EXPECTED
        },
        "source_index_policy": "READ_ONLY_WARNING_RECORDED",
        "gpu_used": False,
    }
    args.output.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
