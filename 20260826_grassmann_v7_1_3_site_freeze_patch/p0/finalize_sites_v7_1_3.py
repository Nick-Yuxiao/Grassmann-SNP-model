from __future__ import annotations

import argparse
import collections
import csv
import gzip
import hashlib
import json
from pathlib import Path


FIELDS = ["chrom", "pos", "ref", "alt", "ac", "an", "af", "maf", "f_missing"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_metrics(path: Path):
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{path}:{line_number}: expected 9 fields, observed {len(fields)}")
            yield dict(zip(FIELDS, fields))


def key(row: dict[str, str]) -> tuple[str, int, str, str]:
    return row["chrom"], int(row["pos"]), row["ref"], row["alt"]


def main() -> None:
    parser = argparse.ArgumentParser(description="Materialize the frozen v7.1.3 chr22 site list")
    parser.add_argument("--donor-sites", type=Path, required=True)
    parser.add_argument("--hgdp-metrics", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--expected-donor-candidates", type=int, default=155586)
    parser.add_argument("--expected-duplicate-positions", type=int, default=368)
    parser.add_argument("--expected-duplicate-records", type=int, default=736)
    parser.add_argument("--expected-final", type=int, default=154850)
    parser.add_argument("--expected-hgdp-monomorphic", type=int, default=104)
    parser.add_argument("--expected-hgdp-an", type=int, default=1536)
    args = parser.parse_args()

    if args.output_dir.exists():
        raise SystemExit(f"refusing to overwrite output directory: {args.output_dir}")

    donor_rows = list(read_metrics(args.donor_sites))
    if len(donor_rows) != args.expected_donor_candidates:
        raise SystemExit(
            f"donor candidate count mismatch: {len(donor_rows)} != {args.expected_donor_candidates}"
        )
    position_counts = collections.Counter((row["chrom"], int(row["pos"])) for row in donor_rows)
    duplicate_positions = {position for position, count in position_counts.items() if count > 1}
    duplicate_records = sum(position_counts[position] for position in duplicate_positions)
    if len(duplicate_positions) != args.expected_duplicate_positions:
        raise SystemExit("duplicate-position count mismatch")
    if duplicate_records != args.expected_duplicate_records:
        raise SystemExit("duplicate-record count mismatch")

    eligible_rows = [
        row for row in donor_rows if (row["chrom"], int(row["pos"])) not in duplicate_positions
    ]
    eligible_by_key = {key(row): row for row in eligible_rows}
    if len(eligible_by_key) != len(eligible_rows):
        raise SystemExit("duplicate exact donor keys remain")
    if len(eligible_rows) != args.expected_final:
        raise SystemExit(f"duplicate-free candidate count mismatch: {len(eligible_rows)}")

    hgdp_by_key: dict[tuple[str, int, str, str], dict[str, str]] = {}
    for row in read_metrics(args.hgdp_metrics):
        row_key = key(row)
        if row_key not in eligible_by_key:
            continue
        if row_key in hgdp_by_key:
            raise SystemExit(f"duplicate exact HGDP key: {row_key}")
        hgdp_by_key[row_key] = row

    missing_keys = set(eligible_by_key) - set(hgdp_by_key)
    if missing_keys:
        raise SystemExit(f"HGDP exact keys missing: {sorted(missing_keys)[:10]}")

    incomplete = []
    for row_key, row in hgdp_by_key.items():
        if int(row["an"]) != args.expected_hgdp_an or float(row["f_missing"]) != 0.0:
            incomplete.append(row_key)
    if incomplete:
        raise SystemExit(f"HGDP incomplete sites: {sorted(incomplete)[:10]}")

    monomorphic = sum(int(hgdp_by_key[row_key]["ac"]) == 0 for row_key in eligible_by_key)
    if monomorphic != args.expected_hgdp_monomorphic:
        raise SystemExit(f"HGDP monomorphic count mismatch: {monomorphic}")

    previous_position = -1
    for row in eligible_rows:
        position = int(row["pos"])
        if position <= previous_position:
            raise SystemExit("final sites are not strictly coordinate increasing")
        previous_position = position

    args.output_dir.mkdir(parents=True, exist_ok=False)
    variants_path = args.output_dir / "FINAL_VARIANTS.v7.1.3.tsv"
    ids_path = args.output_dir / "FINAL_VARIANT_IDS.v7.1.3.txt"
    audit_path = args.output_dir / "FINAL_SITE_FREEZE.v7.1.3.json"

    output_fields = [
        "variant_id", "chrom", "pos", "ref", "alt",
        "donor_ac", "donor_an", "donor_af", "donor_maf", "donor_f_missing",
        "hgdp_ac", "hgdp_an", "hgdp_af", "hgdp_maf", "hgdp_f_missing",
        "hgdp_allele_state",
    ]
    with variants_path.open("w", encoding="utf-8", newline="") as variants_handle, ids_path.open(
        "w", encoding="utf-8", newline=""
    ) as ids_handle:
        writer = csv.DictWriter(
            variants_handle, fieldnames=output_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for donor in eligible_rows:
            row_key = key(donor)
            hgdp = hgdp_by_key[row_key]
            variant_id = f"{donor['chrom']}:{donor['pos']}:{donor['ref']}:{donor['alt']}"
            hgdp_ac = int(hgdp["ac"])
            state = "monomorphic_reference" if hgdp_ac == 0 else "polymorphic"
            writer.writerow({
                "variant_id": variant_id,
                "chrom": donor["chrom"],
                "pos": donor["pos"],
                "ref": donor["ref"],
                "alt": donor["alt"],
                "donor_ac": donor["ac"],
                "donor_an": donor["an"],
                "donor_af": donor["af"],
                "donor_maf": donor["maf"],
                "donor_f_missing": donor["f_missing"],
                "hgdp_ac": hgdp["ac"],
                "hgdp_an": hgdp["an"],
                "hgdp_af": hgdp["af"],
                "hgdp_maf": hgdp["maf"],
                "hgdp_f_missing": hgdp["f_missing"],
                "hgdp_allele_state": state,
            })
            ids_handle.write(variant_id + "\n")

    audit = {
        "schema_version": "1.0",
        "protocol_version": "v7.1.3",
        "status": "PASS",
        "source_filter_semantics": "SOURCE_PREFILTERED_FILTER_UNSET",
        "donor_candidates_before_duplicate_rejection": len(donor_rows),
        "duplicate_positions_rejected": len(duplicate_positions),
        "records_at_duplicate_positions_rejected": duplicate_records,
        "final_L": len(eligible_rows),
        "hgdp_exact_key_matches": len(hgdp_by_key),
        "hgdp_expected_AN": args.expected_hgdp_an,
        "hgdp_incomplete_sites": 0,
        "hgdp_monomorphic_reference_sites_retained": monomorphic,
        "hgdp_polymorphic_sites": len(eligible_rows) - monomorphic,
        "hgdp_AC_used_for_selection": False,
        "length_was_forced": False,
        "snpbag_expected_L": 81920,
        "snpbag_compatibility_status": "DATA_DERIVED_L_DIFFERS",
        "input_sha256": {
            "donor_sites": sha256(args.donor_sites),
            "hgdp_metrics": sha256(args.hgdp_metrics),
        },
    }
    audit_path.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    files = (variants_path, ids_path, audit_path)
    manifest_path = args.output_dir / "FINAL_SITE_FREEZE.v7.1.3.sha256"
    manifest_path.write_text(
        "\n".join(f"{sha256(path)}  {path.name}" for path in files) + "\n", encoding="utf-8"
    )
    print(json.dumps(audit, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
