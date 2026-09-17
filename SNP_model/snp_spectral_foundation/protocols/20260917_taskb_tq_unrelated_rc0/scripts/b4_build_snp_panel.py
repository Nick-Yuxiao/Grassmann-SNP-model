#!/usr/bin/env python3
"""B4: build a thinned genotype panel from PLINK hard-calls for the B arm.

Variants are thinned by physical distance, read directly out of the .bed with a
byte lookup table, filtered on train-only MAF and missingness, and stored
variant-major so the qualification pass can stream them.

Allele orientation note: this cohort's .bed files were converted from BGEN
without an explicit --bgen REF/ALT mode, so A1/A2 may not match the official
REF/ALT. The panel therefore records dosage as "count of the .bim A1 allele" and
nothing here is compared against an external reference.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_split import open_text, read_fam_ids, sha256_file, stable_hex  # noqa: E402

MISSING = np.int8(-1)

# Long-range LD / inversion regions excluded by default (GRCh37 coordinates).
DEFAULT_EXCLUSIONS = [
    ("6", 25_000_000, 34_000_000),   # MHC
    ("8", 8_000_000, 12_000_000),    # 8p23.1 inversion
    ("17", 40_000_000, 45_000_000),  # 17q21.31 inversion
]


def dosage_lookup() -> np.ndarray:
    """Byte -> four A1 dosages. PLINK codes: 00=A1A1, 01=missing, 10=het, 11=A2A2."""
    code_to_dosage = np.array([2, -1, 1, 0], dtype=np.int8)
    table = np.empty((256, 4), dtype=np.int8)
    for byte in range(256):
        for slot in range(4):
            table[byte, slot] = code_to_dosage[(byte >> (2 * slot)) & 0b11]
    return table


def read_bim(path: Path) -> list[dict[str, object]]:
    variants: list[dict[str, object]] = []
    with open_text(path) as handle:
        for index, line in enumerate(handle):
            fields = line.split()
            if len(fields) < 6:
                continue
            try:
                position = int(fields[3])
            except ValueError:
                continue
            variants.append({
                "index": index,
                "chrom": fields[0],
                "variant_id": fields[1],
                "pos": position,
                "a1": fields[4],
                "a2": fields[5],
            })
    if not variants:
        raise ValueError(f"No variants parsed from {path}")
    return variants


def in_exclusion(chrom: str, position: int, regions: list[tuple[str, int, int]]) -> bool:
    key = chrom.replace("chr", "")
    return any(key == region[0] and region[1] <= position <= region[2] for region in regions)


def thin_variants(
    variants: list[dict[str, object]],
    thin_bp: int,
    max_per_chrom: int | None,
    regions: list[tuple[str, int, int]],
) -> list[dict[str, object]]:
    selected: list[dict[str, object]] = []
    per_chrom: dict[str, int] = {}
    last_kept: dict[str, int] = {}
    for variant in sorted(variants, key=lambda v: (str(v["chrom"]), int(v["pos"]))):
        chrom, position = str(variant["chrom"]), int(variant["pos"])
        if len(str(variant["a1"])) != 1 or len(str(variant["a2"])) != 1:
            continue
        if in_exclusion(chrom, position, regions):
            continue
        previous = last_kept.get(chrom)
        if previous is not None and position - previous < thin_bp:
            continue
        if max_per_chrom is not None and per_chrom.get(chrom, 0) >= max_per_chrom:
            continue
        selected.append(variant)
        last_kept[chrom] = position
        per_chrom[chrom] = per_chrom.get(chrom, 0) + 1
    return selected


def choose_samples(
    manifest: Path, fam_ids: list[str], max_samples: int | None, seed: int
) -> list[dict[str, str]]:
    fam_index = {sample_id: row for row, sample_id in enumerate(fam_ids)}
    rows: list[dict[str, str]] = []
    with open_text(manifest) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            sample_id = row["sample_id"]
            if sample_id in fam_index:
                rows.append({**row, "fam_row": str(fam_index[sample_id])})
    if not rows:
        raise SystemExit("No manifest sample appears in the .fam file")
    if max_samples is None or len(rows) <= max_samples:
        return sorted(rows, key=lambda item: int(item["fam_row"]))

    by_split: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        by_split.setdefault(row["split"], []).append(row)
    chosen: list[dict[str, str]] = []
    for split, items in sorted(by_split.items()):
        quota = max(1, round(max_samples * len(items) / len(rows)))
        ordered = sorted(items, key=lambda item: stable_hex(seed, split, item["sample_id"]))
        chosen.extend(ordered[:quota])
    return sorted(chosen, key=lambda item: int(item["fam_row"]))


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bed", type=Path, action="append", required=True,
                        help="PLINK .bed path; .bim and .fam must share the stem. Repeatable.")
    parser.add_argument("--split-manifest", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--thin-bp", type=int, default=50_000)
    parser.add_argument("--max-snps-per-chrom", type=int, default=None)
    parser.add_argument("--max-samples", type=int, default=None,
                        help="Subsample participants, preserving split proportions.")
    parser.add_argument("--maf-min", type=float, default=0.01)
    parser.add_argument("--missing-max", type=float, default=0.05)
    parser.add_argument("--train-split-label", default="train",
                        help="Split whose participants fit MAF and missingness.")
    parser.add_argument("--seed", type=int, default=20260917)
    parser.add_argument("--keep-long-range-ld", action="store_true",
                        help="Do not exclude MHC and the default inversion regions.")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    panel_path = out_dir / "panel.int8.npy"
    if panel_path.exists() and not args.overwrite:
        raise SystemExit(f"Refusing to overwrite: {panel_path}")

    beds = [Path(path).resolve() for path in args.bed]
    fam_ids = read_fam_ids(beds[0].with_suffix(".fam"))
    for bed in beds[1:]:
        if read_fam_ids(bed.with_suffix(".fam")) != fam_ids:
            raise SystemExit(f"Sample order differs between {beds[0].name} and {bed.name}")
    total_samples = len(fam_ids)
    bytes_per_variant = (total_samples + 3) // 4

    samples = choose_samples(args.split_manifest, fam_ids, args.max_samples, args.seed)
    rows = np.array([int(row["fam_row"]) for row in samples], dtype=np.int64)
    is_train = np.array(
        [row["split"] == args.train_split_label for row in samples], dtype=bool
    )
    if not is_train.any():
        raise SystemExit(
            f"No {args.train_split_label} participants selected; "
            "MAF cannot be fitted train-only"
        )

    regions = [] if args.keep_long_range_ld else DEFAULT_EXCLUSIONS
    candidates: list[tuple[Path, dict[str, object]]] = []
    for bed in beds:
        variants = read_bim(bed.with_suffix(".bim"))
        for variant in thin_variants(variants, args.thin_bp, args.max_snps_per_chrom, regions):
            candidates.append((bed, variant))
    if not candidates:
        raise SystemExit("Thinning removed every variant; loosen --thin-bp")

    table = dosage_lookup()
    scratch = np.lib.format.open_memmap(
        out_dir / "panel.tmp.npy", mode="w+", dtype=np.int8,
        shape=(len(candidates), len(samples)),
    )

    kept_rows: list[dict[str, object]] = []
    handles: dict[Path, object] = {}
    dropped = {"maf": 0, "missing": 0, "monomorphic": 0}
    try:
        for bed, variant in candidates:
            handle = handles.get(bed)
            if handle is None:
                handle = bed.open("rb")
                if handle.read(3) != b"\x6c\x1b\x01":
                    raise SystemExit(f"{bed.name} is not a SNP-major PLINK 1 .bed")
                handles[bed] = handle
            handle.seek(3 + int(variant["index"]) * bytes_per_variant)
            buffer = handle.read(bytes_per_variant)
            if len(buffer) != bytes_per_variant:
                raise SystemExit(f"Truncated .bed while reading {variant['variant_id']}")
            decoded = table[np.frombuffer(buffer, dtype=np.uint8)].reshape(-1)[:total_samples]
            dosage = decoded[rows]

            train = dosage[is_train]
            observed = train[train >= 0]
            missing_rate = 1.0 - observed.size / max(train.size, 1)
            if missing_rate > args.missing_max:
                dropped["missing"] += 1
                continue
            if observed.size == 0:
                dropped["monomorphic"] += 1
                continue
            a1_frequency = float(observed.sum()) / (2.0 * observed.size)
            maf = min(a1_frequency, 1.0 - a1_frequency)
            if maf < args.maf_min:
                dropped["maf"] += 1
                continue

            scratch[len(kept_rows), :] = dosage
            kept_rows.append({
                **variant,
                "bed": bed.name,
                "train_a1_frequency": round(a1_frequency, 6),
                "train_maf": round(maf, 6),
                "train_missing_rate": round(missing_rate, 6),
            })
    finally:
        for handle in handles.values():
            handle.close()

    if not kept_rows:
        raise SystemExit("Every candidate variant failed QC")

    panel = np.lib.format.open_memmap(
        panel_path, mode="w+", dtype=np.int8, shape=(len(kept_rows), len(samples))
    )
    chunk = max(1, (1 << 22) // max(len(samples), 1))
    for start in range(0, len(kept_rows), chunk):
        stop = min(start + chunk, len(kept_rows))
        panel[start:stop, :] = scratch[start:stop, :]
    panel.flush()
    del panel, scratch
    (out_dir / "panel.tmp.npy").unlink()

    with (out_dir / "panel_variants.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "panel_row", "bed", "bim_index", "chrom", "pos", "variant_id",
            "a1", "a2", "train_a1_frequency", "train_maf", "train_missing_rate",
        ])
        for panel_row, variant in enumerate(kept_rows):
            writer.writerow([
                panel_row, variant["bed"], variant["index"], variant["chrom"], variant["pos"],
                variant["variant_id"], variant["a1"], variant["a2"],
                variant["train_a1_frequency"], variant["train_maf"], variant["train_missing_rate"],
            ])

    with (out_dir / "panel_samples.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["panel_column", "sample_id", "eid", "split", "stratum"])
        for column, row in enumerate(samples):
            writer.writerow([column, row["sample_id"], row.get("eid", ""), row["split"],
                             row.get("stratum", "")])

    summary = {
        "classification": "TASK_B_GENOTYPE_PANEL",
        "identifiers_included": False,
        "dosage_definition": "count of the .bim A1 allele; REF/ALT orientation unverified",
        "layout": "variant_major_int8_missing_as_minus_one",
        "shape": {"variants": len(kept_rows), "samples": len(samples)},
        "fam_sample_count": total_samples,
        "candidates_considered": len(candidates),
        "dropped": dropped,
        "filters": {
            "thin_bp": args.thin_bp,
            "maf_min": args.maf_min,
            "missing_max": args.missing_max,
            "max_snps_per_chrom": args.max_snps_per_chrom,
            "long_range_ld_excluded": not args.keep_long_range_ld,
            "exclusion_regions_grch37": [] if args.keep_long_range_ld else DEFAULT_EXCLUSIONS,
        },
        "split_sample_counts": {
            split: sum(1 for row in samples if row["split"] == split)
            for split in sorted({row["split"] for row in samples})
        },
        "statistics_fitted_on": f"{args.train_split_label} split only",
        "outputs": {
            "panel": str(panel_path),
            "panel_sha256": sha256_file(panel_path),
            "variants": str(out_dir / "panel_variants.tsv"),
            "samples": str(out_dir / "panel_samples.tsv"),
        },
        "test_opened": False,
    }
    (out_dir / "PANEL_SUMMARY.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "status": "OK",
        "variants_kept": len(kept_rows),
        "samples": len(samples),
        "dropped": dropped,
        "panel": str(panel_path),
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
