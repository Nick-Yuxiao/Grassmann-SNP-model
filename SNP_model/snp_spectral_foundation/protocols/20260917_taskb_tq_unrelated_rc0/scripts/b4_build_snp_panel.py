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
import heapq
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterator

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


def count_bim_rows(path: Path) -> int:
    rows = 0
    with open_text(path) as handle:
        for line in handle:
            if len(line.split()) >= 6:
                rows += 1
    return rows


def stream_bim(path: Path) -> Iterator[dict[str, object]]:
    """Yield variants one at a time.

    An imputed .bim holds millions of rows per chromosome and about 96 million
    across the genome. Materialising them as dictionaries costs tens of gigabytes
    before a single genotype is read, so this streams and lets the caller keep
    only the candidates it can actually use.
    """
    seen = False
    with open_text(path) as handle:
        for index, line in enumerate(handle):
            fields = line.split()
            if len(fields) < 6:
                continue
            try:
                position = int(fields[3])
            except ValueError:
                continue
            seen = True
            yield {
                "index": index,
                "chrom": fields[0],
                "variant_id": fields[1],
                "pos": position,
                "a1": fields[4],
                "a2": fields[5],
            }
    if not seen:
        raise ValueError(f"No variants parsed from {path}")


def in_exclusion(chrom: str, position: int, regions: list[tuple[str, int, int]]) -> bool:
    key = chrom.replace("chr", "")
    return any(key == region[0] and region[1] <= position <= region[2] for region in regions)


def build_windows(
    variants: Iterator[dict[str, object]],
    thin_bp: int,
    max_per_chrom: int | None,
    regions: list[tuple[str, int, int]],
    seed: int,
    max_tries: int,
) -> list[list[dict[str, object]]]:
    """Group biallelic SNPs into fixed-width windows, one kept variant per window.

    Taking the first variant in each window is wrong on an imputed panel: the vast
    majority of imputation v3 sites are rare, so the leftmost variant in a window
    fails the MAF filter roughly nine times out of ten and the window is lost. The
    caller instead walks each window's candidates until one passes QC, so window
    spacing is preserved and the panel actually fills.

    Candidate order inside a window is a seeded hash rather than position, so the
    kept variants are not systematically clustered at window starts. Each window
    retains only the `max_tries` candidates it could ever read, which bounds memory
    at windows x max_tries entries instead of the whole variant table.
    """
    # Per window: a max-heap on the seeded hash, capped at max_tries. Popping the
    # largest keeps exactly the max_tries smallest hashes, the same set a full
    # sort would have taken, so the result does not depend on input order.
    grouped: dict[tuple[str, int], list[tuple[str, int, dict[str, object]]]] = defaultdict(list)
    for variant in variants:
        chrom, position = str(variant["chrom"]), int(variant["pos"])
        if len(str(variant["a1"])) != 1 or len(str(variant["a2"])) != 1:
            continue
        if in_exclusion(chrom, position, regions):
            continue
        window = position // thin_bp
        key = (chrom, window)
        digest = stable_hex(seed, chrom, window, variant["index"])
        heap = grouped[key]
        # heapq is a min-heap, so negate the ordering by pushing the inverted key.
        entry = (_inverted(digest), int(variant["index"]), variant)
        if len(heap) < max_tries:
            heapq.heappush(heap, entry)
        elif entry > heap[0]:
            heapq.heapreplace(heap, entry)

    windows: list[list[dict[str, object]]] = []
    per_chrom: dict[str, int] = {}
    for key in sorted(grouped):
        chrom = key[0]
        if max_per_chrom is not None and per_chrom.get(chrom, 0) >= max_per_chrom:
            continue
        ordered = sorted(grouped[key], key=lambda item: item[0], reverse=True)
        windows.append([item[2] for item in ordered])
        per_chrom[chrom] = per_chrom.get(chrom, 0) + 1
    return windows


def _inverted(digest: str) -> str:
    """Map a hex digest to its order-reversing complement, so a min-heap drops the largest."""
    return "".join("0123456789abcdef"[15 - int(character, 16)] for character in digest)


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
    parser.add_argument("--max-tries-per-window", type=int, default=24,
                        help="Variants read per window before giving that window up. "
                             "Imputed panels are mostly rare variants, so more than one "
                             "candidate per window is normally needed.")
    parser.add_argument("--missing-max", type=float, default=0.05)
    parser.add_argument("--train-split-label", default="train",
                        help="Split whose participants fit MAF and missingness.")
    parser.add_argument("--progress-every", type=int, default=500,
                        help="Write a progress line to stderr every N windows; 0 disables it.")
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

    # Validate every fileset before reading a single genotype. This cohort keeps
    # duplicate copies of some chromosomes and at least one of them is truncated,
    # which otherwise surfaces as a failure an hour into the genotype pass.
    damaged: list[str] = []
    variant_counts: dict[Path, int] = {}
    for bed in beds:
        variant_counts[bed] = count_bim_rows(bed.with_suffix(".bim"))
        expected = 3 + bytes_per_variant * variant_counts[bed]
        actual = bed.stat().st_size
        if actual != expected:
            damaged.append(
                f"{bed}: {actual} bytes, expected {3} + {bytes_per_variant} x "
                f"{variant_counts[bed]} = {expected} (difference {actual - expected})"
            )
    if damaged:
        raise SystemExit(
            "These .bed files do not match their .bim and .fam, so they are incomplete "
            "copies. Point --bed at an intact copy of each:\n  " + "\n  ".join(damaged)
        )
    print(f"[validate] {len(beds)} filesets match their bim and fam", file=sys.stderr, flush=True)

    regions = [] if args.keep_long_range_ld else DEFAULT_EXCLUSIONS
    windows: list[tuple[Path, list[dict[str, object]]]] = []
    for bed in beds:
        for candidates in build_windows(
            stream_bim(bed.with_suffix(".bim")), args.thin_bp, args.max_snps_per_chrom,
            regions, args.seed, args.max_tries_per_window,
        ):
            windows.append((bed, candidates))
        print(f"[windows] {bed.name}: {len(windows)} windows so far",
              file=sys.stderr, flush=True)
    if not windows:
        raise SystemExit("Thinning removed every variant; loosen --thin-bp")

    table = dosage_lookup()
    scratch = np.lib.format.open_memmap(
        out_dir / "panel.tmp.npy", mode="w+", dtype=np.int8,
        shape=(len(windows), len(samples)),
    )

    kept_rows: list[dict[str, object]] = []
    handles: dict[Path, object] = {}
    dropped = {"maf": 0, "missing": 0, "monomorphic": 0}
    variant_reads = 0
    windows_unfilled = 0
    try:
        for window_number, (bed, candidates) in enumerate(windows, start=1):
            handle = handles.get(bed)
            if handle is None:
                handle = bed.open("rb")
                if handle.read(3) != b"\x6c\x1b\x01":
                    raise SystemExit(f"{bed.name} is not a SNP-major PLINK 1 .bed")
                handles[bed] = handle

            if args.progress_every and window_number % args.progress_every == 0:
                print(
                    f"[panel] window {window_number}/{len(windows)} "
                    f"kept={len(kept_rows)} reads={variant_reads}",
                    file=sys.stderr, flush=True,
                )
            filled = False
            for variant in candidates:
                handle.seek(3 + int(variant["index"]) * bytes_per_variant)
                buffer = handle.read(bytes_per_variant)
                if len(buffer) != bytes_per_variant:
                    raise SystemExit(f"Truncated .bed while reading {variant['variant_id']}")
                variant_reads += 1
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
                filled = True
                break
            if not filled:
                windows_unfilled += 1
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
        "filesets": [
            {"bed": str(bed), "variants_in_bim": variant_counts[bed],
             "size_bytes": bed.stat().st_size}
            for bed in beds
        ],
        "windows_considered": len(windows),
        "windows_filled": len(kept_rows),
        "windows_unfilled": windows_unfilled,
        "variant_reads": variant_reads,
        "reads_per_kept_variant": round(variant_reads / max(len(kept_rows), 1), 2),
        "dropped": dropped,
        "filters": {
            "thin_bp": args.thin_bp,
            "maf_min": args.maf_min,
            "missing_max": args.missing_max,
            "max_snps_per_chrom": args.max_snps_per_chrom,
            "max_tries_per_window": args.max_tries_per_window,
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
        "windows_considered": len(windows),
        "windows_unfilled": windows_unfilled,
        "reads_per_kept_variant": round(variant_reads / max(len(kept_rows), 1), 2),
        "dropped": dropped,
        "panel": str(panel_path),
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
