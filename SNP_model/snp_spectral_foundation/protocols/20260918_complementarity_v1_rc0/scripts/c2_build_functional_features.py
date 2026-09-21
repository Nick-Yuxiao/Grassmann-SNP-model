#!/usr/bin/env python3
"""C2: build the functional carrier F from variant annotations and dosage.

A static annotation cannot predict differences between people on its own, so the
annotation is used as a per-variant weight on that person's own dosage:

    F[i, a] = sum_j  standardised_dosage[i, j] * f[j, a]

optionally computed within groups of variants so the carrier keeps some spatial
resolution. Annotation values are phenotype-free by requirement: the manifest must
declare each track phenotype_free, and the script refuses anything else, because a
track derived from trait associations would smuggle the outcome into a predictor.

Only variant positions are used for the join. That matters here: this cohort's
A1/A2 come from a BGEN conversion with no explicit REF/ALT mode, so allele-aware
joins are unreliable while positional ones are not.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from bisect import bisect_right
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_arms import numeric, open_text, read_tsv, sha256_file, standardise_block  # noqa: E402

TRUE_VALUES = {"true", "yes", "1", "t", "y"}


def load_bed(path: Path) -> dict[str, tuple[list[int], list[int], list[float]]]:
    """Read a BED file into per-chromosome sorted starts, ends and scores."""
    staged: dict[str, list[tuple[int, int, float]]] = defaultdict(list)
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split()
            if len(fields) < 3:
                continue
            chrom = fields[0].replace("chr", "")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                continue
            score = 1.0
            if len(fields) >= 5:
                value = numeric(fields[4])
                score = value if np.isfinite(value) else 1.0
            staged[chrom].append((start, end, score))
    indexed: dict[str, tuple[list[int], list[int], list[float]]] = {}
    for chrom, items in staged.items():
        items.sort()
        indexed[chrom] = ([i[0] for i in items], [i[1] for i in items], [i[2] for i in items])
    return indexed


def bed_value(index: dict[str, tuple[list[int], list[int], list[float]]],
              chrom: str, position1: int) -> float:
    """Value of the interval covering a 1-based position, or 0. BED is 0-based half-open."""
    entry = index.get(chrom.replace("chr", ""))
    if entry is None:
        return 0.0
    starts, ends, scores = entry
    position0 = position1 - 1
    slot = bisect_right(starts, position0) - 1
    # Intervals may overlap; walk back a bounded number of them.
    for candidate in range(slot, max(slot - 32, -1), -1):
        if candidate < 0:
            break
        if starts[candidate] <= position0 < ends[candidate]:
            return scores[candidate]
    return 0.0


def load_point_scores(path: Path, value_column: str) -> dict[tuple[str, int], float]:
    values: dict[tuple[str, int], float] = {}
    for row in read_tsv(path):
        chrom = (row.get("chrom") or row.get("chr") or "").replace("chr", "").strip()
        position = row.get("pos") or row.get("position")
        if not chrom or position is None:
            continue
        try:
            key = (chrom, int(position))
        except ValueError:
            continue
        score = numeric(row.get(value_column))
        if np.isfinite(score):
            values[key] = score
    return values


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel-dir", type=Path, required=True)
    parser.add_argument("--split-manifest", type=Path, required=True)
    parser.add_argument("--annotation-manifest", type=Path, required=True,
                        help="TSV with columns track, path, kind (bed|points), value_column, "
                             "phenotype_free, source.")
    parser.add_argument("--group-by", default="chromosome", choices=["none", "chromosome"],
                        help="Compute one aggregate per track, or one per track per chromosome.")
    parser.add_argument("--train-split", default="train")
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    features_path = out_dir / "F_functional.npy"
    if features_path.exists() and not args.overwrite:
        raise SystemExit(f"Refusing to overwrite: {features_path}")

    tracks = list(read_tsv(args.annotation_manifest))
    if not tracks:
        raise SystemExit(f"No tracks listed in {args.annotation_manifest}")
    not_free = [t["track"] for t in tracks
                if (t.get("phenotype_free") or "").strip().lower() not in TRUE_VALUES]
    if not_free:
        raise SystemExit(
            "These tracks are not declared phenotype_free: " + ", ".join(not_free) +
            ". An annotation derived from trait associations would put the outcome into a "
            "predictor, so it cannot be used here."
        )

    panel_dir = args.panel_dir.resolve()
    panel = np.load(panel_dir / "panel.int8.npy", mmap_mode="r")
    variants = list(read_tsv(panel_dir / "panel_variants.tsv"))
    with open_text(panel_dir / "panel_samples.tsv") as handle:
        samples = list(csv.DictReader(handle, delimiter="\t"))
    split_of = {row["sample_id"]: row["split"] for row in read_tsv(args.split_manifest)}
    is_train = np.array([split_of.get(row["sample_id"], "") == args.train_split
                         for row in samples])
    if not is_train.any():
        raise SystemExit(f"No participants in split {args.train_split}")

    weights = np.zeros((len(variants), len(tracks)))
    coverage: list[dict[str, object]] = []
    for column, track in enumerate(tracks):
        path = Path(track["path"])
        if not path.exists():
            raise SystemExit(f"Annotation file missing for {track['track']}: {path}")
        kind = (track.get("kind") or "bed").strip().lower()
        if kind == "bed":
            index = load_bed(path)
            values = np.array([
                bed_value(index, row["chrom"], int(row["pos"])) for row in variants
            ])
        elif kind == "points":
            lookup = load_point_scores(path, (track.get("value_column") or "score").strip())
            values = np.array([
                lookup.get((row["chrom"].replace("chr", ""), int(row["pos"])), 0.0)
                for row in variants
            ])
        else:
            raise SystemExit(f"Unknown kind for {track['track']}: {kind}")
        weights[:, column] = values
        nonzero = int(np.count_nonzero(values))
        coverage.append({
            "track": track["track"],
            "source": track.get("source", ""),
            "kind": kind,
            "variants_annotated": nonzero,
            "coverage": round(nonzero / len(variants), 6),
            "value_mean_over_annotated": round(
                float(values[values != 0].mean()) if nonzero else 0.0, 6
            ),
            "sha256": sha256_file(path),
        })

    usable = [entry for entry in coverage if entry["variants_annotated"] > 0]
    if not usable:
        raise SystemExit(
            "No track annotates a single panel variant. Check the genome build: the panel is "
            "GRCh37 and an annotation on GRCh38 will silently miss everything."
        )

    if args.group_by == "chromosome":
        groups = sorted({row["chrom"] for row in variants}, key=lambda c: int(c))
    else:
        groups = ["ALL"]
    group_of = np.array([
        groups.index(row["chrom"]) if args.group_by == "chromosome" else 0 for row in variants
    ])

    features = np.zeros((len(samples), len(groups) * len(tracks)))
    frequency = np.array([numeric(row["train_a1_frequency"]) for row in variants])
    step = 256
    for start in range(0, len(variants), step):
        stop = min(start + step, len(variants))
        raw = np.asarray(panel[start:stop, :], dtype=np.int8).T
        scaled = standardise_block(raw, frequency[start:stop])
        block_groups = group_of[start:stop]
        block_weights = weights[start:stop]
        for group_index in np.unique(block_groups):
            selector = block_groups == group_index
            contribution = scaled[:, selector] @ block_weights[selector]
            offset = int(group_index) * len(tracks)
            features[:, offset:offset + len(tracks)] += contribution

    # Scale each column on train so the ridge penalty means the same thing per column.
    centre = features[is_train].mean(axis=0)
    scale = features[is_train].std(axis=0)
    scale[scale == 0] = 1.0
    features = ((features - centre) / scale).astype(np.float32)

    keep = features[is_train].std(axis=0) > 0
    column_rows = []
    for column in range(features.shape[1]):
        group_index, track_index = divmod(column, len(tracks))
        column_rows.append({
            "column": column,
            "group": groups[group_index],
            "track": tracks[track_index]["track"],
            "retained": int(bool(keep[column])),
        })
    features = features[:, keep]

    np.save(features_path, features)
    with (out_dir / "F_columns.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(column_rows[0]), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(column_rows)

    summary = {
        "classification": "COMPLEMENTARITY_F_FEATURES",
        "identifiers_included": False,
        "phenotype_read": False,
        "construction": "F[i,a] = sum_j standardised_dosage[i,j] * annotation[j,a]",
        "join": "variant position only; allele orientation is not used",
        "genome_build_expected": "GRCh37",
        "group_by": args.group_by,
        "shape": {"participants": int(features.shape[0]), "features": int(features.shape[1])},
        "tracks": coverage,
        "tracks_with_zero_coverage": [e["track"] for e in coverage if not e["variants_annotated"]],
        "all_tracks_declared_phenotype_free": True,
        "train_split": args.train_split,
        "inputs": {
            "panel_sha256": sha256_file(panel_dir / "panel.int8.npy"),
            "panel_variants_sha256": sha256_file(panel_dir / "panel_variants.tsv"),
            "annotation_manifest_sha256": sha256_file(args.annotation_manifest),
        },
        "outputs": {
            "features": str(features_path),
            "features_sha256": sha256_file(features_path),
        },
    }
    (out_dir / "F_SUMMARY.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "status": "OK",
        "shape": summary["shape"],
        "tracks": [(e["track"], e["coverage"]) for e in coverage],
        "zero_coverage": summary["tracks_with_zero_coverage"],
        "features": str(features_path),
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
