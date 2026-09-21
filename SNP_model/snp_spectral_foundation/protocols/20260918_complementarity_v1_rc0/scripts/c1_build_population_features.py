#!/usr/bin/env python3
"""C1: build the population-structure carrier P from genotype alone.

P is deliberately the cheapest thing that can carry population-genetic structure:
principal components taken within contiguous blocks of the frozen panel, fitted on
training participants only and then projected onto everyone.

Why blocks rather than a global PCA: arm A already contains UKB's genome-wide
PC1-40. A global PCA of the same panel would be close to collinear with those, and
P would have almost nothing left to add by construction. Block-local components
instead carry local haplotype structure, which the global ancestry axes do not
represent. The summary reports the canonical correlations between P and the global
PCs so that overlap is measured rather than assumed.

No phenotype is read here.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_arms import (  # noqa: E402
    canonical_correlations,
    numeric,
    open_text,
    read_tsv,
    sha256_file,
    standardise_block,
)


def load_panel_samples(panel_dir: Path) -> list[dict[str, str]]:
    with open_text(panel_dir / "panel_samples.tsv") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_panel_variants(panel_dir: Path) -> list[dict[str, str]]:
    return list(read_tsv(panel_dir / "panel_variants.tsv"))


def block_ranges(variants: list[dict[str, str]], block_size: int) -> list[tuple[int, int, str]]:
    """Contiguous runs of panel rows, never spanning a chromosome boundary."""
    ranges: list[tuple[int, int, str]] = []
    start = 0
    for index in range(1, len(variants) + 1):
        crossed = index == len(variants) or variants[index]["chrom"] != variants[start]["chrom"]
        if index - start == block_size or crossed:
            ranges.append((start, index, variants[start]["chrom"]))
            start = index
    return [item for item in ranges if item[1] > item[0]]


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel-dir", type=Path, required=True)
    parser.add_argument("--split-manifest", type=Path, required=True)
    parser.add_argument("--covariates", type=Path, default=None,
                        help="Covariates TSV, used only to measure overlap with the global PCs.")
    parser.add_argument("--covariate-id-column", default="eid")
    parser.add_argument("--global-pc-prefix", default="n_22009_0_",
                        help="Column prefix of the cohort's genome-wide PCs.")
    parser.add_argument("--global-pc-count", type=int, default=40)
    parser.add_argument("--block-size", type=int, default=100,
                        help="Panel variants per block.")
    parser.add_argument("--pcs-per-block", type=int, default=3)
    parser.add_argument("--train-split", default="train")
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    features_path = out_dir / "P_population.npy"
    if features_path.exists() and not args.overwrite:
        raise SystemExit(f"Refusing to overwrite: {features_path}")

    panel_dir = args.panel_dir.resolve()
    panel = np.load(panel_dir / "panel.int8.npy", mmap_mode="r")
    samples = load_panel_samples(panel_dir)
    variants = load_panel_variants(panel_dir)
    if panel.shape != (len(variants), len(samples)):
        raise SystemExit("panel matrix does not match its variant and sample tables")

    split_of = {row["sample_id"]: row["split"] for row in read_tsv(args.split_manifest)}
    splits = np.array([split_of.get(row["sample_id"], "") for row in samples])
    is_train = splits == args.train_split
    if not is_train.any():
        raise SystemExit(f"No participants in split {args.train_split}")

    frequency = np.array([numeric(row["train_a1_frequency"]) for row in variants])
    ranges = block_ranges(variants, args.block_size)
    columns_per_block = args.pcs_per_block
    features = np.zeros((len(samples), len(ranges) * columns_per_block), dtype=np.float32)

    column_rows: list[dict[str, object]] = []
    explained: list[float] = []
    for block_index, (start, stop, chrom) in enumerate(ranges):
        raw = np.asarray(panel[start:stop, :], dtype=np.int8).T
        scaled = standardise_block(raw, frequency[start:stop])
        train_block = scaled[is_train]
        covariance = (train_block.T @ train_block) / max(train_block.shape[0] - 1, 1)
        eigenvalues, eigenvectors = np.linalg.eigh(covariance)
        order = np.argsort(eigenvalues)[::-1][:columns_per_block]
        total_variance = float(np.sum(np.maximum(eigenvalues, 0.0)))
        for rank, position in enumerate(order):
            vector = eigenvectors[:, position]
            # Sign of an eigenvector is arbitrary; fix it so reruns agree.
            if vector[np.argmax(np.abs(vector))] < 0:
                vector = -vector
            projected = scaled @ vector
            centre = float(projected[is_train].mean())
            scale = float(projected[is_train].std()) or 1.0
            features[:, block_index * columns_per_block + rank] = (projected - centre) / scale
            share = float(max(eigenvalues[position], 0.0)) / total_variance if total_variance else 0.0
            explained.append(share)
            column_rows.append({
                "column": block_index * columns_per_block + rank,
                "block": block_index,
                "chrom": chrom,
                "panel_row_start": start,
                "panel_row_stop": stop,
                "variants_in_block": stop - start,
                "pc_rank": rank + 1,
                "train_variance_share": round(share, 6),
            })

    np.save(features_path, features)
    with (out_dir / "P_columns.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(column_rows[0]), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(column_rows)
    with (out_dir / "P_samples.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["row", "sample_id", "eid", "split"])
        for row_index, row in enumerate(samples):
            writer.writerow([row_index, row["sample_id"], row.get("eid", ""),
                             split_of.get(row["sample_id"], "")])

    overlap: dict[str, object] = {"measured": False}
    if args.covariates is not None:
        wanted = [f"{args.global_pc_prefix}{i}" for i in range(1, args.global_pc_count + 1)]
        table = {}
        for row in read_tsv(args.covariates):
            key = (row.get(args.covariate_id_column) or "").strip()
            if key:
                table[key] = row
        keys = [row.get("eid", "") for row in samples]
        present = [name for name in wanted if name in next(iter(table.values()), {})]
        if present:
            global_pcs = np.array([
                [numeric(table.get(key, {}).get(name)) for name in present] for key in keys
            ])
            finite = np.isfinite(global_pcs).all(axis=1) & is_train
            if finite.sum() > len(present):
                sample_rows = np.flatnonzero(finite)[:20000]
                overlap = {
                    "measured": True,
                    "global_pc_columns": len(present),
                    "rows_used": int(sample_rows.size),
                    "top_canonical_correlations": canonical_correlations(
                        features[sample_rows], global_pcs[sample_rows], top=5
                    ),
                    "reading": (
                        "Canonical correlations near 1 mean P largely repeats the global "
                        "ancestry axes arm A already carries, and P would add little by "
                        "construction. Values well below 1 mean P carries local structure "
                        "those axes do not represent."
                    ),
                }

    summary = {
        "classification": "COMPLEMENTARITY_P_FEATURES",
        "identifiers_included": False,
        "phenotype_read": False,
        "design": "block-local PCA on the frozen panel, fitted on train only",
        "why_not_global_pca": (
            "arm A already contains the cohort's genome-wide PC1-40, so a global PCA of the "
            "same panel would be near-collinear with it"
        ),
        "shape": {"participants": int(features.shape[0]), "features": int(features.shape[1])},
        "blocks": len(ranges),
        "block_size_variants": args.block_size,
        "pcs_per_block": args.pcs_per_block,
        "median_block_variance_share": round(float(np.median(explained)), 6) if explained else None,
        "standardisation": "(dosage - 2p)/sqrt(2p(1-p)) with train-only p; missing set to 0",
        "train_split": args.train_split,
        "train_participants": int(is_train.sum()),
        "overlap_with_global_pcs": overlap,
        "inputs": {
            "panel_sha256": sha256_file(panel_dir / "panel.int8.npy"),
            "panel_variants_sha256": sha256_file(panel_dir / "panel_variants.tsv"),
            "split_manifest_sha256": sha256_file(args.split_manifest),
        },
        "outputs": {
            "features": str(features_path),
            "features_sha256": sha256_file(features_path),
        },
    }
    (out_dir / "P_SUMMARY.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "status": "OK",
        "shape": summary["shape"],
        "blocks": summary["blocks"],
        "overlap_with_global_pcs": overlap.get("top_canonical_correlations"),
        "features": str(features_path),
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
