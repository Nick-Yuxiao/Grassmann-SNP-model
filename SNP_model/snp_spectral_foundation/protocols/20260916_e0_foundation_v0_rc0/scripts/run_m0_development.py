from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np


CLASSIFICATION = "DEVELOPMENT_ONLY_NON_EVIDENCE"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def gt_to_dosage(gt: tuple[int | None, ...] | None) -> int | None:
    if gt is None or len(gt) != 2 or any(allele is None for allele in gt):
        return None
    if any(allele not in (0, 1) for allele in gt):
        return None
    return int(gt[0]) + int(gt[1])  # type: ignore[arg-type]


def empirical_probabilities(counts: np.ndarray, smoothing: float) -> np.ndarray:
    if counts.shape != (3,) or smoothing <= 0:
        raise ValueError("counts must have shape (3,) and smoothing must be positive")
    return (counts + smoothing) / (counts.sum() + 3.0 * smoothing)


def hwe_probabilities(alt_af: float) -> np.ndarray:
    if not 0.0 <= alt_af <= 1.0:
        raise ValueError("ALT AF must be in [0, 1]")
    return np.asarray(((1.0 - alt_af) ** 2, 2.0 * alt_af * (1.0 - alt_af), alt_af**2))


def stable_rng_seed(mask_seed: int, sample_id: str, block_id: str) -> int:
    payload = f"{mask_seed}\0{sample_id}\0{block_id}".encode("utf-8")
    return int.from_bytes(hashlib.blake2b(payload, digest_size=8).digest(), "little")


def write_tsv_gz(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Development-only M0a/M0b evaluator driven by frozen family/block manifests."
    )
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--family-manifest", type=Path, required=True)
    parser.add_argument("--block-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--evaluation-family-split", default="validation", choices=("validation", "test"))
    parser.add_argument("--evaluation-block-split", default="validation", choices=("validation", "test"))
    parser.add_argument("--mask-rate", type=float, default=0.20)
    parser.add_argument("--mask-seeds", default="20260916")
    parser.add_argument("--smoothing", type=float, default=0.5)
    parser.add_argument("--max-blocks", type=int)
    return parser.parse_args()


def run(args: argparse.Namespace) -> dict[str, object]:
    if not 0.0 < args.mask_rate < 1.0:
        raise ValueError("mask-rate must be between zero and one")
    if args.smoothing <= 0:
        raise ValueError("smoothing must be positive")
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise FileExistsError(f"refusing nonempty output directory: {args.output_dir}")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    family_rows = read_tsv(args.family_manifest)
    block_rows = [row for row in read_tsv(args.block_manifest) if row["split"] == args.evaluation_block_split]
    if args.max_blocks is not None:
        block_rows = block_rows[: args.max_blocks]
    if not block_rows:
        raise ValueError("no evaluation blocks selected")

    train_ids = [row["sample_id"] for row in family_rows if row["split"] == "train"]
    evaluation_rows = [row for row in family_rows if row["split"] == args.evaluation_family_split]
    evaluation_ids = [row["sample_id"] for row in evaluation_rows]
    component_by_sample = {row["sample_id"]: row["component_id"] for row in evaluation_rows}
    if not train_ids or not evaluation_ids:
        raise ValueError("train and evaluation sample sets must both be nonempty")
    seeds = [int(item) for item in args.mask_seeds.split(",") if item.strip()]
    if not seeds:
        raise ValueError("at least one mask seed is required")

    try:
        import pysam
    except ImportError as error:
        raise RuntimeError("pysam==0.24.0 is required") from error

    mask_rows: list[dict[str, object]] = []
    variant_rows: list[dict[str, object]] = []
    cell_sums: dict[tuple[str, str, int], dict[str, float]] = defaultdict(
        lambda: defaultdict(float)
    )
    skipped = defaultdict(int)

    with pysam.VariantFile(str(args.vcf)) as variants:
        available = set(variants.header.samples)
        missing = sorted((set(train_ids) | set(evaluation_ids)) - available)
        if missing:
            raise ValueError(f"manifest samples absent from VCF: {len(missing)}")

        for block in block_rows:
            chrom = block["chrom"]
            start0 = int(block["core_start0"])
            end0 = int(block["core_end0"])
            block_id = block["block_id"]
            rng_by_sample_seed = {
                (sample_id, seed): np.random.default_rng(stable_rng_seed(seed, sample_id, block_id))
                for sample_id in evaluation_ids
                for seed in seeds
            }
            for record in variants.fetch(chrom, start0, end0):
                if record.pos - 1 < start0 or record.pos - 1 >= end0:
                    continue
                if len(record.alleles or ()) != 2 or len(record.ref) != 1 or len(record.alts[0]) != 1:
                    skipped["non_biallelic_snp"] += 1
                    continue
                train_dosages = [gt_to_dosage(record.samples[sample_id].get("GT")) for sample_id in train_ids]
                observed_train = np.asarray([value for value in train_dosages if value is not None], dtype=np.int8)
                if observed_train.size == 0:
                    skipped["no_train_calls"] += 1
                    continue
                counts = np.bincount(observed_train, minlength=3).astype(np.float64)
                empirical = empirical_probabilities(counts, args.smoothing)
                alt_af = float(observed_train.sum() / (2.0 * observed_train.size))
                hwe = hwe_probabilities(alt_af)
                variant_rows.append(
                    {
                        "block_id": block_id,
                        "chrom": chrom,
                        "pos1": record.pos,
                        "ref": record.ref,
                        "alt": record.alts[0],
                        "train_called": int(observed_train.size),
                        "train_count_g0": int(counts[0]),
                        "train_count_g1": int(counts[1]),
                        "train_count_g2": int(counts[2]),
                        "train_alt_af": alt_af,
                        "m0a_p0": empirical[0],
                        "m0a_p1": empirical[1],
                        "m0a_p2": empirical[2],
                        "m0b_p0": hwe[0],
                        "m0b_p1": hwe[1],
                        "m0b_p2": hwe[2],
                    }
                )
                for sample_id in evaluation_ids:
                    target = gt_to_dosage(record.samples[sample_id].get("GT"))
                    if target is None:
                        continue
                    component_id = component_by_sample[sample_id]
                    for seed in seeds:
                        if rng_by_sample_seed[(sample_id, seed)].random() >= args.mask_rate:
                            continue
                        key = (component_id, block_id, seed)
                        stats = cell_sums[key]
                        stats["n_targets"] += 1
                        stats["m0a_ce_sum"] += -math.log(max(float(empirical[target]), 1e-12))
                        stats["m0b_ce_sum"] += -math.log(max(float(hwe[target]), 1e-12))
                        stats["m0a_correct"] += int(int(np.argmax(empirical)) == target)
                        stats["m0b_correct"] += int(int(np.argmax(hwe)) == target)
                        mask_rows.append(
                            {
                                "family_split": args.evaluation_family_split,
                                "block_split": args.evaluation_block_split,
                                "sample_id": sample_id,
                                "component_id": component_id,
                                "block_id": block_id,
                                "chrom": chrom,
                                "pos1": record.pos,
                                "ref": record.ref,
                                "alt": record.alts[0],
                                "mask_seed": seed,
                            }
                        )

    cell_rows: list[dict[str, object]] = []
    for (component_id, block_id, seed), stats in sorted(cell_sums.items()):
        n_targets = int(stats["n_targets"])
        if not n_targets:
            continue
        cell_rows.append(
            {
                "component_id": component_id,
                "block_id": block_id,
                "mask_seed": seed,
                "n_targets": n_targets,
                "m0a_ce": stats["m0a_ce_sum"] / n_targets,
                "m0b_ce": stats["m0b_ce_sum"] / n_targets,
                "m0a_accuracy": stats["m0a_correct"] / n_targets,
                "m0b_accuracy": stats["m0b_correct"] / n_targets,
            }
        )
    if not cell_rows:
        raise ValueError("no masked targets were generated")

    write_tsv_gz(
        args.output_dir / "MASK_MANIFEST.tsv.gz",
        [
            "family_split",
            "block_split",
            "sample_id",
            "component_id",
            "block_id",
            "chrom",
            "pos1",
            "ref",
            "alt",
            "mask_seed",
        ],
        mask_rows,
    )
    write_tsv_gz(
        args.output_dir / "VARIANT_STATS.tsv.gz",
        list(variant_rows[0]),
        variant_rows,
    )
    write_tsv_gz(
        args.output_dir / "M0_CELL_METRICS.tsv.gz",
        list(cell_rows[0]),
        cell_rows,
    )

    summary: dict[str, object] = {
        "classification": CLASSIFICATION,
        "formal_use_allowed": False,
        "evaluation_family_split": args.evaluation_family_split,
        "evaluation_block_split": args.evaluation_block_split,
        "train_sample_count": len(train_ids),
        "evaluation_sample_count": len(evaluation_ids),
        "block_count": len(block_rows),
        "variant_count": len(variant_rows),
        "masked_target_count": len(mask_rows),
        "family_block_seed_cell_count": len(cell_rows),
        "mask_rate": args.mask_rate,
        "mask_seeds": seeds,
        "smoothing": args.smoothing,
        "equal_cell_macro": {
            key: float(np.mean([float(row[key]) for row in cell_rows]))
            for key in ("m0a_ce", "m0b_ce", "m0a_accuracy", "m0b_accuracy")
        },
        "skipped_records": dict(skipped),
        "inputs": {
            "vcf_sha256": sha256(args.vcf),
            "family_manifest_sha256": sha256(args.family_manifest),
            "block_manifest_sha256": sha256(args.block_manifest),
        },
    }
    output_hashes = {}
    for name in ("MASK_MANIFEST.tsv.gz", "VARIANT_STATS.tsv.gz", "M0_CELL_METRICS.tsv.gz"):
        output_hashes[name] = sha256(args.output_dir / name)
    summary["output_sha256"] = output_hashes
    (args.output_dir / "M0_SUMMARY.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return summary


def main() -> None:
    print(json.dumps(run(parse_args()), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
