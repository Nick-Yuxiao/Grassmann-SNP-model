from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_tsv_gz(path: Path):
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Independently validate development M0 artifacts.")
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--family-manifest", type=Path, required=True)
    parser.add_argument("--block-manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def validate(args: argparse.Namespace) -> dict[str, object]:
    summary = json.loads((args.run_dir / "M0_SUMMARY.json").read_text(encoding="utf-8"))
    if summary.get("classification") != "DEVELOPMENT_ONLY_NON_EVIDENCE":
        raise ValueError("unexpected classification")
    if summary.get("formal_use_allowed") is not False:
        raise ValueError("development result must prohibit formal use")

    output_hashes = summary.get("output_sha256", {})
    for name in ("MASK_MANIFEST.tsv.gz", "VARIANT_STATS.tsv.gz", "M0_CELL_METRICS.tsv.gz"):
        observed = sha256(args.run_dir / name)
        if output_hashes.get(name) != observed:
            raise ValueError(f"hash mismatch for {name}")

    family_rows = read_tsv(args.family_manifest)
    evaluation_split = summary["evaluation_family_split"]
    evaluation_rows = [row for row in family_rows if row["split"] == evaluation_split]
    component_by_sample = {row["sample_id"]: row["component_id"] for row in evaluation_rows}
    evaluation_components = set(component_by_sample.values())

    block_rows = read_tsv(args.block_manifest)
    block_split = summary["evaluation_block_split"]
    blocks = {
        row["block_id"]: (row["chrom"], int(row["core_start1"]), int(row["core_end1"]))
        for row in block_rows
        if row["split"] == block_split
    }

    variant_keys: set[tuple[str, str, int, str, str]] = set()
    variant_count = 0
    probability_errors = 0
    for row in read_tsv_gz(args.run_dir / "VARIANT_STATS.tsv.gz"):
        key = (row["block_id"], row["chrom"], int(row["pos1"]), row["ref"], row["alt"])
        if key in variant_keys:
            raise ValueError("duplicate variant-stat key")
        variant_keys.add(key)
        if row["block_id"] not in blocks:
            raise ValueError("variant uses a non-evaluation block")
        chrom, start1, end1 = blocks[row["block_id"]]
        if row["chrom"] != chrom or not start1 <= int(row["pos1"]) <= end1:
            raise ValueError("variant lies outside its core block")
        for prefix in ("m0a", "m0b"):
            values = [float(row[f"{prefix}_p{genotype}"]) for genotype in range(3)]
            if any(value < 0.0 or value > 1.0 for value in values) or not math.isclose(
                sum(values), 1.0, rel_tol=1e-9, abs_tol=1e-9
            ):
                probability_errors += 1
        variant_count += 1
    if probability_errors:
        raise ValueError(f"invalid probability rows: {probability_errors}")

    mask_keys: set[tuple[str, str, int, str, str, int]] = set()
    mask_count = 0
    for row in read_tsv_gz(args.run_dir / "MASK_MANIFEST.tsv.gz"):
        if row["family_split"] != evaluation_split or row["block_split"] != block_split:
            raise ValueError("mask split label mismatch")
        if component_by_sample.get(row["sample_id"]) != row["component_id"]:
            raise ValueError("mask sample/component mismatch")
        variant_key = (row["block_id"], row["chrom"], int(row["pos1"]), row["ref"], row["alt"])
        if variant_key not in variant_keys:
            raise ValueError("mask target lacks frozen variant stats")
        key = (
            row["sample_id"],
            row["block_id"],
            int(row["pos1"]),
            row["ref"],
            row["alt"],
            int(row["mask_seed"]),
        )
        if key in mask_keys:
            raise ValueError("duplicate mask key")
        mask_keys.add(key)
        mask_count += 1

    cell_keys: set[tuple[str, str, int]] = set()
    cell_count = 0
    cell_target_count = 0
    for row in read_tsv_gz(args.run_dir / "M0_CELL_METRICS.tsv.gz"):
        key = (row["component_id"], row["block_id"], int(row["mask_seed"]))
        if key in cell_keys:
            raise ValueError("duplicate cell metric key")
        cell_keys.add(key)
        if row["component_id"] not in evaluation_components or row["block_id"] not in blocks:
            raise ValueError("cell metric crosses the evaluation split")
        n_targets = int(row["n_targets"])
        if n_targets <= 0:
            raise ValueError("cell metric has no targets")
        for field in ("m0a_ce", "m0b_ce"):
            if not math.isfinite(float(row[field])) or float(row[field]) < 0:
                raise ValueError("invalid CE")
        for field in ("m0a_accuracy", "m0b_accuracy"):
            if not 0.0 <= float(row[field]) <= 1.0:
                raise ValueError("invalid accuracy")
        cell_target_count += n_targets
        cell_count += 1

    expected = {
        "variant_count": variant_count,
        "masked_target_count": mask_count,
        "family_block_seed_cell_count": cell_count,
    }
    for field, observed in expected.items():
        if int(summary[field]) != observed:
            raise ValueError(f"summary count mismatch for {field}")
    if cell_target_count != mask_count:
        raise ValueError("cell target counts do not sum to mask count")

    result: dict[str, object] = {
        "status": "PASS",
        "classification": "DEVELOPMENT_ONLY_NON_EVIDENCE",
        "formal_use_allowed": False,
        "identifiers_included": False,
        "checks": {
            "output_hashes_match": True,
            "evaluation_family_split_only": True,
            "evaluation_block_split_only": True,
            "mask_keys_unique": True,
            "variant_probabilities_normalized": True,
            "core_coordinate_membership": True,
            "cell_target_sum_matches_masks": True,
        },
        "counts": expected,
        "input_manifest_sha256": {
            "family": sha256(args.family_manifest),
            "blocks": sha256(args.block_manifest),
        },
    }
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return result


def main() -> None:
    print(json.dumps(validate(parse_args()), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
