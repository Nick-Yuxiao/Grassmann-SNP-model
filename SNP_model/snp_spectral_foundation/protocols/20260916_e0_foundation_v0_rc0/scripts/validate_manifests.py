#!/usr/bin/env python3
"""Validate E0 family and block manifests without exposing participant rows."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"Manifest is empty: {path}")
    return rows


def overlaps(left_start: int, left_end: int, right_start: int, right_end: int) -> bool:
    return left_start < right_end and right_start < left_end


def validate_family(rows: list[dict[str, str]]) -> dict[str, object]:
    required = {"sample_id", "component_id", "component_size", "split"}
    missing = required - set(rows[0])
    if missing:
        raise ValueError(f"Family manifest missing columns: {sorted(missing)}")
    sample_ids = [row["sample_id"] for row in rows]
    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError("Family manifest contains duplicate sample_id values")

    component_splits: dict[str, set[str]] = defaultdict(set)
    component_members: Counter[str] = Counter()
    declared_sizes: dict[str, set[int]] = defaultdict(set)
    allowed_splits = {"train", "validation", "test"}
    for row in rows:
        split = row["split"]
        if split not in allowed_splits:
            raise ValueError(f"Invalid family split: {split}")
        component = row["component_id"]
        component_splits[component].add(split)
        component_members[component] += 1
        declared_sizes[component].add(int(row["component_size"]))
    leaking = [component for component, splits in component_splits.items() if len(splits) != 1]
    if leaking:
        raise ValueError(f"Family components cross splits: {len(leaking)}")
    for component, count in component_members.items():
        if declared_sizes[component] != {count}:
            raise ValueError(f"component_size mismatch for {component}")
    return {
        "sample_count": len(rows),
        "component_count": len(component_members),
        "component_overlap_count": 0,
        "sample_split_counts": dict(sorted(Counter(row["split"] for row in rows).items())),
        "max_component_size": max(component_members.values()),
    }


def validate_blocks(rows: list[dict[str, str]]) -> dict[str, object]:
    required = {
        "block_id",
        "chrom",
        "core_start0",
        "core_end0",
        "core_start1",
        "core_end1",
        "guard_start0",
        "guard_end0",
        "guard_start1",
        "guard_end1",
        "split",
    }
    missing = required - set(rows[0])
    if missing:
        raise ValueError(f"Block manifest missing columns: {sorted(missing)}")
    block_ids = [row["block_id"] for row in rows]
    if len(block_ids) != len(set(block_ids)):
        raise ValueError("Block manifest contains duplicate block_id values")

    allowed_splits = {"train", "validation", "test", "buffer_validation", "buffer_test"}
    parsed = []
    for row in rows:
        split = row["split"]
        if split not in allowed_splits:
            raise ValueError(f"Invalid block split: {split}")
        numeric = {
            key: int(row[key])
            for key in [
                "core_start0",
                "core_end0",
                "core_start1",
                "core_end1",
                "guard_start0",
                "guard_end0",
                "guard_start1",
                "guard_end1",
            ]
        }
        if numeric["core_start0"] < 0 or numeric["core_end0"] <= numeric["core_start0"]:
            raise ValueError(f"Invalid core interval for {row['block_id']}")
        if numeric["core_start1"] != numeric["core_start0"] + 1:
            raise ValueError(f"Core start conversion failed for {row['block_id']}")
        if numeric["core_end1"] != numeric["core_end0"]:
            raise ValueError(f"Core end conversion failed for {row['block_id']}")
        if numeric["guard_start1"] != numeric["guard_start0"] + 1:
            raise ValueError(f"Guard start conversion failed for {row['block_id']}")
        if numeric["guard_end1"] != numeric["guard_end0"]:
            raise ValueError(f"Guard end conversion failed for {row['block_id']}")
        if split in {"validation", "test"} and not (
            numeric["guard_start0"] <= numeric["core_start0"]
            and numeric["guard_end0"] >= numeric["core_end0"]
        ):
            raise ValueError(f"Held-out guard does not contain core {row['block_id']}")
        parsed.append({**row, **numeric})

    by_chrom: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in parsed:
        by_chrom[str(row["chrom"])].append(row)
    for chrom_rows in by_chrom.values():
        ordered = sorted(chrom_rows, key=lambda row: (row["core_start0"], row["core_end0"]))
        for previous, current in zip(ordered, ordered[1:]):
            if int(previous["core_end0"]) > int(current["core_start0"]):
                raise ValueError(
                    f"Core blocks overlap: {previous['block_id']} and {current['block_id']}"
                )

    heldout = [row for row in parsed if row["split"] in {"validation", "test"}]
    primary = [row for row in parsed if row["split"] in {"train", "validation", "test"}]
    for core in primary:
        for holdout in heldout:
            if core["block_id"] == holdout["block_id"] or core["chrom"] != holdout["chrom"]:
                continue
            conflict = core["split"] == "train" or core["split"] != holdout["split"]
            if conflict and overlaps(
                int(core["core_start0"]),
                int(core["core_end0"]),
                int(holdout["guard_start0"]),
                int(holdout["guard_end0"]),
            ):
                raise ValueError(
                    f"Core {core['block_id']} ({core['split']}) overlaps "
                    f"guard {holdout['block_id']} ({holdout['split']})"
                )
    return {
        "block_count": len(parsed),
        "block_split_counts": dict(sorted(Counter(row["split"] for row in parsed).items())),
        "chromosome_count": len(by_chrom),
        "coordinate_conversion_errors": 0,
        "cross_split_guard_overlap_count": 0,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--family-manifest", type=Path, required=True)
    parser.add_argument("--block-manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    report = {
        "status": "PASS",
        "identifiers_included": False,
        "family": validate_family(read_tsv(args.family_manifest)),
        "blocks": validate_blocks(read_tsv(args.block_manifest)),
        "manifest_sha256": {
            "family": sha256_file(args.family_manifest),
            "blocks": sha256_file(args.block_manifest),
        },
    }
    payload = json.dumps(report, ensure_ascii=False, indent=2) + "\n"
    if args.output is not None:
        if args.output.exists():
            raise SystemExit(f"Refusing to overwrite: {args.output}")
        args.output.write_text(payload, encoding="utf-8")
    print(payload, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

