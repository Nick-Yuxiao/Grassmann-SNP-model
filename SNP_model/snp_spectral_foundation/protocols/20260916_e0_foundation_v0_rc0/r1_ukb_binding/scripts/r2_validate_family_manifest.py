#!/usr/bin/env python3
"""Independently validate a frozen E0 family manifest.

This validator does not trust the builder. It re-derives connected components
directly from the kinship file and compares the induced partition against the
manifest, then checks split disjointness, declared sizes, exclusions and ratios.
It prints counts only, never participant identifiers.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable

ALLOWED_SPLITS = {"train", "validation", "test"}
TARGET = {"train": 0.70, "validation": 0.15, "test": 0.15}


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8-sig", errors="replace", newline="")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


class UnionFind:
    def __init__(self, values: Iterable[str]) -> None:
        self.parent = {value: value for value in values}

    def find(self, value: str) -> str:
        root = value
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[value] != root:
            self.parent[value], value = root, self.parent[value]
        return root

    def union(self, left: str, right: str) -> None:
        root_left, root_right = self.find(left), self.find(right)
        if root_left != root_right:
            self.parent[root_right] = root_left


def read_manifest(path: Path) -> list[dict[str, str]]:
    with open_text(path) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"Family manifest is empty: {path}")
    required = {"sample_id", "component_id", "component_size", "split"}
    missing = required - set(rows[0])
    if missing:
        raise ValueError(f"Family manifest missing columns: {sorted(missing)}")
    return rows


def kinship_edges(path: Path, threshold: float, left: int, right: int, value: int, skip_header: bool):
    with open_text(path) as handle:
        for index, line in enumerate(handle):
            if index == 0 and skip_header:
                continue
            fields = line.split()
            if len(fields) <= max(left, right, value):
                continue
            try:
                coefficient = float(fields[value])
            except ValueError:
                continue
            if coefficient >= threshold and fields[left] != fields[right]:
                yield fields[left], fields[right]


def detect_layout(header: list[str]) -> tuple[int, int, int] | None:
    lowered = [token.lower() for token in header]
    for a, b in (("id1", "id2"), ("iid1", "iid2"), ("sample1", "sample2")):
        if a in lowered and b in lowered:
            for value in ("kinship", "kin", "phi"):
                if value in lowered:
                    return lowered.index(a), lowered.index(b), lowered.index(value)
    return None


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--family-manifest", type=Path, required=True)
    parser.add_argument("--kinship", type=Path, required=True)
    parser.add_argument("--kinship-threshold", type=float, default=0.0442)
    parser.add_argument("--kinship-id-columns", type=int, nargs=2, default=None)
    parser.add_argument("--kinship-value-column", type=int, default=None)
    parser.add_argument("--exclude", type=Path, action="append", default=[])
    parser.add_argument("--ratio-tolerance", type=float, default=0.02)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    rows = read_manifest(args.family_manifest)
    errors: list[str] = []

    sample_ids = [row["sample_id"] for row in rows]
    if len(sample_ids) != len(set(sample_ids)):
        errors.append("duplicate sample_id values in family manifest")

    component_splits: dict[str, set[str]] = defaultdict(set)
    component_members: Counter[str] = Counter()
    declared_sizes: dict[str, set[int]] = defaultdict(set)
    for row in rows:
        if row["split"] not in ALLOWED_SPLITS:
            errors.append(f"invalid split label: {row['split']}")
            continue
        component_splits[row["component_id"]].add(row["split"])
        component_members[row["component_id"]] += 1
        declared_sizes[row["component_id"]].add(int(row["component_size"]))

    cross_split = [key for key, value in component_splits.items() if len(value) != 1]
    if cross_split:
        errors.append(f"components crossing splits: {len(cross_split)}")
    size_mismatch = [
        key for key, count in component_members.items() if declared_sizes[key] != {count}
    ]
    if size_mismatch:
        errors.append(f"component_size mismatches: {len(size_mismatch)}")

    excluded: set[str] = set()
    for path in args.exclude:
        with open_text(path) as handle:
            for line in handle:
                token = line.strip().split()
                if token and not token[0].startswith("#"):
                    excluded.add(token[0])
    leaked_exclusions = excluded & set(sample_ids)
    if leaked_exclusions:
        errors.append(f"excluded identifiers present in manifest: {len(leaked_exclusions)}")

    with open_text(args.kinship) as handle:
        header = handle.readline().split()
    detected = detect_layout(header)
    if args.kinship_id_columns is not None and args.kinship_value_column is not None:
        left, right, value = (
            args.kinship_id_columns[0],
            args.kinship_id_columns[1],
            args.kinship_value_column,
        )
        skip_header = detected is not None
    elif detected is not None:
        left, right, value = detected
        skip_header = True
    else:
        raise SystemExit("Could not detect kinship columns; pass them explicitly.")

    split_of = {row["sample_id"]: row["split"] for row in rows}
    union = UnionFind(sample_ids)
    edges_within = 0
    edges_outside = 0
    cross_split_edges = 0
    for first, second in kinship_edges(args.kinship, args.kinship_threshold, left, right, value, skip_header):
        if first in split_of and second in split_of:
            edges_within += 1
            if split_of[first] != split_of[second]:
                cross_split_edges += 1
            union.union(first, second)
        else:
            edges_outside += 1
    if cross_split_edges:
        errors.append(f"related pairs assigned to different splits: {cross_split_edges}")

    rederived: dict[str, list[str]] = defaultdict(list)
    for sample_id in sample_ids:
        rederived[union.find(sample_id)].append(sample_id)

    manifest_groups: dict[str, list[str]] = defaultdict(list)
    for row in rows:
        manifest_groups[row["component_id"]].append(row["sample_id"])

    manifest_partition = {frozenset(members) for members in manifest_groups.values()}
    rederived_partition = {frozenset(members) for members in rederived.values()}
    identical = manifest_partition == rederived_partition
    if not identical:
        errors.append("manifest components differ from kinship-derived components")
    partition_check = {
        "method": "exact_set_comparison",
        "identical": identical,
        "manifest_component_count": len(manifest_partition),
        "rederived_component_count": len(rederived_partition),
    }

    split_counts = Counter(row["split"] for row in rows)
    total = sum(split_counts.values())
    fractions = {name: split_counts.get(name, 0) / total for name in ALLOWED_SPLITS}
    off_target = {
        name: round(fractions[name] - TARGET[name], 6)
        for name in ALLOWED_SPLITS
        if abs(fractions[name] - TARGET[name]) > args.ratio_tolerance
    }
    if off_target:
        errors.append(f"split fractions outside tolerance: {off_target}")

    report = {
        "classification": "E0_FAMILY_MANIFEST_VALIDATION",
        "status": "PASS" if not errors else "FAIL",
        "identifiers_included": False,
        "family_manifest_sha256": sha256_file(args.family_manifest),
        "kinship_sha256": sha256_file(args.kinship),
        "kinship_threshold": args.kinship_threshold,
        "kinship_columns": {"id": [left, right], "value": value},
        "sample_count": len(rows),
        "component_count": len(component_members),
        "max_component_size": max(component_members.values()) if component_members else 0,
        "component_cross_split_count": len(cross_split),
        "related_pairs_within_manifest": edges_within,
        "related_pairs_outside_manifest": edges_outside,
        "related_pairs_crossing_splits": cross_split_edges,
        "partition_check": partition_check,
        "split_sample_counts": dict(sorted(split_counts.items())),
        "split_fractions": {name: round(value, 6) for name, value in sorted(fractions.items())},
        "excluded_identifiers_present": len(leaked_exclusions),
        "errors": errors,
        "block_axis_frozen": False,
        "run_authorized": False,
    }
    text = json.dumps(report, ensure_ascii=False, indent=2) + "\n"
    if args.output is not None:
        if args.output.exists() and not args.overwrite:
            raise SystemExit(f"Refusing to overwrite: {args.output}")
        args.output.write_text(text, encoding="utf-8")
    print(text, end="")
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
