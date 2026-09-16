#!/usr/bin/env python3
"""Freeze the E0 family axis from an authorized kinship file.

Connected components over kinship edges at or above the frozen threshold are the
indivisible split units. Components are assigned to train/validation/test at
70/15/15 inside declared strata with a deterministic, seed-stable rule. The output
manifest schema matches the development manifest so the existing independent
validator keeps working.

No phenotype is read. The summary contains counts only, never identifiers.
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

SPLIT_RATIOS = {"train": 0.70, "validation": 0.15, "test": 0.15}
DEFAULT_THRESHOLD = 0.0442


def stable_hex(*parts: object) -> str:
    payload = "\x1f".join(str(part) for part in parts).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8-sig", errors="replace", newline="")


class UnionFind:
    def __init__(self, values: Iterable[str]) -> None:
        self.parent = {value: value for value in values}
        self.rank = {value: 0 for value in values}

    def find(self, value: str) -> str:
        root = value
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[value] != root:
            self.parent[value], value = root, self.parent[value]
        return root

    def union(self, left: str, right: str) -> None:
        root_left, root_right = self.find(left), self.find(right)
        if root_left == root_right:
            return
        if self.rank[root_left] < self.rank[root_right]:
            root_left, root_right = root_right, root_left
        self.parent[root_right] = root_left
        if self.rank[root_left] == self.rank[root_right]:
            self.rank[root_left] += 1


# ------------------------------------------------------------------- readers


def read_sample_ids(path: Path, column: int | None) -> list[str]:
    suffix = path.suffix.lower()
    ids: list[str] = []
    with open_text(path) as handle:
        for index, line in enumerate(handle):
            stripped = line.strip()
            if not stripped:
                continue
            fields = stripped.split()
            if suffix == ".fam":
                ids.append(fields[1])
            elif suffix == ".psam":
                if stripped.startswith("#"):
                    continue
                ids.append(fields[1] if len(fields) > 1 and fields[0] != fields[1] else fields[0])
            elif suffix == ".sample":
                if index < 2:
                    continue
                ids.append(fields[1] if len(fields) > 1 else fields[0])
            else:
                if index == 0 and column is None and not fields[0].lstrip("-").isdigit():
                    continue
                ids.append(fields[column if column is not None else 0])
    if not ids:
        raise ValueError(f"No sample identifiers parsed from {path}")
    return ids


def read_id_set(paths: Iterable[Path]) -> set[str]:
    values: set[str] = set()
    for path in paths:
        with open_text(path) as handle:
            for line in handle:
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                values.add(stripped.split()[0])
    return values


def detect_kinship_layout(header: list[str]) -> tuple[int, int, int] | None:
    lowered = [token.lower() for token in header]
    for left, right in (("id1", "id2"), ("iid1", "iid2"), ("sample1", "sample2")):
        if left in lowered and right in lowered:
            for value in ("kinship", "kin", "phi"):
                if value in lowered:
                    return lowered.index(left), lowered.index(right), lowered.index(value)
    return None


def read_kinship_edges(
    path: Path,
    threshold: float,
    id_columns: tuple[int, int] | None,
    value_column: int | None,
) -> tuple[list[tuple[str, str]], dict[str, object]]:
    with open_text(path) as handle:
        first = handle.readline()
    header = first.split()
    detected = detect_kinship_layout(header)
    if id_columns is not None and value_column is not None:
        left, right, value_index = id_columns[0], id_columns[1], value_column
        skip_header = detected is not None
        source = "explicit_columns"
    elif detected is not None:
        left, right, value_index = detected
        skip_header = True
        source = "header_detected"
    else:
        raise SystemExit(
            "Could not detect kinship columns. Pass --kinship-id-columns and "
            "--kinship-value-column explicitly."
        )

    edges: list[tuple[str, str]] = []
    rows = 0
    unparsable = 0
    kept = 0
    with open_text(path) as handle:
        for index, line in enumerate(handle):
            if index == 0 and skip_header:
                continue
            fields = line.split()
            if len(fields) <= max(left, right, value_index):
                unparsable += 1
                continue
            try:
                value = float(fields[value_index])
            except ValueError:
                unparsable += 1
                continue
            rows += 1
            if value < threshold:
                continue
            first_id, second_id = fields[left], fields[right]
            if first_id == second_id:
                continue
            kept += 1
            edges.append(tuple(sorted((first_id, second_id))))  # type: ignore[arg-type]
    metadata = {
        "column_source": source,
        "id_columns": [left, right],
        "value_column": value_index,
        "header_fields": header if skip_header else None,
        "pair_rows": rows,
        "unparsable_rows": unparsable,
        "pairs_at_or_above_threshold": kept,
    }
    return sorted(set(edges)), metadata


def read_strata(path: Path, id_column: str, columns: list[str]) -> dict[str, tuple[str, ...]]:
    table: dict[str, tuple[str, ...]] = {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Strata file has no header: {path}")
        missing = [name for name in [id_column, *columns] if name not in reader.fieldnames]
        if missing:
            raise ValueError(f"Strata file missing columns {missing}: {path}")
        for row in reader:
            sample_id = (row.get(id_column) or "").strip()
            if not sample_id:
                continue
            table[sample_id] = tuple((row.get(name) or "NA").strip() or "NA" for name in columns)
    return table


# --------------------------------------------------------------- split logic


def build_components(sample_ids: list[str], edges: list[tuple[str, str]]) -> dict[str, list[str]]:
    union = UnionFind(sample_ids)
    known = set(sample_ids)
    for left, right in edges:
        if left in known and right in known:
            union.union(left, right)
    grouped: dict[str, list[str]] = defaultdict(list)
    for sample_id in sample_ids:
        grouped[union.find(sample_id)].append(sample_id)
    components: dict[str, list[str]] = {}
    for members in grouped.values():
        members = sorted(members)
        components["FC_" + stable_hex(*members)[:16]] = members
    if len(components) != len(grouped):
        raise ValueError("Component identifier collision; widen the component id hash")
    return dict(sorted(components.items()))


def component_stratum(members: list[str], strata: dict[str, tuple[str, ...]]) -> str:
    if not strata:
        return "ALL"
    values = [strata.get(member, ("NA",)) for member in members]
    counts = Counter(values)
    best = max(counts.values())
    winner = sorted(value for value, count in counts.items() if count == best)[0]
    return "|".join(winner)


def assign_component_splits(
    components: dict[str, list[str]], strata_of: dict[str, str], seed: int
) -> dict[str, str]:
    by_stratum: dict[str, list[tuple[str, int]]] = defaultdict(list)
    for component_id, members in components.items():
        by_stratum[strata_of[component_id]].append((component_id, len(members)))

    assignments: dict[str, str] = {}
    for stratum, items in sorted(by_stratum.items()):
        total = sum(size for _, size in items)
        target = {name: ratio * total for name, ratio in SPLIT_RATIOS.items()}
        current = {name: 0 for name in SPLIT_RATIOS}
        ordered = sorted(items, key=lambda item: (-item[1], stable_hex(seed, stratum, item[0])))
        for component_id, size in ordered:

            def squared_error(candidate: str) -> tuple[float, str]:
                after = dict(current)
                after[candidate] += size
                error = sum((after[name] - target[name]) ** 2 for name in SPLIT_RATIOS)
                return error, stable_hex(seed, stratum, component_id, candidate)

            chosen = min(SPLIT_RATIOS, key=squared_error)
            assignments[component_id] = chosen
            current[chosen] += size
    return assignments


# ---------------------------------------------------------------------- main


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-source", type=Path, required=True,
                        help=".fam/.psam/.sample/.txt listing every candidate participant.")
    parser.add_argument("--sample-id-column", type=int, default=None,
                        help="0-based column for plain text sample lists.")
    parser.add_argument("--kinship", type=Path, required=True)
    parser.add_argument("--kinship-threshold", type=float, default=DEFAULT_THRESHOLD)
    parser.add_argument("--kinship-id-columns", type=int, nargs=2, default=None,
                        metavar=("LEFT", "RIGHT"))
    parser.add_argument("--kinship-value-column", type=int, default=None)
    parser.add_argument("--exclude", type=Path, action="append", default=[],
                        help="File of identifiers to drop (withdrawals, QC exclusions); repeatable.")
    parser.add_argument("--strata-file", type=Path, default=None)
    parser.add_argument("--strata-id-column", default="sample_id")
    parser.add_argument("--strata-column", action="append", default=[])
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = out_dir / "family_manifest.tsv"
    summary_path = out_dir / "FAMILY_SUMMARY.json"
    for path in (manifest_path, summary_path):
        if path.exists() and not args.overwrite:
            raise SystemExit(f"Refusing to overwrite: {path}")

    raw_ids = read_sample_ids(args.sample_source, args.sample_id_column)
    duplicate_count = len(raw_ids) - len(set(raw_ids))
    excluded = read_id_set(args.exclude)
    negative_ids = {value for value in raw_ids if value.startswith("-")}
    sample_ids = sorted(set(raw_ids) - excluded - negative_ids)
    if not sample_ids:
        raise SystemExit("No samples remain after exclusions")

    edges, kinship_metadata = read_kinship_edges(
        args.kinship,
        args.kinship_threshold,
        tuple(args.kinship_id_columns) if args.kinship_id_columns else None,
        args.kinship_value_column,
    )
    known = set(sample_ids)
    usable_edges = [edge for edge in edges if edge[0] in known and edge[1] in known]

    strata_table: dict[str, tuple[str, ...]] = {}
    if args.strata_file is not None:
        if not args.strata_column:
            raise SystemExit("--strata-file requires at least one --strata-column")
        strata_table = read_strata(args.strata_file, args.strata_id_column, args.strata_column)

    components = build_components(sample_ids, usable_edges)
    strata_of = {
        component_id: component_stratum(members, strata_table)
        for component_id, members in components.items()
    }
    assignments = assign_component_splits(components, strata_of, args.seed)

    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["sample_id", "component_id", "component_size", "split", "stratum", "split_seed"]
        )
        for component_id, members in components.items():
            split = assignments[component_id]
            for member in members:
                writer.writerow(
                    [member, component_id, len(members), split, strata_of[component_id], args.seed]
                )

    sizes = Counter(len(members) for members in components.values())
    split_samples: Counter[str] = Counter()
    split_components: Counter[str] = Counter()
    per_stratum: dict[str, Counter[str]] = defaultdict(Counter)
    for component_id, members in components.items():
        split = assignments[component_id]
        split_samples[split] += len(members)
        split_components[split] += 1
        per_stratum[strata_of[component_id]][split] += len(members)

    total = sum(split_samples.values())
    summary = {
        "classification": "E0_FAMILY_MANIFEST_SUMMARY",
        "identifiers_included": False,
        "phenotype_read": False,
        "seed": args.seed,
        "kinship": {
            "path": str(args.kinship),
            "sha256": sha256_file(args.kinship),
            "threshold": args.kinship_threshold,
            **kinship_metadata,
            "edges_within_sample_list": len(usable_edges),
            "edges_dropped_outside_sample_list": len(edges) - len(usable_edges),
        },
        "samples": {
            "source_path": str(args.sample_source),
            "source_sha256": sha256_file(args.sample_source),
            "rows_read": len(raw_ids),
            "duplicate_rows": duplicate_count,
            "negative_style_ids_dropped": len(negative_ids),
            "explicit_exclusions_supplied": len(excluded),
            "retained": len(sample_ids),
        },
        "components": {
            "count": len(components),
            "multi_member_count": sum(1 for members in components.values() if len(members) > 1),
            "largest_size": max(sizes) if sizes else 0,
            "size_histogram": dict(sorted(sizes.items())),
            "cross_split_count": 0,
        },
        "split": {
            "sample_counts": dict(sorted(split_samples.items())),
            "component_counts": dict(sorted(split_components.items())),
            "sample_fractions": {
                name: round(count / total, 6) for name, count in sorted(split_samples.items())
            },
            "target_fractions": SPLIT_RATIOS,
            "strata_used": sorted(set(strata_of.values())),
            "per_stratum_sample_counts": {
                stratum: dict(sorted(counts.items())) for stratum, counts in sorted(per_stratum.items())
            },
        },
        "outputs": {
            "family_manifest": str(manifest_path),
            "family_manifest_sha256": sha256_file(manifest_path),
        },
        "formal_use_allowed": False,
        "notes": [
            "Block axis is not frozen here; R2 also requires block_manifest.tsv before RUN_AUTHORIZED.",
            "AF, PCs, standardization and any reference panel must be fitted on train components only.",
        ],
    }
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(
        {
            "status": "OK",
            "family_manifest": str(manifest_path),
            "family_manifest_sha256": summary["outputs"]["family_manifest_sha256"],
            "summary": str(summary_path),
            "samples": len(sample_ids),
            "components": len(components),
            "split_sample_counts": dict(sorted(split_samples.items())),
        },
        ensure_ascii=False,
        indent=2,
    ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
