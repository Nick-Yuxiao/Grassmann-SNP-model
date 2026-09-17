#!/usr/bin/env python3
"""Shared deterministic split and IO helpers for the Task B qualification package.

The stratified assignment is byte-for-byte the same rule used by the E0 R2 family
manifest builder, so a split produced here and a split produced there agree given
the same components, strata and seed.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Iterator

SPLIT_RATIOS = {"train": 0.70, "validation": 0.15, "test": 0.15}


def stable_hex(*parts: object) -> str:
    payload = "\x1f".join(str(part) for part in parts).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def open_text(path: Path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8-sig", errors="replace", newline="")


def read_tsv(path: Path) -> Iterator[dict[str, str]]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"No header in {path}")
        for row in reader:
            yield row


def read_id_set(paths: Iterable[Path]) -> set[str]:
    values: set[str] = set()
    for path in paths:
        with open_text(Path(path)) as handle:
            for line in handle:
                token = line.strip().split()
                if token and not token[0].startswith("#"):
                    values.add(token[0])
    return values


def read_fam_ids(path: Path) -> list[str]:
    ids: list[str] = []
    with open_text(Path(path)) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) >= 2:
                ids.append(fields[1])
    if not ids:
        raise ValueError(f"No sample identifiers parsed from {path}")
    return ids


def assign_splits(units: dict[str, int], strata_of: dict[str, str], seed: int) -> dict[str, str]:
    """Assign indivisible units to train/validation/test inside each stratum.

    `units` maps unit id to its size in individuals. The rule is greedy in a
    deterministic order: larger units first, ties broken by a seeded hash, each
    unit going to whichever split minimises the squared deviation from the target
    counts. Identical inputs always give identical output.
    """
    by_stratum: dict[str, list[tuple[str, int]]] = defaultdict(list)
    for unit_id, size in units.items():
        by_stratum[strata_of[unit_id]].append((unit_id, size))

    assignments: dict[str, str] = {}
    for stratum, items in sorted(by_stratum.items()):
        total = sum(size for _, size in items)
        target = {name: ratio * total for name, ratio in SPLIT_RATIOS.items()}
        current = {name: 0 for name in SPLIT_RATIOS}
        ordered = sorted(items, key=lambda item: (-item[1], stable_hex(seed, stratum, item[0])))
        for unit_id, size in ordered:

            def squared_error(candidate: str) -> tuple[float, str]:
                after = dict(current)
                after[candidate] += size
                error = sum((after[name] - target[name]) ** 2 for name in SPLIT_RATIOS)
                return error, stable_hex(seed, stratum, unit_id, candidate)

            chosen = min(SPLIT_RATIOS, key=squared_error)
            assignments[unit_id] = chosen
            current[chosen] += size
    return assignments


def split_counts(assignments: dict[str, str], units: dict[str, int]) -> dict[str, int]:
    counts: Counter[str] = Counter()
    for unit_id, split in assignments.items():
        counts[split] += units[unit_id]
    return dict(sorted(counts.items()))
