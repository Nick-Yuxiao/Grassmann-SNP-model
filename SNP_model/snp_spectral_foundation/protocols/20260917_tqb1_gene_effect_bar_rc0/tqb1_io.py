"""Input binding and readers for TQ-B1. Plain text plus PLINK 1 .bed only."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from tqb1_core import bed_bytes_per_variant


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def dump_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def read_tsv(path: Path) -> tuple[list[str], list[list[str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        return header, [row for row in reader if row]


def read_fam(path: Path) -> list[str]:
    ids: list[str] = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        for line in handle:
            parts = line.split()
            if not parts:
                continue
            if len(parts) < 6:
                raise ValueError(f".fam line has {len(parts)} fields, expected 6")
            ids.append(parts[1])
    if len(set(ids)) != len(ids):
        raise ValueError(".fam contains duplicate IIDs")
    return ids


def read_bim(path: Path) -> tuple[list[str], np.ndarray, list[str]]:
    chroms: list[str] = []
    positions: list[int] = []
    variant_ids: list[str] = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        for line in handle:
            parts = line.split()
            if not parts:
                continue
            if len(parts) < 6:
                raise ValueError(f".bim line has {len(parts)} fields, expected 6")
            chroms.append(parts[0])
            variant_ids.append(parts[1])
            positions.append(int(parts[3]))
    if len(set(variant_ids)) != len(variant_ids):
        raise ValueError(".bim contains duplicate variant IDs")
    return chroms, np.asarray(positions, dtype=np.int64), variant_ids


def check_bed_size(bed: Path, n_samples: int, n_variants: int) -> None:
    expected = 3 + n_variants * bed_bytes_per_variant(n_samples)
    actual = bed.stat().st_size
    if actual != expected:
        raise ValueError(
            f".bed size {actual} does not match .bim/.fam: expected {expected} bytes for "
            f"{n_variants} variants x {n_samples} samples"
        )


def numeric_table(path: Path, id_column: str = "sample_id") -> tuple[list[str], list[str], np.ndarray]:
    header, rows = read_tsv(path)
    if header[0] != id_column:
        raise ValueError(f"{path.name}: first column must be '{id_column}', found '{header[0]}'")
    names = header[1:]
    ids = [r[0] for r in rows]
    values = np.empty((len(rows), len(names)), dtype=np.float64)
    for i, row in enumerate(rows):
        if len(row) != len(header):
            raise ValueError(f"{path.name}: row {i + 2} has {len(row)} fields, expected {len(header)}")
        for j, cell in enumerate(row[1:]):
            if cell == "" or cell.upper() in {"NA", "NAN"}:
                raise ValueError(f"{path.name}: missing value at row {i + 2}, column '{names[j]}'")
            values[i, j] = float(cell)
    if not np.isfinite(values).all():
        raise ValueError(f"{path.name}: contains non-finite values")
    return ids, names, values


def read_samples(path: Path) -> tuple[list[str], list[str], list[str]]:
    header, rows = read_tsv(path)
    if header[0] != "sample_id":
        raise ValueError("samples file must start with a sample_id column")
    index = {name: i for i, name in enumerate(header)}
    ids = [r[0] for r in rows]
    if len(set(ids)) != len(ids):
        raise ValueError("samples file contains duplicate sample_id values")
    groups = [r[index["group"]] if "group" in index else r[0] for r in rows]
    splits = [r[index["split"]] if "split" in index else "" for r in rows]
    return ids, groups, splits


def read_genes(path: Path) -> list[dict[str, object]]:
    header, rows = read_tsv(path)
    for required in ("gene_id", "chrom", "tss"):
        if required not in header:
            raise ValueError(f"genes file needs a '{required}' column")
    index = {name: i for i, name in enumerate(header)}
    genes = [
        {
            "gene_id": r[index["gene_id"]],
            "chrom": r[index["chrom"]],
            "tss": int(r[index["tss"]]),
        }
        for r in rows
    ]
    ids = [g["gene_id"] for g in genes]
    if len(set(ids)) != len(ids):
        raise ValueError("genes file contains duplicate gene_id values")
    return genes


def read_burden(path: Path) -> dict[tuple[str, str], tuple[float, float]]:
    header, rows = read_tsv(path)
    for required in ("gene_id", "trait", "beta", "se"):
        if required not in header:
            raise ValueError(f"burden file needs a '{required}' column")
    index = {name: i for i, name in enumerate(header)}
    table: dict[tuple[str, str], tuple[float, float]] = {}
    for r in rows:
        key = (r[index["gene_id"]], r[index["trait"]])
        if key in table:
            raise ValueError(f"burden file has duplicate rows for {key}")
        beta = float(r[index["beta"]])
        se = float(r[index["se"]])
        if not (np.isfinite(beta) and np.isfinite(se)) or se <= 0:
            raise ValueError(f"burden file has an invalid beta/se for {key}")
        table[key] = (beta, se)
    return table


def assign_splits(
    ids: list[str], groups: list[str], supplied: list[str], fractions: dict[str, float], seed: int
) -> np.ndarray:
    """Split by group so related individuals never straddle a boundary."""
    if any(s for s in supplied):
        missing = [i for i, s in zip(ids, supplied) if s not in {"train", "validation", "evaluation"}]
        if missing:
            raise ValueError(f"{len(missing)} rows have an unrecognised split label")
        return np.asarray(supplied)
    total = sum(fractions.values())
    if abs(total - 1.0) > 1e-9:
        raise ValueError("split fractions must sum to 1")
    unique = sorted(set(groups))
    keyed = sorted(
        unique, key=lambda g: hashlib.sha256(f"{seed}:{g}".encode()).hexdigest()
    )
    n = len(keyed)
    n_train = int(round(fractions["train"] * n))
    n_validation = int(round(fractions["validation"] * n))
    label_of_group = {}
    for i, group in enumerate(keyed):
        if i < n_train:
            label_of_group[group] = "train"
        elif i < n_train + n_validation:
            label_of_group[group] = "validation"
        else:
            label_of_group[group] = "evaluation"
    return np.asarray([label_of_group[g] for g in groups])
