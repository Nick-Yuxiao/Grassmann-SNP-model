#!/usr/bin/env python3
"""Build deterministic development-only family and genomic-block manifests.

This script is deliberately not a UKB adapter. It exercises the split mechanics on
public genotype data while preserving the formal protocol boundary. Formal runs must
bind an authorized kinship source and an LD/genetic-map-derived block table.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import hashlib
import json
import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


SPLIT_RATIOS = {"train": 0.70, "validation": 0.15, "test": 0.15}


def stable_hex(*parts: object) -> str:
    payload = "\x1f".join(str(part) for part in parts).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


class UnionFind:
    def __init__(self, values: Iterable[str]) -> None:
        self.parent = {value: value for value in values}
        self.rank = {value: 0 for value in values}

    def find(self, value: str) -> str:
        parent = self.parent[value]
        if parent != value:
            self.parent[value] = self.find(parent)
        return self.parent[value]

    def union(self, left: str, right: str) -> None:
        root_left, root_right = self.find(left), self.find(right)
        if root_left == root_right:
            return
        if self.rank[root_left] < self.rank[root_right]:
            root_left, root_right = root_right, root_left
        self.parent[root_right] = root_left
        if self.rank[root_left] == self.rank[root_right]:
            self.rank[root_left] += 1


@dataclass(frozen=True)
class Sample:
    sample_id: str
    population: str
    super_population: str
    sex: str


@dataclass(frozen=True)
class CoreBlock:
    block_id: str
    chrom: str
    start0: int
    end0: int

    def __post_init__(self) -> None:
        if self.start0 < 0 or self.end0 <= self.start0:
            raise ValueError(f"Invalid 0-based half-open block: {self}")


def read_panel(path: Path) -> dict[str, Sample]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows: dict[str, Sample] = {}
        for raw in reader:
            row = {(key or "").strip(): (value or "").strip() for key, value in raw.items()}
            sample_id = row.get("sample") or row.get("Sample")
            if not sample_id:
                continue
            if sample_id in rows:
                raise ValueError(f"Duplicate sample in panel: {sample_id}")
            rows[sample_id] = Sample(
                sample_id=sample_id,
                population=row.get("pop") or row.get("Population") or "UNKNOWN",
                super_population=row.get("super_pop") or row.get("SuperPopulation") or "UNKNOWN",
                sex=row.get("gender") or row.get("Gender") or row.get("sex") or "UNKNOWN",
            )
    if not rows:
        raise ValueError(f"No samples parsed from {path}")
    return rows


def relationship_edges(path: Path, allowed_samples: set[str]) -> list[tuple[str, str]]:
    edges: set[tuple[str, str]] = set()
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for raw in reader:
            row = {(key or "").strip(): (value or "").strip() for key, value in raw.items()}
            sample_id = row.get("Sample") or row.get("sample")
            reason = row.get("Reason for exclusion") or row.get("reason") or ""
            if not sample_id or sample_id not in allowed_samples or ":" not in reason:
                continue
            references = [item.strip() for item in reason.split(":", 1)[1].split(",")]
            for related in references:
                if related and related in allowed_samples and related != sample_id:
                    edges.add(tuple(sorted((sample_id, related))))
    return sorted(edges)


def build_components(sample_ids: Iterable[str], edges: Iterable[tuple[str, str]]) -> dict[str, list[str]]:
    sample_ids = sorted(sample_ids)
    union_find = UnionFind(sample_ids)
    for left, right in edges:
        union_find.union(left, right)
    by_root: dict[str, list[str]] = defaultdict(list)
    for sample_id in sample_ids:
        by_root[union_find.find(sample_id)].append(sample_id)
    components: dict[str, list[str]] = {}
    for members in by_root.values():
        members = sorted(members)
        component_id = "FC_" + stable_hex(*members)[:16]
        components[component_id] = members
    return dict(sorted(components.items()))


def component_stratum(members: list[str], samples: dict[str, Sample]) -> str:
    values = [samples[sample_id].super_population for sample_id in members]
    counts = Counter(values)
    best_count = max(counts.values())
    return sorted(value for value, count in counts.items() if count == best_count)[0]


def assign_component_splits(
    components: dict[str, list[str]], samples: dict[str, Sample], seed: int
) -> dict[str, str]:
    by_stratum: dict[str, list[tuple[str, list[str]]]] = defaultdict(list)
    for component_id, members in components.items():
        by_stratum[component_stratum(members, samples)].append((component_id, members))

    assignments: dict[str, str] = {}
    for stratum, items in sorted(by_stratum.items()):
        total = sum(len(members) for _, members in items)
        target = {name: ratio * total for name, ratio in SPLIT_RATIOS.items()}
        current = {name: 0 for name in SPLIT_RATIOS}
        ordered = sorted(items, key=lambda item: (-len(item[1]), stable_hex(seed, stratum, item[0])))
        for component_id, members in ordered:
            size = len(members)

            def squared_error(candidate: str) -> tuple[float, str]:
                after = dict(current)
                after[candidate] += size
                error = sum((after[name] - target[name]) ** 2 for name in SPLIT_RATIOS)
                return error, stable_hex(seed, stratum, component_id, candidate)

            chosen = min(SPLIT_RATIOS, key=squared_error)
            assignments[component_id] = chosen
            current[chosen] += size
    return assignments


def fixed_bp_blocks(chrom: str, contig_length: int, block_bp: int) -> list[CoreBlock]:
    if contig_length <= 0 or block_bp <= 0:
        raise ValueError("contig_length and block_bp must be positive")
    blocks = []
    for index, start0 in enumerate(range(0, contig_length, block_bp)):
        end0 = min(start0 + block_bp, contig_length)
        blocks.append(CoreBlock(f"DEV_{chrom}_{index:04d}", chrom, start0, end0))
    return blocks


def overlaps(left_start: int, left_end: int, right_start: int, right_end: int) -> bool:
    return left_start < right_end and right_start < left_end


def assign_blocks_with_guards(
    blocks: list[CoreBlock],
    guard_bp: int,
    seed: int,
    *,
    custom_guards: dict[str, tuple[int, int]] | None = None,
    block_source: str = "DEVELOPMENT_ONLY_FIXED_BP",
) -> list[dict[str, object]]:
    if guard_bp < 0:
        raise ValueError("guard_bp must be non-negative")
    blocks = sorted(blocks, key=lambda block: (block.chrom, block.start0, block.end0, block.block_id))
    for previous, current in zip(blocks, blocks[1:]):
        if previous.chrom == current.chrom and previous.end0 > current.start0:
            raise ValueError(f"Overlapping core blocks: {previous.block_id}, {current.block_id}")

    status = {block.block_id: "unassigned" for block in blocks}
    guard_bounds: dict[str, tuple[int, int]] = {}

    def select_holdouts(split: str, ratio: float) -> None:
        target = max(1, round(len(blocks) * ratio))
        candidate_windows: list[list[CoreBlock]] = []
        for start in range(0, len(blocks) - target + 1):
            window = blocks[start : start + target]
            if len({block.chrom for block in window}) != 1:
                continue
            if all(status[block.block_id] == "unassigned" for block in window):
                candidate_windows.append(window)
        if not candidate_windows:
            raise ValueError(f"No contiguous run can hold {target} {split} blocks with guard isolation")
        selected_window = min(
            candidate_windows,
            key=lambda window: stable_hex(seed, split, window[0].block_id, window[-1].block_id),
        )
        for block in selected_window:
            if custom_guards is None:
                guard_start = max(0, block.start0 - guard_bp)
                guard_end = block.end0 + guard_bp
            else:
                guard_start, guard_end = custom_guards[block.block_id]
                if guard_start > block.start0 or guard_end < block.end0:
                    raise ValueError(f"Custom guard does not contain core block {block.block_id}")
            status[block.block_id] = split
            guard_bounds[block.block_id] = (guard_start, guard_end)
        for block in selected_window:
            guard_start, guard_end = guard_bounds[block.block_id]
            for neighbor in blocks:
                if neighbor.chrom != block.chrom or status[neighbor.block_id] != "unassigned":
                    continue
                if overlaps(neighbor.start0, neighbor.end0, guard_start, guard_end):
                    status[neighbor.block_id] = f"buffer_{split}"

    select_holdouts("test", SPLIT_RATIOS["test"])
    select_holdouts("validation", SPLIT_RATIOS["validation"])
    for block in blocks:
        if status[block.block_id] == "unassigned":
            status[block.block_id] = "train"

    rows: list[dict[str, object]] = []
    for block in blocks:
        split = status[block.block_id]
        if split in {"test", "validation"}:
            guard_start0, guard_end0 = guard_bounds[block.block_id]
        else:
            guard_start0, guard_end0 = block.start0, block.end0
        rows.append(
            {
                "block_id": block.block_id,
                "chrom": block.chrom,
                "core_start0": block.start0,
                "core_end0": block.end0,
                "core_start1": block.start0 + 1,
                "core_end1": block.end0,
                "guard_start0": guard_start0,
                "guard_end0": guard_end0,
                "guard_start1": guard_start0 + 1,
                "guard_end1": guard_end0,
                "split": split,
                "coordinate_contract": "start0/end0=0-based-half-open;start1/end1=1-based-inclusive",
                "block_source": block_source,
            }
        )
    validate_block_rows(rows)
    return rows


def read_genetic_map(path: Path, chrom: str) -> tuple[list[int], list[float]]:
    wanted = chrom.removeprefix("chr")
    positions: list[int] = []
    map_cm: list[float] = []
    if path.name.endswith(".gz"):
        handle_context = gzip.open(path, "rt", encoding="utf-8")
    else:
        handle_context = path.open("r", encoding="utf-8")
    with handle_context as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            observed = (row.get("Chromosome") or "").removeprefix("chr")
            if observed != wanted:
                continue
            positions.append(int(row["Position(bp)"]))
            map_cm.append(float(row["Map(cM)"]))
    if len(positions) < 2:
        raise ValueError(f"Fewer than two map rows found for chromosome {chrom}")
    if any(left >= right for left, right in zip(positions, positions[1:])):
        raise ValueError("Genetic-map positions are not strictly increasing")
    if any(left > right for left, right in zip(map_cm, map_cm[1:])):
        raise ValueError("Genetic-map cM values are decreasing")
    return positions, map_cm


def fixed_cm_blocks(
    chrom: str, positions1: list[int], map_cm: list[float], block_cm: float
) -> list[CoreBlock]:
    if block_cm <= 0:
        raise ValueError("block_cm must be positive")
    boundaries1 = [positions1[0]]
    threshold = (math.floor(map_cm[0] / block_cm) + 1) * block_cm
    for position1, value_cm in zip(positions1[1:], map_cm[1:]):
        if value_cm >= threshold:
            if position1 > boundaries1[-1]:
                boundaries1.append(position1)
            while value_cm >= threshold:
                threshold += block_cm
    boundaries1.append(positions1[-1] + 1)
    blocks = []
    for index, (start1, next_start1) in enumerate(zip(boundaries1, boundaries1[1:])):
        start0 = start1 - 1
        end0 = next_start1 - 1
        if end0 > start0:
            blocks.append(CoreBlock(f"DEV_CM_{chrom}_{index:04d}", chrom, start0, end0))
    return blocks


def interpolate_position1_at_cm(
    positions1: list[int], map_cm: list[float], target_cm: float
) -> float:
    if target_cm <= map_cm[0]:
        return float(positions1[0])
    if target_cm >= map_cm[-1]:
        return float(positions1[-1])
    high = bisect.bisect_left(map_cm, target_cm)
    low = high - 1
    low_cm, high_cm = map_cm[low], map_cm[high]
    if high_cm == low_cm:
        return float(positions1[high])
    fraction = (target_cm - low_cm) / (high_cm - low_cm)
    return positions1[low] + fraction * (positions1[high] - positions1[low])


def interpolate_cm_at_position1(
    positions1: list[int], map_cm: list[float], target_position1: int
) -> float:
    if target_position1 <= positions1[0]:
        return map_cm[0]
    if target_position1 >= positions1[-1]:
        return map_cm[-1]
    high = bisect.bisect_left(positions1, target_position1)
    low = high - 1
    fraction = (target_position1 - positions1[low]) / (positions1[high] - positions1[low])
    return map_cm[low] + fraction * (map_cm[high] - map_cm[low])


def cm_guard_bounds(
    blocks: list[CoreBlock],
    positions1: list[int],
    map_cm: list[float],
    guard_cm: float,
) -> dict[str, tuple[int, int]]:
    if guard_cm < 0:
        raise ValueError("guard_cm must be non-negative")
    bounds = {}
    for block in blocks:
        start_cm = interpolate_cm_at_position1(positions1, map_cm, block.start0 + 1)
        end_cm = interpolate_cm_at_position1(positions1, map_cm, block.end0)
        guard_start1 = math.floor(
            interpolate_position1_at_cm(positions1, map_cm, start_cm - guard_cm)
        )
        guard_end1 = math.ceil(
            interpolate_position1_at_cm(positions1, map_cm, end_cm + guard_cm)
        )
        guard_start0 = min(block.start0, max(0, guard_start1 - 1))
        guard_end0 = max(block.end0, guard_end1)
        bounds[block.block_id] = (guard_start0, guard_end0)
    return bounds


def validate_block_rows(rows: list[dict[str, object]]) -> None:
    for row in rows:
        if int(row["core_start1"]) != int(row["core_start0"]) + 1:
            raise ValueError(f"Core coordinate conversion failed for {row['block_id']}")
        if int(row["core_end1"]) != int(row["core_end0"]):
            raise ValueError(f"Core end conversion failed for {row['block_id']}")
        if int(row["guard_start1"]) != int(row["guard_start0"]) + 1:
            raise ValueError(f"Guard coordinate conversion failed for {row['block_id']}")
        if int(row["guard_end1"]) != int(row["guard_end0"]):
            raise ValueError(f"Guard end conversion failed for {row['block_id']}")

    heldout = [row for row in rows if row["split"] in {"validation", "test"}]
    train = [row for row in rows if row["split"] == "train"]
    for train_row in train:
        for heldout_row in heldout:
            if train_row["chrom"] != heldout_row["chrom"]:
                continue
            if overlaps(
                int(train_row["core_start0"]),
                int(train_row["core_end0"]),
                int(heldout_row["guard_start0"]),
                int(heldout_row["guard_end0"]),
            ):
                raise ValueError(
                    f"Train core {train_row['block_id']} overlaps held-out guard {heldout_row['block_id']}"
                )


def read_vcf_header(path: Path, contig: str) -> tuple[list[str], int, dict[str, object]]:
    try:
        import pysam
    except ImportError as error:
        raise RuntimeError("pysam==0.24.0 is required to inspect the VCF header") from error
    with pysam.VariantFile(str(path)) as variants:
        samples = list(variants.header.samples)
        if contig not in variants.header.contigs:
            raise ValueError(f"Contig {contig!r} absent from VCF header")
        length = variants.header.contigs[contig].length
        if length is None:
            raise ValueError(f"Contig {contig!r} has no declared length")
        records_checked = 0
        called_gt = 0
        phased_called_gt = 0
        for record in variants.fetch(contig):
            records_checked += 1
            for call in record.samples.values():
                genotype = call.get("GT")
                if genotype and all(allele is not None for allele in genotype):
                    called_gt += 1
                    phased_called_gt += int(call.phased)
            if records_checked >= 100:
                break
    phase_audit = {
        "records_checked": records_checked,
        "called_gt_checked": called_gt,
        "phased_called_gt": phased_called_gt,
        "phased_called_fraction": phased_called_gt / called_gt if called_gt else None,
    }
    return samples, int(length), phase_audit


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--relationships", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--contig", default="22")
    parser.add_argument("--genome-build", required=True)
    parser.add_argument("--seed", type=int, default=20260916)
    parser.add_argument("--genetic-map", type=Path)
    parser.add_argument("--development-block-cm", type=float, default=0.5)
    parser.add_argument("--development-guard-cm", type=float, default=1.0)
    parser.add_argument("--development-block-bp", type=int, default=1_000_000)
    parser.add_argument("--development-guard-bp", type=int, default=1_000_000)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise SystemExit(f"Refusing non-empty output directory: {args.output_dir}")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    samples = read_panel(args.panel)
    vcf_samples, contig_length, phase_audit = read_vcf_header(args.vcf, args.contig)
    missing_panel = sorted(set(vcf_samples) - set(samples))
    extra_panel = sorted(set(samples) - set(vcf_samples))
    if missing_panel:
        raise SystemExit(f"VCF samples missing from panel: {len(missing_panel)}")
    active_samples = {sample_id: samples[sample_id] for sample_id in vcf_samples}
    edges = relationship_edges(args.relationships, set(active_samples))
    components = build_components(active_samples, edges)
    assignments = assign_component_splits(components, active_samples, args.seed)

    member_to_component = {
        member: component_id for component_id, members in components.items() for member in members
    }
    family_rows = []
    for sample_id in sorted(active_samples):
        sample = active_samples[sample_id]
        component_id = member_to_component[sample_id]
        family_rows.append(
            {
                "sample_id": sample_id,
                "component_id": component_id,
                "component_size": len(components[component_id]),
                "split": assignments[component_id],
                "population": sample.population,
                "super_population": sample.super_population,
                "sex": sample.sex,
                "split_seed": args.seed,
            }
        )

    genetic_map_summary = None
    if args.genetic_map is not None:
        positions1, map_cm = read_genetic_map(args.genetic_map, args.contig)
        blocks = fixed_cm_blocks(args.contig, positions1, map_cm, args.development_block_cm)
        custom_guards = cm_guard_bounds(blocks, positions1, map_cm, args.development_guard_cm)
        block_rows = assign_blocks_with_guards(
            blocks,
            0,
            args.seed,
            custom_guards=custom_guards,
            block_source="DEVELOPMENT_ONLY_FIXED_CM",
        )
        genetic_map_summary = {
            "path_name": args.genetic_map.name,
            "sha256": sha256_file(args.genetic_map),
            "row_count": len(positions1),
            "position1_min": positions1[0],
            "position1_max": positions1[-1],
            "map_cm_min": map_cm[0],
            "map_cm_max": map_cm[-1],
            "development_block_cm": args.development_block_cm,
            "development_guard_cm": args.development_guard_cm,
        }
    else:
        blocks = fixed_bp_blocks(args.contig, contig_length, args.development_block_bp)
        block_rows = assign_blocks_with_guards(blocks, args.development_guard_bp, args.seed)

    write_tsv(
        args.output_dir / "family_manifest.tsv",
        family_rows,
        [
            "sample_id",
            "component_id",
            "component_size",
            "split",
            "population",
            "super_population",
            "sex",
            "split_seed",
        ],
    )
    write_tsv(
        args.output_dir / "block_manifest.tsv",
        block_rows,
        [
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
            "coordinate_contract",
            "block_source",
        ],
    )

    component_splits = defaultdict(set)
    for row in family_rows:
        component_splits[row["component_id"]].add(row["split"])
    overlap_count = sum(len(values) > 1 for values in component_splits.values())
    summary = {
        "classification": "DEVELOPMENT_ONLY_NON_EVIDENCE",
        "genome_build": args.genome_build,
        "contig": args.contig,
        "contig_length": contig_length,
        "vcf_sample_count": len(vcf_samples),
        "vcf_phase_audit": phase_audit,
        "panel_extra_sample_count": len(extra_panel),
        "relationship_edge_count": len(edges),
        "family_component_count": len(components),
        "family_component_overlap_count": overlap_count,
        "sample_split_counts": dict(sorted(Counter(row["split"] for row in family_rows).items())),
        "block_split_counts": dict(sorted(Counter(row["split"] for row in block_rows).items())),
        "genetic_map": genetic_map_summary,
        "coordinate_contract": "VCF POS is 1-based; manifest start0/end0 is 0-based half-open; start1/end1 is 1-based inclusive",
        "formal_use_allowed": False,
        "formal_blocking_reasons": [
            "development blocks are fixed-width engineering units rather than train-derived LD blocks",
            "development guard does not yet take max(receptive field, adjacent LD block, 1 cM)",
            "1000G relationship metadata is not UKB kinship",
        ],
        "inputs": {
            "vcf_sha256": sha256_file(args.vcf),
            "panel_sha256": sha256_file(args.panel),
            "relationships_sha256": sha256_file(args.relationships),
        },
    }
    (args.output_dir / "RUN_SUMMARY.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
