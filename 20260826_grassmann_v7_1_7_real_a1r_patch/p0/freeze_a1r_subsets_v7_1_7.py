from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def rank_key(seed: int, population: str, sample_id: str) -> str:
    return hashlib.sha256(f"{seed}\t{population}\t{sample_id}".encode()).hexdigest()


def allocate(counts: dict[str, int], target: int) -> dict[str, int]:
    total = sum(counts.values())
    raw = {p: counts[p] * target / total for p in counts}
    quota = {p: math.floor(raw[p]) for p in counts}
    remaining = target - sum(quota.values())
    order = sorted(counts, key=lambda p: (-(raw[p] - quota[p]), p))
    for p in order[:remaining]:
        quota[p] += 1
    if sum(quota.values()) != target or any(quota[p] > counts[p] for p in counts):
        raise ValueError("invalid Hamilton allocation")
    return quota


def write_subset(path: Path, rows: list[dict[str, str]]) -> None:
    fields = ["sample_id", "population", "cohort", "superpopulation"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({field: row[field] for field in fields} for row in rows)
    path.with_suffix(".samples.txt").write_text(
        "".join(f"{row['sample_id']}\n" for row in rows), encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Freeze nested real-1KGP A1-R donor subsets")
    parser.add_argument("--donor-train", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=71107)
    parser.add_argument("--expected-total", type=int, default=2247)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise SystemExit(f"refusing to overwrite output directory: {args.output_dir}")

    with args.donor_train.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample_id", "population", "cohort", "superpopulation"}
        if not reader.fieldnames or not required <= set(reader.fieldnames):
            raise SystemExit(f"missing fields: {sorted(required - set(reader.fieldnames or []))}")
        rows = [dict(row) for row in reader]
    if len(rows) != args.expected_total or len({r["sample_id"] for r in rows}) != len(rows):
        raise SystemExit("donor count or uniqueness mismatch")
    if any(r["cohort"] != "1KGP" or not r["population"] for r in rows):
        raise SystemExit("non-1KGP or unlabelled donor row")

    by_pop: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_pop[row["population"]].append(row)
    ranked = {
        pop: sorted(group, key=lambda r: rank_key(args.seed, pop, r["sample_id"]))
        for pop, group in by_pop.items()
    }
    counts = {pop: len(group) for pop, group in ranked.items()}
    targets = {"25P_562": 562, "50P_1124": 1124, "100P_2247": 2247}
    # Allocate 50% from the full population counts, then allocate 25% from the
    # 50% quotas. This house-monotone construction guarantees nesting while
    # retaining largest-remainder population stratification at both stages.
    quota_50 = allocate(counts, targets["50P_1124"])
    quota_25 = allocate(quota_50, targets["25P_562"])
    quotas = {
        "25P_562": quota_25,
        "50P_1124": quota_50,
        "100P_2247": counts,
    }
    selected: dict[str, list[dict[str, str]]] = {}
    for name, target in targets.items():
        chosen = [row for pop in sorted(ranked) for row in ranked[pop][:quotas[name][pop]]]
        selected[name] = sorted(chosen, key=lambda r: r["sample_id"])
        if len(chosen) != target:
            raise SystemExit(f"target mismatch for {name}")
    ids = {name: {r["sample_id"] for r in group} for name, group in selected.items()}
    if not ids["25P_562"] < ids["50P_1124"] < ids["100P_2247"]:
        raise SystemExit("nested subset invariant failed")

    args.output_dir.mkdir(parents=True, exist_ok=False)
    output_paths: list[Path] = []
    for name, group in selected.items():
        path = args.output_dir / f"DONOR_TRAIN_{name}.tsv"
        write_subset(path, group)
        output_paths.extend([path, path.with_suffix(".samples.txt")])
    audit = {
        "schema_version": "1.0",
        "protocol_version": "v7.1.7",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PASS",
        "seed": args.seed,
        "ranking": "sha256(seed_TAB_population_TAB_sample_id)",
        "allocation": "nested_Hamilton_50_from_100_then_25_from_50_by_population",
        "source_sha256": sha256_file(args.donor_train),
        "counts": {name: len(group) for name, group in selected.items()},
        "nested": True,
        "population_counts": {
            name: dict(sorted(Counter(r["population"] for r in group).items()))
            for name, group in selected.items()
        },
        "sample_file_sha256": {
            name: sha256_file(args.output_dir / f"DONOR_TRAIN_{name}.samples.txt")
            for name in targets
        },
        "hgdp_used": False,
    }
    audit_path = args.output_dir / "A1R_SUBSET_FREEZE.v7.1.7.json"
    audit_path.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    output_paths.append(audit_path)
    manifest = args.output_dir / "A1R_SUBSET_FREEZE.v7.1.7.sha256"
    manifest.write_text(
        "".join(f"{sha256_file(path)}  {path.name}\n" for path in sorted(output_paths)), encoding="utf-8"
    )
    print(json.dumps(audit, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
