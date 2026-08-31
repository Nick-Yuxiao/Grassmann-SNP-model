from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


SAMPLE_ID_FIELD = "s"
POPULATION_FIELD = "hgdp_tgp_meta.Population"
SUPERPOPULATION_FIELD = "hgdp_tgp_meta.Genetic.region"


def truth(value: str) -> bool:
    return value.strip().lower() == "true"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_bcf_samples(path: Path) -> list[str]:
    values = [line.strip() for line in path.read_text(encoding="utf-8-sig").splitlines() if line.strip()]
    if not values or len(values) != len(set(values)):
        raise ValueError("BCF sample list is empty or contains duplicate IDs")
    return values


def read_metadata(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            SAMPLE_ID_FIELD,
            POPULATION_FIELD,
            SUPERPOPULATION_FIELD,
            "subsets.tgp",
            "subsets.hgdp",
            "high_quality",
            "release",
            "sample_filters.release_related",
            "sample_filters.all_samples_related",
        }
        if not reader.fieldnames or not required <= set(reader.fieldnames):
            missing = sorted(required - set(reader.fieldnames or []))
            raise ValueError(f"metadata missing required columns: {missing}")
        result: dict[str, dict[str, str]] = {}
        for row in reader:
            sample_id = row[SAMPLE_ID_FIELD].strip()
            if not sample_id or sample_id in result:
                raise ValueError(f"missing or duplicate metadata sample ID: {sample_id!r}")
            result[sample_id] = row
    return result


def cohort(row: dict[str, str]) -> str:
    is_tgp = truth(row["subsets.tgp"])
    is_hgdp = truth(row["subsets.hgdp"])
    if is_tgp == is_hgdp:
        raise ValueError(f"sample {row[SAMPLE_ID_FIELD]} has invalid cohort flags")
    return "1KGP" if is_tgp else "HGDP"


def ranking_key(seed: int, population: str, sample_id: str) -> str:
    value = f"{seed}\t{population}\t{sample_id}".encode("utf-8")
    return hashlib.sha256(value).hexdigest()


def write_samples(path: Path, rows: list[dict[str, str]]) -> None:
    fields = ["sample_id", "population", "cohort", "superpopulation"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({field: row[field] for field in fields} for row in rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Freeze v7.1.1 source-specific sample partitions")
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--bcf-samples", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=71001)
    parser.add_argument("--validation-fraction", type=float, default=0.10)
    parser.add_argument("--expected-tgp", type=int, default=2496)
    parser.add_argument("--expected-hgdp", type=int, default=768)
    args = parser.parse_args()

    if not 0 < args.validation_fraction < 0.5:
        raise SystemExit("validation fraction must be between 0 and 0.5")
    if args.output_dir.exists():
        raise SystemExit(f"refusing to overwrite existing output directory: {args.output_dir}")

    bcf_samples = read_bcf_samples(args.bcf_samples)
    metadata = read_metadata(args.metadata)
    missing = sorted(set(bcf_samples) - set(metadata))
    if missing:
        raise SystemExit(f"BCF samples missing metadata: {missing[:20]}")

    release_rows: list[dict[str, str]] = []
    for sample_id in bcf_samples:
        source = metadata[sample_id]
        if not truth(source["release"]):
            continue
        if not truth(source["high_quality"]):
            raise SystemExit(f"release sample is not high quality: {sample_id}")
        if truth(source["sample_filters.release_related"]):
            raise SystemExit(f"release sample is marked release_related: {sample_id}")
        population = source[POPULATION_FIELD].strip()
        superpopulation = source[SUPERPOPULATION_FIELD].strip()
        if not population or not superpopulation:
            raise SystemExit(f"release sample lacks population metadata: {sample_id}")
        release_rows.append(
            {
                "sample_id": sample_id,
                "population": population,
                "cohort": cohort(source),
                "superpopulation": superpopulation,
                "all_samples_related": str(truth(source["sample_filters.all_samples_related"])).lower(),
            }
        )

    tgp = [row for row in release_rows if row["cohort"] == "1KGP"]
    hgdp = [row for row in release_rows if row["cohort"] == "HGDP"]
    if len(tgp) != args.expected_tgp or len(hgdp) != args.expected_hgdp:
        raise SystemExit(
            f"release count mismatch: observed 1KGP={len(tgp)} HGDP={len(hgdp)}; "
            f"expected 1KGP={args.expected_tgp} HGDP={args.expected_hgdp}"
        )

    by_population: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in tgp:
        by_population[row["population"]].append(row)

    validation_ids: set[str] = set()
    allocation: dict[str, dict[str, int]] = {}
    for population, rows in sorted(by_population.items()):
        if len(rows) < 2:
            raise SystemExit(f"population cannot be represented in both donor partitions: {population}")
        ranked = sorted(rows, key=lambda row: ranking_key(args.seed, population, row["sample_id"]))
        n_validation = math.floor(len(ranked) * args.validation_fraction + 0.5)
        n_validation = min(len(ranked) - 1, max(1, n_validation))
        validation_ids.update(row["sample_id"] for row in ranked[:n_validation])
        allocation[population] = {"all": len(ranked), "train": len(ranked) - n_validation, "validation": n_validation}

    donor_validation = sorted(
        (row for row in tgp if row["sample_id"] in validation_ids), key=lambda row: row["sample_id"]
    )
    donor_train = sorted(
        (row for row in tgp if row["sample_id"] not in validation_ids), key=lambda row: row["sample_id"]
    )
    hgdp_primary = sorted(hgdp, key=lambda row: row["sample_id"])
    if set(row["sample_id"] for row in donor_train) & set(row["sample_id"] for row in donor_validation):
        raise SystemExit("donor train/validation overlap")

    args.output_dir.mkdir(parents=True, exist_ok=False)
    write_samples(args.output_dir / "1KGP_DONOR_TRAIN.tsv", donor_train)
    write_samples(args.output_dir / "1KGP_DONOR_VALIDATION.tsv", donor_validation)
    write_samples(args.output_dir / "HGDP_PRIMARY.tsv", hgdp_primary)
    write_samples(args.output_dir / "HGDP_SNPBAG_CALIBRATION.tsv", [])
    write_samples(args.output_dir / "SAMPLE_POPULATION_COHORT.tsv", sorted(release_rows, key=lambda row: row["sample_id"]))

    calibration = {
        "schema_version": "1.0",
        "protocol_version": "v7.1.1",
        "status": "NOT_COMPARABLE_SOURCE_MISMATCH",
        "current_sample_count": 0,
        "published_expected_sample_count": 216,
        "published_expected_population_count": 54,
        "surrogate_used": False,
        "reason": "audited source has 925 HGDP BCF samples and 52 population labels; exact published IDs/source unavailable",
    }
    (args.output_dir / "CALIBRATION_STATUS.json").write_text(
        json.dumps(calibration, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    summary = {
        "schema_version": "1.0",
        "protocol_version": "v7.1.1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_counts": {"bcf": len(bcf_samples), "1KGP_release": len(tgp), "HGDP_release": len(hgdp)},
        "partitions": {
            "donor_train": len(donor_train),
            "donor_validation": len(donor_validation),
            "hgdp_primary": len(hgdp_primary),
            "hgdp_snpbag_calibration": 0,
        },
        "population_counts": {
            "1KGP": len({row["population"] for row in tgp}),
            "HGDP": len({row["population"] for row in hgdp}),
        },
        "superpopulation_counts": {
            "1KGP": len({row["superpopulation"] for row in tgp}),
            "HGDP": len({row["superpopulation"] for row in hgdp}),
        },
        "donor_allocation_by_population": allocation,
        "release_all_samples_related": {
            "1KGP": sum(row["all_samples_related"] == "true" for row in tgp),
            "HGDP": sum(row["all_samples_related"] == "true" for row in hgdp),
            "HGDP_ids": sorted(row["sample_id"] for row in hgdp if row["all_samples_related"] == "true"),
            "interpretation": "retained representatives; release_related is false",
        },
        "split": {
            "seed": args.seed,
            "validation_fraction": args.validation_fraction,
            "ranking": "sha256(seed_TAB_population_TAB_sample_id)",
        },
        "input_sha256": {
            "metadata": sha256_file(args.metadata),
            "bcf_samples": sha256_file(args.bcf_samples),
        },
        "status": "PASS",
    }
    (args.output_dir / "SPLIT_SUMMARY.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    output_files = sorted(path for path in args.output_dir.iterdir() if path.is_file())
    manifest_lines = [f"{sha256_file(path)}  {path.name}" for path in output_files]
    (args.output_dir / "SAMPLE_FREEZE.sha256").write_text("\n".join(manifest_lines) + "\n", encoding="utf-8")

    print(json.dumps({
        "status": "PASS",
        "donor_train": len(donor_train),
        "donor_validation": len(donor_validation),
        "hgdp_primary": len(hgdp_primary),
        "hgdp_snpbag_calibration": 0,
        "calibration_status": calibration["status"],
        "output_dir": str(args.output_dir),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
