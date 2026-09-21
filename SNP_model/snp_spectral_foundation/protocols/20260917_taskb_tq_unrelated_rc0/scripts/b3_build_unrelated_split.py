#!/usr/bin/env python3
"""B3: freeze a relatedness-disjoint split on the unrelated subset.

Every retained participant has UKB field 22021 == 0, meaning no third-degree or
closer relative was identified among participants. Each retained participant is
therefore a singleton relatedness component, and an individual-level split is
relatedness-disjoint by construction.

What this buys and what it costs is stated in the summary: the design answers
"a new unrelated individual", not "a new family". It is not a substitute for the
pairwise kinship file when a family-level claim is wanted.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_split import (  # noqa: E402
    BRIDGE_RESERVED_RATIOS,
    SPLIT_RATIOS,
    assign_splits,
    parse_ratio_arguments,
    open_text,
    read_fam_ids,
    read_id_set,
    read_tsv,
    sha256_file,
)

DEFAULT_KINSHIP_COLUMN = "n_22021_0_0"
DEFAULT_ID_COLUMN = "eid"

# UKB data-coding for field 22021. Verify against the showcase before freezing;
# the observed histogram in the summary is what actually justifies the filter.
FLAG_CODING_REFERENCE = {
    "0": "no kinship found -> keep",
    "1": "at least one relative identified -> drop",
    "10": "ten or more third-degree relatives identified -> drop",
    "-1": "participant excluded from the kinship inference process -> drop, relatedness UNKNOWN",
    "": "value absent -> drop, relatedness UNKNOWN",
}


def array_from_batch(value: str) -> str:
    try:
        batch = int(float(value))
    except (TypeError, ValueError):
        return "NA"
    if batch < 0:
        return "UKBiLEVE"
    return "Axiom" if batch > 0 else "NA"


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--covariates", type=Path, required=True,
                        help="TSV with the identifier column plus strata columns.")
    parser.add_argument("--covariate-id-column", default="eid")
    parser.add_argument("--genotype-id-column", default="iid",
                        help="Column in the covariates TSV holding the genotype fam IID.")
    parser.add_argument("--fam", type=Path, required=True,
                        help="Genotype .fam file; participants absent from it are dropped.")
    parser.add_argument("--kinship-flag", type=Path, required=True,
                        help="TSV with the identifier column and UKB field 22021.")
    parser.add_argument("--kinship-flag-id-column", default=DEFAULT_ID_COLUMN)
    parser.add_argument("--kinship-flag-column", default=DEFAULT_KINSHIP_COLUMN)
    parser.add_argument("--keep-flag-value", action="append", default=["0"],
                        help="Values of field 22021 kept as unrelated; default keeps only 0.")
    parser.add_argument("--qc-flag", action="append", default=[],
                        metavar="COLUMN",
                        help="Column that must be blank/zero to keep the participant "
                             "(e.g. n_22027_0_0 heterozygosity outlier); repeatable.")
    parser.add_argument("--exclude", type=Path, action="append", default=[],
                        help="File of identifiers to drop (withdrawals); repeatable.")
    parser.add_argument("--strata-column", action="append", default=[],
                        help="Covariate column used for stratification; repeatable.")
    parser.add_argument("--batch-column", default="n_22000_0_0",
                        help="Batch column turned into an array label and used as a stratum.")
    parser.add_argument("--split", action="append", default=[], metavar="NAME=FRACTION",
                        help="Split name and fraction; repeatable. Fractions must sum to 1.")
    parser.add_argument("--reserve-bridge-holdout", action="store_true",
                        help="Shorthand for train=0.55 validation=0.15 tq_test=0.15 "
                             "bridge_holdout=0.15, keeping a Bridge test the qualification "
                             "never opens.")
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = out_dir / "tq_split_manifest.tsv"
    summary_path = out_dir / "SPLIT_SUMMARY.json"
    for path in (manifest_path, summary_path):
        if path.exists() and not args.overwrite:
            raise SystemExit(f"Refusing to overwrite: {path}")

    if args.split and args.reserve_bridge_holdout:
        raise SystemExit("Pass either --split or --reserve-bridge-holdout, not both")
    if args.split:
        ratios = parse_ratio_arguments(args.split)
    elif args.reserve_bridge_holdout:
        ratios = dict(BRIDGE_RESERVED_RATIOS)
    else:
        ratios = dict(SPLIT_RATIOS)

    keep_values = {value.strip() for value in args.keep_flag_value}
    flag_of: dict[str, str] = {}
    qc_of: dict[str, dict[str, str]] = {}
    for row in read_tsv(args.kinship_flag):
        identifier = (row.get(args.kinship_flag_id_column) or "").strip()
        if not identifier:
            continue
        raw = (row.get(args.kinship_flag_column) or "").strip()
        try:
            normalized = str(int(float(raw))) if raw else ""
        except ValueError:
            normalized = raw
        flag_of[identifier] = normalized
        if args.qc_flag:
            qc_of[identifier] = {name: (row.get(name) or "").strip() for name in args.qc_flag}
    if not flag_of:
        raise SystemExit(f"No rows parsed from {args.kinship_flag}")

    fam_ids = set(read_fam_ids(args.fam))
    excluded = read_id_set(args.exclude)

    counters = Counter()
    flag_histogram: Counter[str] = Counter()
    dropped_by_flag_value: Counter[str] = Counter()
    kept: dict[str, dict[str, str]] = {}
    for row in read_tsv(args.covariates):
        counters["covariate_rows"] += 1
        identifier = (row.get(args.covariate_id_column) or "").strip()
        genotype_id = (row.get(args.genotype_id_column) or "").strip()
        if not identifier or not genotype_id:
            counters["dropped_missing_identifier"] += 1
            continue
        if genotype_id.startswith("-") or identifier.startswith("-"):
            counters["dropped_negative_identifier"] += 1
            continue
        if identifier in excluded or genotype_id in excluded:
            counters["dropped_explicit_exclusion"] += 1
            continue
        if genotype_id not in fam_ids:
            counters["dropped_not_in_fam"] += 1
            continue
        flag = flag_of.get(identifier)
        if flag is None:
            counters["dropped_no_kinship_flag"] += 1
            flag_histogram["<absent from flag file>"] += 1
            continue
        flag_histogram[flag if flag != "" else "<blank>"] += 1
        if flag not in keep_values:
            counters["dropped_related_or_unassessed"] += 1
            dropped_by_flag_value[flag if flag != "" else "<blank>"] += 1
            continue
        failed_qc = False
        for name in args.qc_flag:
            value = qc_of.get(identifier, {}).get(name, "")
            if value not in {"", "0", "0.0"}:
                failed_qc = True
                break
        if failed_qc:
            counters["dropped_qc_flag"] += 1
            continue
        if genotype_id in kept:
            counters["dropped_duplicate_genotype_id"] += 1
            continue
        kept[genotype_id] = {"eid": identifier, **row}
        counters["retained"] += 1

    if not kept:
        raise SystemExit("No participants retained; check the flag column and keep values")

    strata_columns = list(args.strata_column)
    strata_of: dict[str, str] = {}
    for genotype_id, row in kept.items():
        parts = []
        for column in strata_columns:
            value = (row.get(column) or "").strip()
            parts.append(value if value else "NA")
        if args.batch_column:
            parts.append(array_from_batch(row.get(args.batch_column, "")))
        strata_of[genotype_id] = "|".join(parts) if parts else "ALL"

    units = {genotype_id: 1 for genotype_id in kept}
    assignments = assign_splits(units, strata_of, args.seed, ratios)

    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["sample_id", "eid", "component_id", "component_size", "split", "stratum", "split_seed"]
        )
        for genotype_id in sorted(kept):
            writer.writerow([
                genotype_id,
                kept[genotype_id]["eid"],
                f"UN_{genotype_id}",
                1,
                assignments[genotype_id],
                strata_of[genotype_id],
                args.seed,
            ])

    per_stratum: dict[str, Counter[str]] = defaultdict(Counter)
    overall: Counter[str] = Counter()
    for genotype_id, split in assignments.items():
        overall[split] += 1
        per_stratum[strata_of[genotype_id]][split] += 1
    total = sum(overall.values())

    summary = {
        "classification": "TASK_B_UNRELATED_SPLIT",
        "identifiers_included": False,
        "design": "unrelated_subset_singleton_components",
        "relatedness_axis": {
            "source": "UKB field 22021 (genetic kinship to other participants)",
            "kept_values": sorted(keep_values),
            "observed_value_counts": dict(sorted(flag_histogram.items())),
            "dropped_by_value": dict(sorted(dropped_by_flag_value.items())),
            "coding_reference": FLAG_CODING_REFERENCE,
            "unknown_relatedness_is_dropped_not_kept": True,
            "pairwise_edges_available": False,
            "claim_scope": "new unrelated individual, not new family",
            "caveat": (
                "Field 22021 is a per-participant summary, not a pairwise edge list. "
                "Relatedness-disjointness holds only as strongly as UKB's own kinship "
                "inference; a family-level claim still needs the pairwise file."
            ),
        },
        "inputs": {
            "covariates": str(args.covariates),
            "covariates_sha256": sha256_file(args.covariates),
            "kinship_flag": str(args.kinship_flag),
            "kinship_flag_sha256": sha256_file(args.kinship_flag),
            "fam": str(args.fam),
            "fam_sample_count": len(fam_ids),
            "explicit_exclusions_supplied": len(excluded),
        },
        "filtering": dict(sorted(counters.items())),
        "split": {
            "seed": args.seed,
            "sample_counts": dict(sorted(overall.items())),
            "sample_fractions": {k: round(v / total, 6) for k, v in sorted(overall.items())},
            "target_fractions": ratios,
            "strata_columns": strata_columns + ([args.batch_column] if args.batch_column else []),
            "strata_count": len(per_stratum),
            "per_stratum_sample_counts": {
                stratum: dict(sorted(counts.items()))
                for stratum, counts in sorted(per_stratum.items())
            },
        },
        "outputs": {
            "tq_split_manifest": str(manifest_path),
            "tq_split_manifest_sha256": sha256_file(manifest_path),
        },
        "test_opened": False,
    }
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": "OK",
        "retained": total,
        "split_sample_counts": summary["split"]["sample_counts"],
        "strata_count": len(per_stratum),
        "manifest": str(manifest_path),
        "manifest_sha256": summary["outputs"]["tq_split_manifest_sha256"],
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
