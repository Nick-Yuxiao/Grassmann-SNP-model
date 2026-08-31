from __future__ import annotations

import argparse
import collections
import gzip
import json
import math
from pathlib import Path


def count_rows(path: Path) -> int:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return sum(1 for line in handle if line.strip())


def quantiles(values: list[float], probs: list[float], names: list[str]) -> dict[str, float]:
    if not values:
        return {}
    ordered = sorted(values)
    result: list[float] = []
    for probability in probs:
        index = (len(ordered) - 1) * probability
        lower = math.floor(index)
        upper = math.ceil(index)
        if lower == upper:
            result.append(float(ordered[lower]))
        else:
            weight = index - lower
            result.append(float(ordered[lower] * (1 - weight) + ordered[upper] * weight))
    return dict(zip(names, result))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--all-metrics", type=Path, required=True)
    parser.add_argument("--maf-sites", type=Path, required=True)
    parser.add_argument("--source-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    positions: collections.Counter[tuple[str, int]] = collections.Counter()
    missing: list[float] = []
    mafs: list[float] = []

    with gzip.open(args.maf_sites, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise SystemExit(f"line {line_number}: expected 9 fields, observed {len(fields)}")
            chrom, pos, _ref, _alt, _ac, _an, _af, maf, f_missing = fields
            positions[(chrom, int(pos))] += 1
            mafs.append(float(maf))
            missing.append(float(f_missing))

    duplicate_sizes = [value for value in positions.values() if value > 1]
    all_count = count_rows(args.all_metrics)
    status = "PASS" if all_count > 0 and mafs else "FAIL_EMPTY"
    payload = {
        "schema_version": "1.0",
        "protocol_version": "v7.1.2",
        "stage": "SOURCE_PREFILTERED_FILTER_UNSET_biallelic_SNP_donor_MAF_gt_0.01",
        "source_bcf_sha256": args.source_sha256,
        "filter_semantics": "VCF_FILTER_DOT_UNSET_NOT_LITERAL_PASS",
        "pass_label_used": False,
        "source_prefiltered_claim_is_provenance_only": True,
        "unset_filter_biallelic_snp_records": all_count,
        "maf_gt_0.01_records": len(mafs),
        "unique_positions": len(positions),
        "duplicate_positions": len(duplicate_sizes),
        "records_at_duplicate_positions": sum(duplicate_sizes),
        "donor_missing_fraction": {
            "equal_0": sum(value == 0 for value in missing),
            "le_0.001": sum(value <= 0.001 for value in missing),
            "le_0.01": sum(value <= 0.01 for value in missing),
            "gt_0.01": sum(value > 0.01 for value in missing),
            "quantiles": quantiles(
                missing,
                [0, 0.5, 0.9, 0.95, 0.99, 1],
                ["q0", "q50", "q90", "q95", "q99", "q100"],
            ),
        },
        "maf_quantiles": quantiles(
            mafs,
            [0, 0.1, 0.5, 0.9, 1],
            ["q0", "q10", "q50", "q90", "q100"],
        ),
        "status": status,
    }
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))
    if status != "PASS":
        raise SystemExit("empty stage-1 result")


if __name__ == "__main__":
    main()
