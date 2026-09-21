#!/usr/bin/env python3
"""C1b: the full overlap spectrum between P and the global ancestry PCs.

c1 records only the top five canonical correlations, which is not enough to say
how much of P arm A already carries. With 654 block-local components against 40
global PCs there are 40 canonical pairs, and what matters is how fast they decay:
a couple near 1 means P repeats the leading ancestry axes and little else, while a
long plateau near 1 would mean P is largely redundant with arm A.

This reads the P matrix c1 already wrote, so it does not touch the panel and runs
in seconds. Train rows only, so it stays consistent with how P was fitted.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_arms import canonical_correlations, numeric, read_tsv, sha256_file  # noqa: E402


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--population-dir", type=Path, required=True,
                        help="c1 --out-dir, holding P_population.npy and P_samples.tsv.")
    parser.add_argument("--covariates", type=Path, required=True)
    parser.add_argument("--covariate-id-column", default="eid")
    parser.add_argument("--global-pc-prefix", default="n_22009_0_")
    parser.add_argument("--global-pc-count", type=int, default=40)
    parser.add_argument("--train-split", default="train")
    parser.add_argument("--rows", type=int, default=50000,
                        help="Train rows sampled for the decomposition.")
    parser.add_argument("--seed", type=int, default=20260921)
    parser.add_argument("--out", type=Path, default=None,
                        help="Where to write the report; defaults into --population-dir.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    population_dir = args.population_dir.resolve()
    features = np.load(population_dir / "P_population.npy", mmap_mode="r")

    samples = list(read_tsv(population_dir / "P_samples.tsv"))
    if len(samples) != features.shape[0]:
        raise SystemExit("P_samples.tsv does not match P_population.npy")

    wanted = [f"{args.global_pc_prefix}{i}" for i in range(1, args.global_pc_count + 1)]
    table: dict[str, dict[str, str]] = {}
    for row in read_tsv(args.covariates):
        key = (row.get(args.covariate_id_column) or "").strip()
        if key:
            table[key] = row
    present = [name for name in wanted if name in next(iter(table.values()), {})]
    if not present:
        raise SystemExit(f"No columns matching {args.global_pc_prefix}* in {args.covariates}")

    eligible = []
    for index, row in enumerate(samples):
        if row.get("split") != args.train_split:
            continue
        record = table.get(row.get("eid", ""))
        if record is None:
            continue
        values = [numeric(record.get(name)) for name in present]
        if all(np.isfinite(values)):
            eligible.append((index, values))
    if len(eligible) <= len(present):
        raise SystemExit(f"Only {len(eligible)} usable train rows; cannot decompose")

    rng = np.random.default_rng(args.seed)
    if len(eligible) > args.rows:
        picked = rng.choice(len(eligible), size=args.rows, replace=False)
        picked.sort()
        eligible = [eligible[i] for i in picked]

    rows = np.array([index for index, _ in eligible])
    global_pcs = np.array([values for _, values in eligible])
    block = np.asarray(features[rows], dtype=np.float64)

    spectrum = canonical_correlations(block, global_pcs, top=len(present))
    above = {
        f"above_{cut}": int(sum(1 for value in spectrum if value >= cut))
        for cut in (0.99, 0.95, 0.9, 0.8, 0.5)
    }

    report = {
        "classification": "COMPLEMENTARITY_P_OVERLAP_SPECTRUM",
        "identifiers_included": False,
        "phenotype_read": False,
        "question": "how much of P does arm A already carry through PC1-40",
        "p_columns": int(features.shape[1]),
        "global_pc_columns": len(present),
        "canonical_pairs": len(spectrum),
        "train_rows_used": int(rows.size),
        "spectrum": spectrum,
        "pairs_above_threshold": above,
        "redundant_fraction_of_P": round(above["above_0.95"] / features.shape[1], 6),
        "reading": (
            "Each value is one canonical pair. Pairs near 1 are directions of P that "
            "arm A already carries; the count of those against P's column total is the "
            "share of P that is redundant. A fast decay means the overlap is confined to "
            "a few ancestry axes and the rest of P is local structure arm A cannot "
            "represent. This is a diagnostic, not the verdict: whether P adds is decided "
            "by P_minus_A in c3."
        ),
        "inputs": {
            "features_sha256": sha256_file(population_dir / "P_population.npy"),
            "covariates": str(args.covariates),
        },
        "seed": args.seed,
    }

    out_path = args.out or (population_dir / "P_OVERLAP_SPECTRUM.json")
    out_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n",
                        encoding="utf-8")
    with out_path.with_suffix(".csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["pair", "canonical_correlation"])
        for rank, value in enumerate(spectrum, start=1):
            writer.writerow([rank, value])

    print(json.dumps({
        "status": "OK",
        "p_columns": report["p_columns"],
        "canonical_pairs": report["canonical_pairs"],
        "pairs_above_threshold": above,
        "redundant_fraction_of_P": report["redundant_fraction_of_P"],
        "first_ten": spectrum[:10],
        "last_five": spectrum[-5:],
        "report": str(out_path),
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
