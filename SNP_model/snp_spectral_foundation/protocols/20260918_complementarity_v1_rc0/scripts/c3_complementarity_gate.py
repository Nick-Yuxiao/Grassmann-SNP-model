#!/usr/bin/env python3
"""C3: the Concept Gate — A / B / P / F / PF under one decoder.

The gate asks two questions that Task Qualification cannot answer:

    delta F | P  =  R2(PF) - R2(P)   does functional annotation add anything
                                     once population structure is already there
    delta P | F  =  R2(PF) - R2(F)   and the other way round

Arms:
    A   covariates only
    B   A plus the train-fitted marginal-effect score at the FROZEN p-threshold
    P   A plus block-local population components   (c1)
    F   A plus functional annotation aggregates    (c2)
    PF  A plus both

Every arm uses the same ridge decoder and the same penalty grid, and every arm
picks its own penalty by K-fold cross-validation INSIDE the training rows. The
evaluation split never takes part in fitting or in penalty selection.

The B arm is recomputed at the threshold Task Qualification already froze. No
threshold is re-selected here, so this run adds no further selection on top of
the one validation already carries.

Sealed splits are not openable from this script at all: there is no flag for it.
The gate is a development instrument, and its verdict is developmental. Validation
has by now been used twice — once to pick the B threshold, once for this gate — so
a PASS here is a reason to build the confirmatory run, not a confirmatory result.

Exports are aggregate only: per-arm R2 and per-comparison deltas with CIs. No
participant row, pseudonymous or otherwise, leaves this script.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_arms import (  # noqa: E402
    LAMBDA_GRID,
    fit_ridge,
    numeric,
    open_text,
    out_of_sample_r2,
    paired_bootstrap,
    pearson,
    read_tsv,
    select_penalty_by_inner_cv,
    sha256_file,
)

# Splits this script must never read. There is deliberately no override flag.
SEALED_SPLITS = {"tq_test", "test", "bridge_holdout", "holdout"}

COMPARISONS = [
    ("PF_minus_P", "PF", "P", "delta F | P: functional annotation beyond population structure"),
    ("PF_minus_F", "PF", "F", "delta P | F: population structure beyond functional annotation"),
    ("PF_minus_A", "PF", "A", "both carriers beyond covariates"),
    ("PF_minus_B", "PF", "B", "both carriers beyond the frozen dosage score"),
    ("P_minus_A", "P", "A", "population structure beyond covariates"),
    ("F_minus_A", "F", "A", "functional annotation beyond covariates"),
    ("B_minus_A", "B", "A", "frozen dosage score beyond covariates, recomputed here"),
]


def normal_sf(z: np.ndarray) -> np.ndarray:
    return np.array([math.erfc(abs(float(value)) / math.sqrt(2.0)) for value in z])


def read_table(path: Path, id_column: str, columns: list[str]) -> dict[str, dict[str, str]]:
    table: dict[str, dict[str, str]] = {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit(f"No header in {path}")
        missing = [name for name in [id_column, *columns] if name not in reader.fieldnames]
        if missing:
            raise SystemExit(f"{path.name} is missing columns: {missing}")
        for row in reader:
            key = (row.get(id_column) or "").strip()
            if key:
                table[key] = row
    return table


def build_covariate_design(
    keys: list[str],
    covariates: dict[str, dict[str, str]],
    numeric_columns: list[str],
    categorical_columns: list[str],
    train_mask: np.ndarray,
    squared_columns: list[str],
) -> tuple[np.ndarray, list[str]]:
    """Arm A, built exactly as Task Qualification builds it."""
    blocks = [np.ones((len(keys), 1))]
    names = ["intercept"]
    for column in numeric_columns:
        values = np.array([numeric(covariates[key].get(column, "")) for key in keys])
        finite = values[train_mask][np.isfinite(values[train_mask])]
        fill = float(finite.mean()) if finite.size else 0.0
        values = np.where(np.isfinite(values), values, fill)
        centre, scale = float(values[train_mask].mean()), float(values[train_mask].std())
        scaled = (values - centre) / (scale if scale > 0 else 1.0)
        blocks.append(scaled.reshape(-1, 1))
        names.append(column)
        if column in squared_columns:
            blocks.append((scaled ** 2).reshape(-1, 1))
            names.append(f"{column}^2")
    for column in categorical_columns:
        raw = [str(covariates[key].get(column, "") or "NA").strip() for key in keys]
        levels = sorted({raw[i] for i in range(len(raw)) if train_mask[i]})
        for level in levels[1:]:  # first level is the reference
            blocks.append(
                np.array([1.0 if value == level else 0.0 for value in raw]).reshape(-1, 1)
            )
            names.append(f"{column}={level}")
    return np.hstack(blocks), names


def transform_phenotype(y: np.ndarray, train_mask: np.ndarray, method: str) -> np.ndarray:
    train_values = y[train_mask]
    if method == "none":
        return y
    if method == "zscore":
        centre, scale = float(train_values.mean()), float(train_values.std())
        clipped = np.clip(y, centre - 5 * scale, centre + 5 * scale)
        return (clipped - centre) / (scale if scale > 0 else 1.0)
    if method == "rint":
        ordered = np.sort(train_values)
        ranks = np.searchsorted(ordered, y, side="left") + 0.5
        quantiles = np.clip(ranks / (ordered.size + 1.0), 1e-6, 1 - 1e-6)
        grid = np.linspace(-6.0, 6.0, 24001)
        cdf = 0.5 * (1.0 + np.array([math.erf(value / math.sqrt(2.0)) for value in grid]))
        return np.interp(quantiles, cdf, grid)
    raise SystemExit(f"Unknown --phenotype-transform: {method}")


def frozen_dosage_score(
    panel: np.ndarray,
    columns: np.ndarray,
    design: np.ndarray,
    y: np.ndarray,
    train_mask: np.ndarray,
    threshold: float,
    block: int,
) -> tuple[np.ndarray, int]:
    """Rebuild arm B's score at one frozen threshold. Nothing is selected here.

    Same two passes as Task Qualification: train-only marginal effects on
    covariate-residualised dosage, then one score over the variants that clear the
    already-frozen p-value.
    """
    gram = design[train_mask].T @ design[train_mask]
    gram.flat[:: gram.shape[0] + 1] += 1e-8 * float(np.trace(gram)) / gram.shape[0]
    gram_inverse = np.linalg.inv(gram)
    design_train = design[train_mask]

    beta_a = np.linalg.solve(gram, design_train.T @ y[train_mask])
    residual_train = y[train_mask] - design_train @ beta_a

    def residualise(raw: np.ndarray) -> np.ndarray:
        coefficients = gram_inverse @ (design_train.T @ raw[train_mask])
        return raw - design @ coefficients

    n_variants = panel.shape[0]
    betas = np.zeros(n_variants)
    pvalues = np.ones(n_variants)
    frequencies = np.zeros(n_variants)
    degrees = max(int(train_mask.sum()) - design.shape[1] - 1, 1)

    for start in range(0, n_variants, block):
        stop = min(start + block, n_variants)
        raw = np.array(panel[start:stop, :][:, columns], dtype=np.float64).T
        for index in range(raw.shape[1]):
            observed = raw[train_mask, index]
            observed = observed[observed >= 0]
            mean = float(observed.mean()) if observed.size else 0.0
            frequencies[start + index] = mean / 2.0
            missing = raw[:, index] < 0
            if missing.any():
                raw[missing, index] = mean
        centred = residualise(raw)
        block_train = centred[train_mask]
        denominator = np.einsum("ij,ij->j", block_train, block_train)
        denominator[denominator <= 0] = np.nan
        block_beta = (block_train.T @ residual_train) / denominator
        sse = float(residual_train @ residual_train) - (block_beta ** 2) * denominator
        standard_error = np.sqrt(np.maximum(sse / degrees, 0) / denominator)
        with np.errstate(invalid="ignore", divide="ignore"):
            z = np.where(standard_error > 0, block_beta / standard_error, 0.0)
        betas[start:stop] = np.nan_to_num(block_beta)
        pvalues[start:stop] = np.nan_to_num(normal_sf(z), nan=1.0)

    score = np.zeros(len(columns))
    selected = 0
    for start in range(0, n_variants, block):
        stop = min(start + block, n_variants)
        keep = pvalues[start:stop] < threshold
        if not keep.any():
            continue
        raw = np.array(panel[start:stop, :][:, columns], dtype=np.float64).T
        for index in range(raw.shape[1]):
            missing = raw[:, index] < 0
            if missing.any():
                raw[missing, index] = 2.0 * frequencies[start + index]
        centred = residualise(raw)
        score += centred[:, keep] @ betas[start:stop][keep]
        selected += int(keep.sum())
    if selected == 0:
        raise SystemExit(
            f"The frozen threshold {threshold:g} selects no variant. "
            "Check that the threshold and the panel come from the same qualification run."
        )
    return score, selected


def standardise_on_train(block: np.ndarray, train_mask: np.ndarray) -> np.ndarray:
    centre = block[train_mask].mean(axis=0)
    scale = block[train_mask].std(axis=0)
    scale[scale == 0] = 1.0
    return (block - centre) / scale


def load_feature_block(
    features_path: Path,
    samples_path: Path | None,
    panel_samples: list[dict[str, str]],
    usable: np.ndarray,
    label: str,
) -> np.ndarray:
    block = np.load(features_path)
    if block.shape[0] != len(panel_samples):
        raise SystemExit(
            f"{label} has {block.shape[0]} rows but the panel has {len(panel_samples)} "
            "participants; they must be in panel order."
        )
    if samples_path is not None and samples_path.exists():
        listed = [row["sample_id"] for row in read_tsv(samples_path)]
        expected = [row["sample_id"] for row in panel_samples]
        if listed != expected:
            raise SystemExit(f"{label} row order does not match panel_samples.tsv")
    return np.asarray(block[usable], dtype=np.float64)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--panel-dir", type=Path, required=True)
    parser.add_argument("--covariates", type=Path, required=True)
    parser.add_argument("--covariate-id-column", default="eid")
    parser.add_argument("--covariate-column", action="append", default=[])
    parser.add_argument("--categorical-column", action="append", default=[])
    parser.add_argument("--square-column", action="append", default=[])
    parser.add_argument("--phenotype", type=Path, required=True)
    parser.add_argument("--phenotype-id-column", default="eid")
    parser.add_argument("--trait-column", required=True)
    parser.add_argument("--phenotype-transform", default="zscore",
                        choices=["none", "zscore", "rint"])
    parser.add_argument("--population-features", type=Path, required=True,
                        help="P_population.npy from c1.")
    parser.add_argument("--functional-features", type=Path, required=True,
                        help="F_functional.npy from c2.")
    parser.add_argument("--tq-results", type=Path, default=None,
                        help="TQ_RESULTS.json for this trait; supplies the frozen threshold.")
    parser.add_argument("--frozen-threshold", type=float, default=None,
                        help="Frozen p-threshold, if not taken from --tq-results.")
    parser.add_argument("--split-manifest", type=Path, default=None,
                        help="Recorded for provenance only.")
    parser.add_argument("--train-split", default="train")
    parser.add_argument("--eval-split", default="validation",
                        help="Split the gate reports on. Sealed splits are refused.")
    parser.add_argument("--inner-folds", type=int, default=5)
    parser.add_argument("--bootstrap", type=int, default=2000)
    parser.add_argument("--block", type=int, default=256)
    parser.add_argument("--seed", type=int, default=20260918)
    parser.add_argument("--max-design-gb", type=float, default=24.0,
                        help="Refuse to build a design larger than this. Lower the number of "
                             "P columns with c1 --block-size if you hit it.")
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.eval_split in SEALED_SPLITS:
        raise SystemExit(
            f"--eval-split {args.eval_split} is sealed. The Concept Gate runs on validation; "
            "opening a sealed split needs the confirmatory protocol, not this script."
        )
    if (args.tq_results is None) == (args.frozen_threshold is None):
        raise SystemExit("Give exactly one of --tq-results or --frozen-threshold")

    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    result_path = out_dir / "GATE_RESULTS.json"
    if result_path.exists() and not args.overwrite:
        raise SystemExit(f"Refusing to overwrite: {result_path}")

    threshold_source = "--frozen-threshold"
    threshold = args.frozen_threshold
    tq_sha = None
    if args.tq_results is not None:
        record = json.loads(args.tq_results.read_text(encoding="utf-8"))
        threshold = float(record["selected_p_threshold"])
        threshold_source = f"{args.tq_results.name}:selected_p_threshold"
        tq_sha = sha256_file(args.tq_results)
        if record.get("trait_column") not in (None, args.trait_column):
            raise SystemExit(
                f"{args.tq_results.name} was written for trait {record.get('trait_column')}, "
                f"not {args.trait_column}"
            )

    panel_dir = args.panel_dir.resolve()
    panel = np.load(panel_dir / "panel.int8.npy", mmap_mode="r")
    with open_text(panel_dir / "panel_samples.tsv") as handle:
        panel_samples = list(csv.DictReader(handle, delimiter="\t"))
    if panel.shape[1] != len(panel_samples):
        raise SystemExit("panel_samples.tsv does not match the panel matrix")

    covariates = read_table(args.covariates, args.covariate_id_column,
                            args.covariate_column + args.categorical_column)
    phenotypes = read_table(args.phenotype, args.phenotype_id_column, [args.trait_column])

    usable: list[int] = []
    keys: list[str] = []
    values: list[float] = []
    dropped = {"sealed_split": 0, "other_split": 0, "no_covariates": 0,
               "no_phenotype": 0, "missing_trait": 0}
    for column, row in enumerate(panel_samples):
        split = row.get("split", "")
        if split in SEALED_SPLITS:
            dropped["sealed_split"] += 1
            continue
        if split not in (args.train_split, args.eval_split):
            dropped["other_split"] += 1
            continue
        key = (row.get("eid") or "").strip()
        if key not in covariates:
            dropped["no_covariates"] += 1
            continue
        if key not in phenotypes:
            dropped["no_phenotype"] += 1
            continue
        value = numeric(phenotypes[key].get(args.trait_column, ""))
        if not math.isfinite(value):
            dropped["missing_trait"] += 1
            continue
        usable.append(column)
        keys.append(key)
        values.append(value)
    if len(usable) < 1000:
        raise SystemExit(f"Only {len(usable)} usable participants; refusing to run the gate")

    columns = np.array(usable, dtype=np.int64)
    split_labels = np.array([panel_samples[index]["split"] for index in usable])
    train_mask = split_labels == args.train_split
    eval_mask = split_labels == args.eval_split
    if not train_mask.any() or not eval_mask.any():
        raise SystemExit(f"Both {args.train_split} and {args.eval_split} participants are required")

    y = transform_phenotype(np.array(values, dtype=np.float64), train_mask,
                            args.phenotype_transform)
    train_mean = float(y[train_mask].mean())

    design_a, design_names = build_covariate_design(
        keys, covariates, args.covariate_column, args.categorical_column,
        train_mask, args.square_column,
    )

    population = standardise_on_train(
        load_feature_block(args.population_features,
                           args.population_features.parent / "P_samples.tsv",
                           panel_samples, columns, "P_population.npy"),
        train_mask,
    )
    functional = standardise_on_train(
        load_feature_block(args.functional_features, None, panel_samples, columns,
                           "F_functional.npy"),
        train_mask,
    )

    widest = design_a.shape[1] + population.shape[1] + functional.shape[1]
    estimated_gb = len(usable) * widest * 8 / (1 << 30)
    if estimated_gb > args.max_design_gb:
        raise SystemExit(
            f"The PF design would need about {estimated_gb:.1f} GB "
            f"({len(usable)} rows x {widest} columns). Rebuild P with a larger "
            "c1 --block-size, or raise --max-design-gb if the machine really has it."
        )

    score, selected_variants = frozen_dosage_score(
        panel, columns, design_a, y, train_mask, threshold, args.block
    )
    score = standardise_on_train(score.reshape(-1, 1), train_mask)

    # Built one at a time and released: holding all five at once would cost about
    # twice the PF design in memory for no reason.
    arm_blocks = {
        "A": [],
        "B": [score],
        "P": [population],
        "F": [functional],
        "PF": [population, functional],
    }

    arm_rows: list[dict[str, object]] = []
    predictions: dict[str, np.ndarray] = {}
    penalty_traces: dict[str, list[dict[str, float]]] = {}
    for name, extra in arm_blocks.items():
        design = np.hstack([design_a, *extra]) if extra else design_a
        penalty, trace = select_penalty_by_inner_cv(
            design[train_mask], y[train_mask], args.inner_folds, args.seed, LAMBDA_GRID
        )
        train_block = design[train_mask]
        beta = fit_ridge(train_block.T @ train_block, train_block.T @ y[train_mask], penalty)
        prediction = design @ beta
        predictions[name] = prediction
        penalty_traces[name] = trace
        arm_rows.append({
            "trait": args.trait_column,
            "arm": name,
            "r2": round(out_of_sample_r2(y[eval_mask], prediction[eval_mask], train_mean), 6),
            "pearson": round(pearson(y[eval_mask], prediction[eval_mask]), 6),
            "n_features": int(design.shape[1] - 1),
            "penalty": penalty,
            "penalty_at_grid_edge": int(penalty in (LAMBDA_GRID[0], LAMBDA_GRID[-1])),
            "n_eval": int(eval_mask.sum()),
            "eval_split": args.eval_split,
        })
        if design is not design_a:
            del design

    comparison_rows: list[dict[str, object]] = []
    for index, (label, left, right, reading) in enumerate(COMPARISONS):
        result = paired_bootstrap(
            y[eval_mask], predictions[left][eval_mask], predictions[right][eval_mask],
            train_mean, args.bootstrap, args.seed + index,
        )
        comparison_rows.append({
            "trait": args.trait_column,
            "comparison": label,
            "delta_r2": round(result["delta_r2"], 6),
            "ci_low": round(result["ci_low"], 6),
            "ci_high": round(result["ci_high"], 6),
            "replicates": result["replicates"],
            "seed": result["bootstrap_seed"],
            "positive_and_ci_excludes_zero": int(
                result["delta_r2"] > 0 and result["ci_low"] > 0
            ),
            "reading": reading,
        })
    by_label = {row["comparison"]: row for row in comparison_rows}

    def passed(label: str) -> bool:
        return bool(by_label[label]["positive_and_ci_excludes_zero"])

    functional_adds = passed("PF_minus_P")
    population_adds = passed("PF_minus_F")
    verdict = (
        "BOTH_CARRIERS_ADD" if functional_adds and population_adds
        else "ONLY_FUNCTIONAL_ADDS" if functional_adds
        else "ONLY_POPULATION_ADDS" if population_adds
        else "NEITHER_CARRIER_ADDS"
    )

    report = {
        "classification": "COMPLEMENTARITY_CONCEPT_GATE",
        "status": "DEVELOPMENTAL_NOT_CONFIRMATORY",
        "identifiers_included": False,
        "individual_level_rows_exported": False,
        "question": "does functional annotation add beyond population structure, and conversely",
        "not_answered": [
            "whether self-supervised pretraining reduces phenotype sample complexity",
            "whether any learned representation beats raw dosage",
            "anything confirmatory: validation now carries two selections",
        ],
        "why_developmental": (
            "validation chose the B threshold during Task Qualification and is the evaluation "
            "split here as well. A PASS is grounds to pre-register a confirmatory run on a "
            "sealed split, not a confirmed effect."
        ),
        "trait_column": args.trait_column,
        "phenotype_transform": args.phenotype_transform,
        "splits": {
            "train": args.train_split,
            "eval": args.eval_split,
            "sealed_and_never_read": sorted(SEALED_SPLITS),
        },
        "decoder": {
            "model": "ridge, intercept unpenalised",
            "penalty_grid": LAMBDA_GRID,
            "penalty_selected_by": f"{args.inner_folds}-fold CV inside train",
            "penalty_never_selected_on": args.eval_split,
            "note": (
                "arm A is ridge here and was ordinary least squares in Task Qualification, so "
                "its R2 need not match the qualification number exactly"
            ),
        },
        "frozen_threshold": {
            "value": threshold,
            "source": threshold_source,
            "tq_results_sha256": tq_sha,
            "variants_selected": selected_variants,
            "reselected_here": False,
        },
        "participants": {
            "usable": len(usable),
            "train": int(train_mask.sum()),
            "eval": int(eval_mask.sum()),
            "dropped": dropped,
        },
        "design_terms_arm_A": design_names,
        "feature_blocks": {
            "P": {"columns": int(population.shape[1]),
                  "sha256": sha256_file(args.population_features)},
            "F": {"columns": int(functional.shape[1]),
                  "sha256": sha256_file(args.functional_features)},
        },
        "inputs": {
            "panel_sha256": sha256_file(panel_dir / "panel.int8.npy"),
            "split_manifest_sha256": (
                sha256_file(args.split_manifest) if args.split_manifest else None
            ),
        },
        "arms": arm_rows,
        "comparisons": comparison_rows,
        "verdict": verdict,
        "verdict_rule": "delta_r2 > 0 and paired 95% CI lower bound > 0",
        "seed": args.seed,
    }
    result_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n",
                           encoding="utf-8")
    (out_dir / "PENALTY_TRACE.json").write_text(
        json.dumps(penalty_traces, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )

    arm_fields = ["trait", "arm", "r2", "pearson", "n_features", "penalty",
                  "penalty_at_grid_edge", "n_eval", "eval_split"]
    with (out_dir / "arms.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=arm_fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(arm_rows)
    comparison_fields = ["trait", "comparison", "delta_r2", "ci_low", "ci_high",
                         "replicates", "seed", "positive_and_ci_excludes_zero", "reading"]
    with (out_dir / "comparisons.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=comparison_fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(comparison_rows)

    print(json.dumps({
        "status": "OK",
        "trait": args.trait_column,
        "verdict": verdict,
        "developmental": True,
        "arms": [{k: row[k] for k in ("arm", "r2", "n_features", "penalty")} for row in arm_rows],
        "key_comparisons": [
            {k: by_label[label][k] for k in ("comparison", "delta_r2", "ci_low", "ci_high")}
            for label in ("PF_minus_P", "PF_minus_F")
        ],
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
