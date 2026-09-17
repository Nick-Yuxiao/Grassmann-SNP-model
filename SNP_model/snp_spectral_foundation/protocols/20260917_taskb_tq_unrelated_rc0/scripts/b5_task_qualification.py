#!/usr/bin/env python3
"""B5: Task Qualification for the phenotype bridge — arm A vs arm B.

Arm A is covariates only. Arm B is covariates plus raw additive dosage, entered
through a train-fitted linear predictor over the frozen panel. The question is
narrow and pre-registered: on held-out unrelated participants, does raw genotype
add out-of-sample phenotype information over covariates alone?

Qualification is not the bridge. A PASS says the trait carries detectable
genotype signal on this panel and split, which is the precondition for asking
whether a learned representation beats raw dosage. It says nothing about any
representation, about Grassmann, or about causality.

Firewall: validation selects the single p-value threshold; the test split stays
closed until --allow-test is passed, and the run writes a marker so a second
opening is visible.
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
from lib_split import open_text, sha256_file  # noqa: E402

P_THRESHOLDS = [1.0, 0.5, 0.1, 0.05, 0.01, 1e-3, 1e-4, 1e-5, 1e-6]
DEFAULT_BLOCK = 256


def normal_sf(z: np.ndarray) -> np.ndarray:
    return np.array([math.erfc(abs(float(value)) / math.sqrt(2.0)) for value in z])


def fit_linear(design: np.ndarray, target: np.ndarray) -> np.ndarray:
    gram = design.T @ design
    ridge = 1e-8 * float(np.trace(gram)) / max(gram.shape[0], 1)
    gram.flat[:: gram.shape[0] + 1] += ridge
    return np.linalg.solve(gram, design.T @ target)


def out_of_sample_r2(y: np.ndarray, prediction: np.ndarray, train_mean: float) -> float:
    residual = float(np.sum((y - prediction) ** 2))
    total = float(np.sum((y - train_mean) ** 2))
    return 1.0 - residual / total if total > 0 else float("nan")


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


def numeric(value: str) -> float:
    text = (value or "").strip()
    if not text:
        return float("nan")
    try:
        return float(text)
    except ValueError:
        return float("nan")


def build_design(
    keys: list[str],
    covariates: dict[str, dict[str, str]],
    numeric_columns: list[str],
    categorical_columns: list[str],
    train_mask: np.ndarray,
    squared_columns: list[str],
) -> tuple[np.ndarray, list[str]]:
    blocks = [np.ones((len(keys), 1))]
    names = ["intercept"]
    for column in numeric_columns:
        values = np.array([numeric(covariates[key].get(column, "")) for key in keys])
        train_values = values[train_mask]
        finite = train_values[np.isfinite(train_values)]
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
        levels = sorted({raw[index] for index in range(len(raw)) if train_mask[index]})
        for level in levels[1:]:  # first level is the reference
            blocks.append(np.array([1.0 if value == level else 0.0 for value in raw]).reshape(-1, 1))
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
        # Acklam-free inverse normal via erfinv identity on a fine grid.
        grid = np.linspace(-6.0, 6.0, 24001)
        cdf = 0.5 * (1.0 + np.array([math.erf(value / math.sqrt(2.0)) for value in grid]))
        return np.interp(quantiles, cdf, grid)
    raise SystemExit(f"Unknown --phenotype-transform: {method}")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
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
    parser.add_argument("--sesoi-delta-r2", type=float, default=0.005,
                        help="Pre-registered smallest incremental out-of-sample R2 of interest.")
    parser.add_argument("--bootstrap", type=int, default=2000)
    parser.add_argument("--block", type=int, default=DEFAULT_BLOCK,
                        help="Variants held in memory per pass; lower it if RAM is tight.")
    parser.add_argument("--seed", type=int, default=20260917)
    parser.add_argument("--allow-test", action="store_true",
                        help="Open the test split. Without it only validation is reported.")
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    result_path = out_dir / "TQ_RESULTS.json"
    marker_path = out_dir / "TEST_OPENED.marker"
    if result_path.exists() and not args.overwrite:
        raise SystemExit(f"Refusing to overwrite: {result_path}")
    if args.allow_test and marker_path.exists() and not args.overwrite:
        raise SystemExit(
            f"Test split already opened for this run directory: {marker_path}. "
            "A second opening needs a new output directory and a recorded reason."
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
    dropped = {"no_covariates": 0, "no_phenotype": 0, "missing_trait": 0}
    for column, row in enumerate(panel_samples):
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
        raise SystemExit(f"Only {len(usable)} usable participants; refusing to qualify on this")

    columns = np.array(usable, dtype=np.int64)
    split_labels = np.array([panel_samples[index]["split"] for index in usable])
    train_mask = split_labels == "train"
    validation_mask = split_labels == "validation"
    test_mask = split_labels == "test"
    if not train_mask.any() or not validation_mask.any():
        raise SystemExit("Train and validation participants are both required")

    y = transform_phenotype(np.array(values, dtype=np.float64), train_mask,
                            args.phenotype_transform)
    design, design_names = build_design(
        keys, covariates, args.covariate_column, args.categorical_column,
        train_mask, args.square_column,
    )

    beta_a = fit_linear(design[train_mask], y[train_mask])
    prediction_a = design @ beta_a
    train_mean = float(y[train_mask].mean())
    residual = y - prediction_a

    # Pass 1: train-only marginal effects on covariate-residualised dosage.
    gram = design[train_mask].T @ design[train_mask]
    gram.flat[:: gram.shape[0] + 1] += 1e-8 * float(np.trace(gram)) / gram.shape[0]
    gram_inverse = np.linalg.inv(gram)
    design_train = design[train_mask]
    residual_train = residual[train_mask]
    n_variants = panel.shape[0]
    betas = np.zeros(n_variants)
    pvalues = np.ones(n_variants)
    train_frequency = np.zeros(n_variants)

    def residualise(block: np.ndarray) -> np.ndarray:
        coefficients = gram_inverse @ (design_train.T @ block[train_mask])
        return block - design @ coefficients

    degrees = int(train_mask.sum()) - design.shape[1] - 1
    for start in range(0, n_variants, args.block):
        stop = min(start + args.block, n_variants)
        raw = np.array(panel[start:stop, :][:, columns], dtype=np.float64).T
        missing = raw < 0
        frequency = np.zeros(raw.shape[1])
        for index in range(raw.shape[1]):
            observed = raw[train_mask, index]
            observed = observed[observed >= 0]
            frequency[index] = observed.mean() if observed.size else 0.0
            if missing[:, index].any():
                raw[missing[:, index], index] = frequency[index]
        train_frequency[start:stop] = frequency / 2.0
        centred = residualise(raw)
        block_train = centred[train_mask]
        denominator = np.einsum("ij,ij->j", block_train, block_train)
        denominator[denominator <= 0] = np.nan
        block_beta = (block_train.T @ residual_train) / denominator
        # sse = sum(r^2) - beta^2 * sum(g^2); avoids materialising an n x block residual
        sse = float(residual_train @ residual_train) - (block_beta ** 2) * denominator
        standard_error = np.sqrt(np.maximum(sse / max(degrees, 1), 0) / denominator)
        with np.errstate(invalid="ignore", divide="ignore"):
            z = np.where(standard_error > 0, block_beta / standard_error, 0.0)
        betas[start:stop] = np.nan_to_num(block_beta)
        pvalues[start:stop] = np.nan_to_num(normal_sf(z), nan=1.0)

    # Pass 2: one score per threshold, weights and residualisation both train-fitted.
    scores = np.zeros((len(P_THRESHOLDS), len(columns)))
    selected_counts = np.zeros(len(P_THRESHOLDS), dtype=int)
    for start in range(0, n_variants, args.block):
        stop = min(start + args.block, n_variants)
        raw = np.array(panel[start:stop, :][:, columns], dtype=np.float64).T
        for index in range(raw.shape[1]):
            column_missing = raw[:, index] < 0
            if column_missing.any():
                raw[column_missing, index] = 2.0 * train_frequency[start + index]
        centred = residualise(raw)
        block_beta = betas[start:stop]
        block_p = pvalues[start:stop]
        for level, threshold in enumerate(P_THRESHOLDS):
            keep = block_p < threshold
            if not keep.any():
                continue
            scores[level] += centred[:, keep] @ block_beta[keep]
            selected_counts[level] += int(keep.sum())

    def evaluate(mask: np.ndarray, score: np.ndarray) -> tuple[float, float, np.ndarray]:
        augmented = np.hstack([design, score.reshape(-1, 1)])
        beta_b = fit_linear(augmented[train_mask], y[train_mask])
        prediction_b = augmented @ beta_b
        return (
            out_of_sample_r2(y[mask], prediction_a[mask], train_mean),
            out_of_sample_r2(y[mask], prediction_b[mask], train_mean),
            prediction_b,
        )

    validation_rows = []
    best_level, best_gain = None, -math.inf
    for level, threshold in enumerate(P_THRESHOLDS):
        if selected_counts[level] == 0:
            continue
        r2_a, r2_b, _ = evaluate(validation_mask, scores[level])
        gain = r2_b - r2_a
        validation_rows.append({
            "p_threshold": threshold,
            "variants_selected": int(selected_counts[level]),
            "r2_A": round(r2_a, 6),
            "r2_B": round(r2_b, 6),
            "delta_r2": round(gain, 6),
        })
        if gain > best_gain:
            best_level, best_gain = level, gain
    if best_level is None:
        raise SystemExit("No p-value threshold selected any variant")

    report: dict[str, object] = {
        "classification": "TASK_B_TASK_QUALIFICATION",
        "identifiers_included": False,
        "question": "does raw additive dosage add out-of-sample phenotype information over covariates",
        "not_answered": [
            "whether any learned representation beats raw dosage",
            "biological specificity or causality",
            "anything about Grassmann geometry",
        ],
        "trait_column": args.trait_column,
        "phenotype_transform": args.phenotype_transform,
        "panel": {
            "variants": int(n_variants),
            "panel_sha256": sha256_file(panel_dir / "panel.int8.npy"),
        },
        "participants": {
            "usable": len(usable),
            "dropped": dropped,
            "train": int(train_mask.sum()),
            "validation": int(validation_mask.sum()),
            "test": int(test_mask.sum()),
        },
        "design_terms": design_names,
        "validation_scan": validation_rows,
        "selected_p_threshold": P_THRESHOLDS[best_level],
        "selected_variant_count": int(selected_counts[best_level]),
        "sesoi_delta_r2": args.sesoi_delta_r2,
        "test_opened": bool(args.allow_test),
    }

    if args.allow_test:
        if not test_mask.any():
            raise SystemExit("--allow-test given but the split has no test participants")
        r2_a, r2_b, prediction_b = evaluate(test_mask, scores[best_level])
        delta = r2_b - r2_a
        rng = np.random.default_rng(args.seed)
        indices = np.flatnonzero(test_mask)
        deltas = np.empty(args.bootstrap)
        for replicate in range(args.bootstrap):
            draw = rng.choice(indices, size=indices.size, replace=True)
            deltas[replicate] = (
                out_of_sample_r2(y[draw], prediction_b[draw], train_mean)
                - out_of_sample_r2(y[draw], prediction_a[draw], train_mean)
            )
        low, high = np.percentile(deltas, [2.5, 97.5])

        permuted = scores[best_level].copy()
        permuted[indices] = rng.permutation(permuted[indices])
        _, r2_negative, _ = evaluate(test_mask, permuted)
        negative_delta = r2_negative - r2_a

        passed = bool(delta >= args.sesoi_delta_r2 and low > 0)
        report["test"] = {
            "r2_A_covariates_only": round(r2_a, 6),
            "r2_B_covariates_plus_dosage": round(r2_b, 6),
            "delta_r2": round(delta, 6),
            "delta_r2_ci95": [round(float(low), 6), round(float(high), 6)],
            "bootstrap_replicates": args.bootstrap,
            "negative_control_permuted_score_delta_r2": round(negative_delta, 6),
            "negative_control_reading": (
                "the train-fitted score is scrambled against test outcomes; a real signal "
                "drives this clearly negative. A value at or above delta_r2 means the score "
                "and the outcome were misaligned and the run is invalid."
            ),
        }
        report["verdict"] = (
            "QUALIFIED" if passed
            else "STATISTICALLY_POSITIVE_PRACTICALLY_TIED" if low > 0
            else "NOT_QUALIFIED"
        )
        marker_path.write_text(
            json.dumps({"trait": args.trait_column, "panel": str(panel_dir)}, ensure_ascii=False)
            + "\n",
            encoding="utf-8",
        )
    else:
        report["verdict"] = "VALIDATION_ONLY_TEST_NOT_OPENED"

    result_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    printable = {key: report[key] for key in ("trait_column", "verdict", "selected_p_threshold",
                                              "selected_variant_count")}
    printable["participants"] = report["participants"]
    if "test" in report:
        printable["test"] = report["test"]
    else:
        printable["best_validation"] = validation_rows[
            max(range(len(validation_rows)), key=lambda i: validation_rows[i]["delta_r2"])
        ]
    print(json.dumps(printable, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
