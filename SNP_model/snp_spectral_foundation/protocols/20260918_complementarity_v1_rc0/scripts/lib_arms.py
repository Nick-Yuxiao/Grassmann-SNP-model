#!/usr/bin/env python3
"""Shared pieces for the complementarity arms: design, ridge, evaluation.

Every arm uses the same decoder, the same regularisation grid and the same
participants, so a difference between arms is a difference in the information the
features carry rather than in how they were fitted.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import math
from pathlib import Path
from typing import Iterator

import numpy as np

LAMBDA_GRID = [1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0, 1e3, 1e4, 1e5]


def open_text(path: Path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8-sig", errors="replace", newline="")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> Iterator[dict[str, str]]:
    with open_text(path) as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def numeric(value: str | None) -> float:
    text = (value or "").strip()
    if not text:
        return float("nan")
    try:
        return float(text)
    except ValueError:
        return float("nan")


def standardise_block(
    dosage: np.ndarray, a1_frequency: np.ndarray, missing_value: int = -1
) -> np.ndarray:
    """Centre and scale dosage with train-only allele frequencies.

    Missing entries become zero, that is the train mean after centring. The
    frequencies come from the panel's own train-only column, so no validation or
    test participant contributes to the scaling.
    """
    scaled = dosage.astype(np.float64, copy=True)
    missing = dosage == missing_value
    mean = 2.0 * a1_frequency
    scale = np.sqrt(np.maximum(2.0 * a1_frequency * (1.0 - a1_frequency), 1e-12))
    scaled -= mean
    scaled /= scale
    scaled[missing] = 0.0
    return scaled


def fit_ridge(gram: np.ndarray, moment: np.ndarray, penalty: float) -> np.ndarray:
    """Solve (X'X + penalty*I) b = X'y, leaving the intercept column unpenalised."""
    regularised = gram.copy()
    diagonal = np.arange(1, gram.shape[0])
    regularised[diagonal, diagonal] += penalty
    return np.linalg.solve(regularised, moment)


def out_of_sample_r2(y: np.ndarray, prediction: np.ndarray, train_mean: float) -> float:
    residual = float(np.sum((y - prediction) ** 2))
    total = float(np.sum((y - train_mean) ** 2))
    return 1.0 - residual / total if total > 0 else float("nan")


def pearson(y: np.ndarray, prediction: np.ndarray) -> float:
    if y.size < 2:
        return float("nan")
    centred_y, centred_p = y - y.mean(), prediction - prediction.mean()
    denominator = math.sqrt(float(centred_y @ centred_y) * float(centred_p @ centred_p))
    return float(centred_y @ centred_p) / denominator if denominator > 0 else float("nan")


def select_penalty_by_inner_cv(
    design: np.ndarray,
    y: np.ndarray,
    folds: int,
    seed: int,
    grid: list[float] | None = None,
) -> tuple[float, list[dict[str, float]]]:
    """Choose the ridge penalty by K-fold cross-validation inside the training rows.

    The penalty is never chosen on the evaluation split. That keeps the split used
    for the gate honest for every arm, since each arm picks its own penalty the
    same way and none of them sees the evaluation rows while doing it.
    """
    grid = grid or LAMBDA_GRID
    rng = np.random.default_rng(seed)
    assignment = rng.permutation(np.arange(y.size) % folds)

    total_gram = design.T @ design
    total_moment = design.T @ y
    fold_grams, fold_moments, fold_rows = [], [], []
    for fold in range(folds):
        mask = assignment == fold
        block = design[mask]
        fold_grams.append(block.T @ block)
        fold_moments.append(block.T @ y[mask])
        fold_rows.append(mask)

    scored: list[dict[str, float]] = []
    for penalty in grid:
        errors = 0.0
        for fold in range(folds):
            mask = fold_rows[fold]
            beta = fit_ridge(total_gram - fold_grams[fold], total_moment - fold_moments[fold],
                             penalty)
            residual = y[mask] - design[mask] @ beta
            errors += float(residual @ residual)
        scored.append({"penalty": penalty, "cv_sse": errors})
    best = min(scored, key=lambda item: item["cv_sse"])
    return float(best["penalty"]), scored


def paired_bootstrap(
    y: np.ndarray,
    prediction_left: np.ndarray,
    prediction_right: np.ndarray,
    train_mean: float,
    replicates: int,
    seed: int,
) -> dict[str, float]:
    """Percentile CI for R2(left) - R2(right), resampling the same rows for both."""
    rng = np.random.default_rng(seed)
    deltas = np.empty(replicates)
    size = y.size
    for replicate in range(replicates):
        draw = rng.integers(0, size, size)
        deltas[replicate] = (
            out_of_sample_r2(y[draw], prediction_left[draw], train_mean)
            - out_of_sample_r2(y[draw], prediction_right[draw], train_mean)
        )
    low, high = np.percentile(deltas, [2.5, 97.5])
    return {
        "delta_r2": float(
            out_of_sample_r2(y, prediction_left, train_mean)
            - out_of_sample_r2(y, prediction_right, train_mean)
        ),
        "ci_low": float(low),
        "ci_high": float(high),
        "replicates": replicates,
        "bootstrap_seed": seed,
    }


def canonical_correlations(left: np.ndarray, right: np.ndarray, top: int = 5) -> list[float]:
    """Largest canonical correlations between two feature blocks.

    Used to report how much a new block overlaps one already in the model, for
    instance local block PCs against the global ancestry PCs that arm A carries.
    """
    def whiten(block: np.ndarray) -> np.ndarray:
        centred = block - block.mean(axis=0)
        u, s, _ = np.linalg.svd(centred, full_matrices=False)
        keep = s > (s.max() * 1e-10 if s.size and s.max() > 0 else 0)
        return u[:, keep]

    left_basis, right_basis = whiten(left), whiten(right)
    if left_basis.size == 0 or right_basis.size == 0:
        return []
    values = np.linalg.svd(left_basis.T @ right_basis, compute_uv=False)
    return [round(float(value), 6) for value in values[:top]]
