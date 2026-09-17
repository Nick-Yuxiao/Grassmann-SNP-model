"""Asset-free core of the TQ-G1 gene-layer readout Pilot.

Everything in this module is importable and testable without the chr18
GEUVADIS assets. ``run_tqg1.py`` binds the assets and calls into here, so the
statistical contract can be unit tested on synthetic data before a real run.
"""

from __future__ import annotations

import hashlib
import math

import numpy as np


# ---------------------------------------------------------------------------
# deterministic selection
# ---------------------------------------------------------------------------


def deterministic_gene_order(gene_ids: list[str], salt: str) -> list[str]:
    """Order genes by a fixed hash so the cap never selects on phenotype signal."""
    keyed = [(hashlib.sha256(f"{salt}:{gene}".encode()).hexdigest(), gene) for gene in gene_ids]
    keyed.sort()
    return [gene for _, gene in keyed]


def uniform_thin(indices: np.ndarray, count: int) -> np.ndarray:
    """Keep ``count`` position-ordered entries spread uniformly over ``indices``."""
    if count < 1:
        raise ValueError("count must be positive")
    if len(indices) < count:
        raise ValueError("cannot thin below the requested count")
    picks = np.linspace(0, len(indices) - 1, count)
    chosen = np.unique(np.rint(picks).astype(np.int64))
    # Rounding collisions are possible on short inputs; fill deterministically.
    if len(chosen) < count:
        remaining = np.setdiff1d(np.arange(len(indices)), chosen, assume_unique=False)
        chosen = np.sort(np.concatenate([chosen, remaining[: count - len(chosen)]]))
    return indices[chosen[:count]]


# ---------------------------------------------------------------------------
# genotype summaries
# ---------------------------------------------------------------------------


def af_missingness_dosage(
    genotype: np.ndarray, train_idx: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Train-only ALT allele frequency, missingness and AF-imputed dosage."""
    train = genotype[train_idx]
    observed = train >= 0
    count = observed.sum(axis=0)
    af = np.divide(
        np.where(observed, train, 0).sum(axis=0),
        2.0 * count,
        out=np.full(genotype.shape[1], np.nan),
        where=count > 0,
    )
    missingness = 1.0 - count / float(len(train_idx))
    dosage = genotype.astype(np.float32)
    missing = dosage < 0
    dosage[missing] = np.broadcast_to((2.0 * af).astype(np.float32), dosage.shape)[missing]
    return af, missingness, dosage


# ---------------------------------------------------------------------------
# annotation-free gene pooling
# ---------------------------------------------------------------------------


def tss_distance_weights(positions: np.ndarray, tss_1based: int, kernel_bp: float) -> np.ndarray:
    """Exponential cis-distance kernel. Pure geometry: no eQTL, no annotation."""
    if kernel_bp <= 0:
        raise ValueError("kernel_bp must be positive")
    distance = np.abs(positions.astype(np.float64) - float(tss_1based))
    weights = np.exp(-distance / float(kernel_bp))
    total = float(weights.sum())
    if total <= 0:
        raise ValueError("degenerate distance kernel")
    return weights / total


def gene_pool_features(
    hidden: np.ndarray, positions: np.ndarray, tss_1based: int, kernel_bp: float
) -> np.ndarray:
    """Pool ``[n, blocks, snps, d]`` hidden states into one gene-level vector.

    Three annotation-free views are concatenated: a uniform cis mean, a
    TSS-distance-weighted mean, and a per-block mean. The routing weights are
    functions of genomic distance only, so a later burden-test comparison is not
    contaminated by prior knowledge of which genes matter.
    """
    n, blocks, snps, d = hidden.shape
    if len(positions) != blocks * snps:
        raise ValueError("positions must cover every SNP slot in the window")
    flat = hidden.reshape(n, blocks * snps, d).astype(np.float64)
    uniform = flat.mean(axis=1)
    weights = tss_distance_weights(positions, tss_1based, kernel_bp)
    weighted = np.einsum("nsd,s->nd", flat, weights)
    block_mean = hidden.astype(np.float64).mean(axis=2).reshape(n, blocks * d)
    return np.concatenate([uniform, weighted, block_mean], axis=1).astype(np.float32)


# ---------------------------------------------------------------------------
# ridge with a frozen alpha rule
# ---------------------------------------------------------------------------


def fit_ridge(
    x: np.ndarray,
    y_train: np.ndarray,
    train_idx: np.ndarray,
    val_idx: np.ndarray,
    alphas: list[float],
    y_val: np.ndarray,
):
    """Dual ridge; alpha is chosen only by validation MSE, ties take larger alpha.

    Returns a predictor closure over arbitrary standardized row sets plus the
    fitted scaling, so a perturbed design matrix can reuse the identical model.
    """
    x_train = x[train_idx].astype(np.float64)
    mean = x_train.mean(axis=0)
    sd = x_train.std(axis=0)
    keep = sd >= 1e-8
    if not keep.any():
        raise ValueError("every feature is constant on dev_train")
    mean_keep = mean[keep]
    sd_keep = sd[keep]
    z_train = (x_train[:, keep] - mean_keep) / sd_keep
    z_val = (x[val_idx].astype(np.float64)[:, keep] - mean_keep) / sd_keep

    y_mean = float(y_train.mean())
    y_sd = float(y_train.std(ddof=0))
    if y_sd <= 0:
        raise ValueError("training outcome has zero variance")
    y_train_std = (y_train - y_mean) / y_sd

    gram = z_train @ z_train.T
    val_cross = z_val @ z_train.T
    identity = np.eye(len(train_idx))
    best: tuple[float, float, np.ndarray] | None = None
    for alpha in alphas:
        dual = np.linalg.solve(gram + float(alpha) * identity, y_train_std)
        prediction = y_mean + y_sd * (val_cross @ dual)
        mse = float(np.mean((y_val - prediction) ** 2))
        if best is None or mse < best[0] - 1e-12 or (abs(mse - best[0]) <= 1e-12 and alpha > best[1]):
            best = (mse, float(alpha), dual)
    assert best is not None
    mse, alpha, dual = best

    def predict(rows: np.ndarray) -> np.ndarray:
        z = (rows.astype(np.float64)[:, keep] - mean_keep) / sd_keep
        return y_mean + y_sd * (z @ z_train.T @ dual)

    info = {
        "alpha": alpha,
        "validation_mse": mse,
        "input_features": int(x.shape[1]),
        "nonconstant_features": int(keep.sum()),
        "train_mean": y_mean,
    }
    return predict, info


def r2_against_train_mean(truth: np.ndarray, prediction: np.ndarray, train_mean: float) -> float:
    denominator = float(np.sum((truth - train_mean) ** 2))
    if denominator <= 0:
        raise ValueError("nonpositive R2 denominator")
    return 1.0 - float(np.sum((truth - prediction) ** 2)) / denominator


# ---------------------------------------------------------------------------
# two-way clustered bootstrap over individuals and genes
# ---------------------------------------------------------------------------


def two_way_clustered_bootstrap(
    truth: np.ndarray,
    predictions: dict[str, np.ndarray],
    train_means: np.ndarray,
    contrast: tuple[str, str],
    seed: int,
    replicates: int,
    batch: int = 250,
) -> dict[str, object]:
    """Resample evaluation individuals AND genes; genes are the replication unit.

    ``truth`` and each ``predictions[arm]`` are ``[n_genes, n_individuals]``.
    """
    n_genes, n_individuals = truth.shape
    for arm, value in predictions.items():
        if value.shape != truth.shape:
            raise ValueError(f"prediction shape mismatch for arm {arm}")
    if train_means.shape != (n_genes,):
        raise ValueError("train_means must be one value per gene")

    centred = truth - train_means[:, None]
    point: dict[str, float] = {}
    per_gene: dict[str, np.ndarray] = {}
    for arm, value in predictions.items():
        denominator = np.sum(centred**2, axis=1)
        if np.any(denominator <= 0):
            raise ValueError("nonpositive R2 denominator")
        gene_r2 = 1.0 - np.sum((truth - value) ** 2, axis=1) / denominator
        per_gene[arm] = gene_r2
        point[arm] = float(gene_r2.mean())

    rng = np.random.default_rng(seed)
    draws = {arm: np.zeros(replicates) for arm in predictions}
    done = 0
    while done < replicates:
        take = min(batch, replicates - done)
        gene_sample = rng.integers(0, n_genes, size=(take, n_genes))
        individual_sample = rng.integers(0, n_individuals, size=(take, n_individuals))
        yb = truth[gene_sample[:, :, None], individual_sample[:, None, :]]
        mb = train_means[gene_sample]
        denominator = np.sum((yb - mb[:, :, None]) ** 2, axis=2)
        valid = denominator > 0
        for arm, value in predictions.items():
            pb = value[gene_sample[:, :, None], individual_sample[:, None, :]]
            gene_r2 = np.where(
                valid, 1.0 - np.sum((yb - pb) ** 2, axis=2) / np.where(valid, denominator, 1.0), np.nan
            )
            draws[arm][done : done + take] = np.nanmean(gene_r2, axis=1)
        done += take

    left, right = contrast
    delta_draws = draws[left] - draws[right]
    delta = point[left] - point[right]
    low = float(np.quantile(delta_draws, 0.025))
    high = float(np.quantile(delta_draws, 0.975))
    return {
        "n_genes": int(n_genes),
        "n_individuals": int(n_individuals),
        "macro_r2": point,
        "per_gene_r2": {arm: [float(v) for v in values] for arm, values in per_gene.items()},
        "contrast": f"{left}-{right}",
        "macro_delta_r2": float(delta),
        "paired_bootstrap_ci95": [low, high],
        "pass": bool(delta > 0 and low > 0),
        "bootstrap_replicates": int(replicates),
        "bootstrap_seed": int(seed),
        "resampling_units": ["evaluation_individual", "gene"],
    }


# ---------------------------------------------------------------------------
# rank statistics and the matched null
# ---------------------------------------------------------------------------


def rankdata_average(values: np.ndarray) -> np.ndarray:
    """Average-tie ranks; avoids a scipy dependency in the frozen runner."""
    values = np.asarray(values, dtype=np.float64)
    order = np.argsort(values, kind="stable")
    ranks = np.empty(len(values), dtype=np.float64)
    sorted_values = values[order]
    i = 0
    while i < len(values):
        j = i
        while j + 1 < len(values) and sorted_values[j + 1] == sorted_values[i]:
            j += 1
        ranks[order[i : j + 1]] = 0.5 * (i + j) + 1.0
        i = j + 1
    return ranks


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    ra = rankdata_average(a)
    rb = rankdata_average(b)
    ra = ra - ra.mean()
    rb = rb - rb.mean()
    denominator = math.sqrt(float(np.sum(ra * ra) * np.sum(rb * rb)))
    if denominator <= 0:
        return 0.0
    return float(np.sum(ra * rb) / denominator)


def decile_strata(values: np.ndarray, bins: int = 10) -> np.ndarray:
    ranks = rankdata_average(values)
    return np.minimum((ranks - 1) / len(values) * bins, bins - 1).astype(np.int64)


def matched_null_spearman(
    score: np.ndarray,
    target: np.ndarray,
    strata: np.ndarray,
    seed: int,
    replicates: int,
) -> dict[str, object]:
    """Permute the gene score only within matched strata.

    Strata hold cis SNP count, MAF spectrum and expression variance roughly
    fixed, so a positive result cannot be produced by those nuisance variables
    alone.
    """
    observed = spearman(score, target)
    rng = np.random.default_rng(seed)
    groups = [np.flatnonzero(strata == s) for s in np.unique(strata)]
    null = np.zeros(replicates)
    for r in range(replicates):
        permuted = score.copy()
        for group in groups:
            if len(group) > 1:
                permuted[group] = score[rng.permutation(group)]
        null[r] = spearman(permuted, target)
    greater = float((np.sum(null >= observed) + 1) / (replicates + 1))
    two_sided = float((np.sum(np.abs(null) >= abs(observed)) + 1) / (replicates + 1))
    return {
        "observed_spearman": observed,
        "matched_null_mean": float(null.mean()),
        "matched_null_sd": float(null.std(ddof=1)),
        "p_one_sided_greater": greater,
        "p_two_sided": two_sided,
        "permutation_replicates": int(replicates),
        "permutation_seed": int(seed),
        "n_strata": int(len(groups)),
    }


def paired_gene_bootstrap_difference(
    score_a: np.ndarray,
    score_b: np.ndarray,
    target: np.ndarray,
    seed: int,
    replicates: int,
) -> dict[str, object]:
    """Bootstrap genes to bound ``spearman(a, T) - spearman(b, T)``."""
    observed = spearman(score_a, target) - spearman(score_b, target)
    rng = np.random.default_rng(seed)
    n = len(target)
    draws = np.zeros(replicates)
    for r in range(replicates):
        pick = rng.integers(0, n, size=n)
        draws[r] = spearman(score_a[pick], target[pick]) - spearman(score_b[pick], target[pick])
    low = float(np.quantile(draws, 0.025))
    high = float(np.quantile(draws, 0.975))
    return {
        "observed_difference": float(observed),
        "bootstrap_ci95": [low, high],
        "pass": bool(observed > 0 and low > 0),
        "bootstrap_replicates": int(replicates),
        "bootstrap_seed": int(seed),
    }
