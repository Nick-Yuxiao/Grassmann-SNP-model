"""Asset-free core of TQ-B1: the gene-effect bar.

Pure numpy. No torch, no GPU, no compiled extension. Everything here is
importable and unit-testable without any real cohort data.

The one question this package answers: how well does a purely additive cis
model's gene perturbation score already predict an independent WES burden
effect? That number is the bar any foundation model has to clear later.
"""

from __future__ import annotations

import hashlib
import math

import numpy as np

# ---------------------------------------------------------------------------
# PLINK 1 .bed reader (SNP-major), seek-based so a cis window costs one read
# ---------------------------------------------------------------------------

BED_MAGIC = b"\x6c\x1b\x01"
# PLINK1 2-bit codes -> dosage of A1 (the .bim column 5 allele). 0b01 is missing.
_BED_LOOKUP = np.array([2, -1, 1, 0], dtype=np.int8)


def bed_bytes_per_variant(n_samples: int) -> int:
    return (n_samples + 3) // 4


def read_bed_variants(
    path, n_samples: int, variant_indices: np.ndarray, n_variants: int
) -> np.ndarray:
    """Read the named variants as an ``[n_samples, len(variant_indices)]`` int8.

    Values are A1 dosage in {0, 1, 2} with -1 for missing. Variants are read by
    direct seek, so a 256-SNP cis window never touches the rest of the file.
    """
    stride = bed_bytes_per_variant(n_samples)
    out = np.empty((n_samples, len(variant_indices)), dtype=np.int8)
    with open(path, "rb") as handle:
        magic = handle.read(3)
        if magic != BED_MAGIC:
            raise ValueError(
                "not a SNP-major PLINK 1 .bed file (magic bytes "
                f"{magic!r}; individual-major files must be converted first)"
            )
        for column, variant in enumerate(variant_indices):
            v = int(variant)
            if not 0 <= v < n_variants:
                raise ValueError(f"variant index {v} out of range")
            handle.seek(3 + v * stride)
            raw = handle.read(stride)
            if len(raw) != stride:
                raise ValueError(f"truncated .bed at variant {v}")
            codes = np.unpackbits(np.frombuffer(raw, dtype=np.uint8), bitorder="little")
            codes = codes.reshape(-1, 2)
            packed = codes[:, 0] + 2 * codes[:, 1]
            out[:, column] = _BED_LOOKUP[packed[:n_samples]]
    return out


# ---------------------------------------------------------------------------
# deterministic, signal-free selection
# ---------------------------------------------------------------------------


def deterministic_order(keys: list[str], salt: str) -> list[str]:
    keyed = [(hashlib.sha256(f"{salt}:{k}".encode()).hexdigest(), k) for k in keys]
    keyed.sort()
    return [k for _, k in keyed]


def uniform_thin(indices: np.ndarray, count: int) -> np.ndarray:
    if count < 1:
        raise ValueError("count must be positive")
    if len(indices) < count:
        raise ValueError("cannot thin below the requested count")
    picks = np.rint(np.linspace(0, len(indices) - 1, count)).astype(np.int64)
    chosen = np.unique(picks)
    if len(chosen) < count:
        remaining = np.setdiff1d(np.arange(len(indices)), chosen)
        chosen = np.sort(np.concatenate([chosen, remaining[: count - len(chosen)]]))
    return indices[chosen[:count]]


def af_missingness_dosage(
    genotype: np.ndarray, train_idx: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Train-only allele frequency, missingness, and AF-imputed dosage."""
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
# primal ridge: O(n p^2), so n can be 500,000 while p stays a few hundred
# ---------------------------------------------------------------------------


class PrimalRidge:
    """Ridge by normal equations, fit once and evaluated at many alphas.

    The dual form used for small cohorts needs an n-by-n gram matrix and is
    unusable at UKB scale. Here the cost is driven by the feature count, not
    the sample count, so the same code runs on 500 or 500,000 individuals.
    """

    def __init__(
        self, x_train: np.ndarray, y_train: np.ndarray, penalty_mask: np.ndarray | None = None
    ) -> None:
        x = np.asarray(x_train, dtype=np.float64)
        y = np.asarray(y_train, dtype=np.float64)
        if len(x) != len(y):
            raise ValueError("design and outcome lengths differ")
        self.mean = x.mean(axis=0)
        sd = x.std(axis=0)
        self.keep = sd >= 1e-8
        if not self.keep.any():
            raise ValueError("every feature is constant on the training split")
        self.sd = np.where(self.keep, sd, 1.0)
        z = (x[:, self.keep] - self.mean[self.keep]) / self.sd[self.keep]
        self.y_mean = float(y.mean())
        self.y_sd = float(y.std(ddof=0))
        if self.y_sd <= 0:
            raise ValueError("training outcome has zero variance")
        y_std = (y - self.y_mean) / self.y_sd
        self.xtx = z.T @ z
        self.xty = z.T @ y_std
        self.p = int(self.keep.sum())
        # Covariates must be unpenalised. If one alpha shrinks covariates and
        # genotypes together, adding null dosage columns forces a larger alpha,
        # flattens the covariate coefficients, and makes the genetic arm lose to
        # the covariate-only arm for reasons that have nothing to do with
        # genetics. With a zero entry here, alpha -> infinity drives the
        # penalised coefficients to zero and the arm degenerates exactly onto
        # the covariate-only fit, so it can never be structurally worse.
        if penalty_mask is None:
            self.penalty = np.ones(self.p)
        else:
            mask = np.asarray(penalty_mask, dtype=np.float64)
            if mask.shape != (x.shape[1],):
                raise ValueError(
                    f"penalty_mask has {mask.shape} entries, expected one per input "
                    f"feature ({x.shape[1]})"
                )
            self.penalty = mask[self.keep]
        if np.any(self.penalty < 0):
            raise ValueError("penalty_mask entries must be non-negative")
        self.n = int(len(x))
        self._coefficients: dict[float, np.ndarray] = {}

    def coefficients(self, alpha: float) -> np.ndarray:
        key = float(alpha)
        if key not in self._coefficients:
            self._coefficients[key] = np.linalg.solve(
                self.xtx + key * np.diag(self.penalty), self.xty
            )
        return self._coefficients[key]

    def standardize(self, x: np.ndarray) -> np.ndarray:
        x = np.asarray(x, dtype=np.float64)
        return (x[:, self.keep] - self.mean[self.keep]) / self.sd[self.keep]

    def predict(self, x: np.ndarray, alpha: float) -> np.ndarray:
        return self.y_mean + self.y_sd * (self.standardize(x) @ self.coefficients(alpha))

    def select_alpha(self, x_val: np.ndarray, y_val: np.ndarray, alphas: list[float]):
        """Smallest validation MSE; ties take the larger alpha."""
        z = self.standardize(x_val)
        best: tuple[float, float] | None = None
        for alpha in alphas:
            prediction = self.y_mean + self.y_sd * (z @ self.coefficients(float(alpha)))
            mse = float(np.mean((np.asarray(y_val, dtype=np.float64) - prediction) ** 2))
            if best is None or mse < best[0] - 1e-12 or (abs(mse - best[0]) <= 1e-12 and alpha > best[1]):
                best = (mse, float(alpha))
        assert best is not None
        return best[1], best[0]


def r2_against_train_mean(truth: np.ndarray, prediction: np.ndarray, train_mean: float) -> float:
    denominator = float(np.sum((truth - train_mean) ** 2))
    if denominator <= 0:
        raise ValueError("nonpositive R2 denominator")
    return 1.0 - float(np.sum((truth - prediction) ** 2)) / denominator


def perturbation_score(
    model: PrimalRidge,
    design: np.ndarray,
    perturbed: np.ndarray,
    alpha: float,
) -> tuple[np.ndarray, float]:
    """Held-out prediction change when a gene's cis genotypes are removed.

    Columns that are identical in both designs cancel exactly, so the score
    isolates the cis contribution and never picks up covariate variance.
    """
    delta = model.predict(design, alpha) - model.predict(perturbed, alpha)
    if len(delta) < 2:
        raise ValueError("need at least two held-out individuals")
    return delta, float(np.std(delta, ddof=1))


# ---------------------------------------------------------------------------
# rank statistics, matched null, clustered bootstrap
# ---------------------------------------------------------------------------


def rankdata_average(values: np.ndarray) -> np.ndarray:
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
    ra = rankdata_average(a) - rankdata_average(a).mean()
    rb = rankdata_average(b) - rankdata_average(b).mean()
    denominator = math.sqrt(float(np.sum(ra * ra) * np.sum(rb * rb)))
    return 0.0 if denominator <= 0 else float(np.sum(ra * rb) / denominator)


def decile_strata(values: np.ndarray, bins: int = 10) -> np.ndarray:
    ranks = rankdata_average(values)
    return np.minimum((ranks - 1) / len(values) * bins, bins - 1).astype(np.int64)


def matched_null_spearman(
    score: np.ndarray, target: np.ndarray, strata: np.ndarray, seed: int, replicates: int
) -> dict[str, object]:
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
    return {
        "observed_spearman": observed,
        "matched_null_mean": float(null.mean()),
        "matched_null_sd": float(null.std(ddof=1)),
        "p_one_sided_greater": float((np.sum(null >= observed) + 1) / (replicates + 1)),
        "p_two_sided": float((np.sum(np.abs(null) >= abs(observed)) + 1) / (replicates + 1)),
        "permutation_replicates": int(replicates),
        "permutation_seed": int(seed),
        "n_strata": int(len(groups)),
    }


def gene_bootstrap_spearman(
    score: np.ndarray, target: np.ndarray, seed: int, replicates: int
) -> dict[str, object]:
    rng = np.random.default_rng(seed)
    n = len(target)
    draws = np.empty(replicates)
    for r in range(replicates):
        pick = rng.integers(0, n, size=n)
        draws[r] = spearman(score[pick], target[pick])
    return {
        "observed_spearman": spearman(score, target),
        "bootstrap_ci95": [float(np.quantile(draws, 0.025)), float(np.quantile(draws, 0.975))],
        "bootstrap_replicates": int(replicates),
        "bootstrap_seed": int(seed),
    }


def clustered_macro_bootstrap(
    truth: np.ndarray,
    predictions: dict[str, np.ndarray],
    train_means: np.ndarray,
    contrast: tuple[str, str],
    seed: int,
    replicates: int,
    batch: int = 200,
) -> dict[str, object]:
    """Macro R2 over genes with resampling of both individuals and genes."""
    n_genes, n_individuals = truth.shape
    for arm, value in predictions.items():
        if value.shape != truth.shape:
            raise ValueError(f"prediction shape mismatch for arm {arm}")
    centred = truth - train_means[:, None]
    denominator = np.sum(centred**2, axis=1)
    if np.any(denominator <= 0):
        raise ValueError("nonpositive R2 denominator")
    per_gene = {
        arm: 1.0 - np.sum((truth - value) ** 2, axis=1) / denominator
        for arm, value in predictions.items()
    }
    point = {arm: float(v.mean()) for arm, v in per_gene.items()}

    rng = np.random.default_rng(seed)
    draws = {arm: np.zeros(replicates) for arm in predictions}
    done = 0
    while done < replicates:
        take = min(batch, replicates - done)
        genes = rng.integers(0, n_genes, size=(take, n_genes))
        people = rng.integers(0, n_individuals, size=(take, n_individuals))
        yb = truth[genes[:, :, None], people[:, None, :]]
        mb = train_means[genes]
        denom = np.sum((yb - mb[:, :, None]) ** 2, axis=2)
        ok = denom > 0
        for arm, value in predictions.items():
            pb = value[genes[:, :, None], people[:, None, :]]
            r2 = np.where(ok, 1.0 - np.sum((yb - pb) ** 2, axis=2) / np.where(ok, denom, 1.0), np.nan)
            draws[arm][done : done + take] = np.nanmean(r2, axis=1)
        done += take

    left, right = contrast
    delta_draws = draws[left] - draws[right]
    delta = point[left] - point[right]
    low = float(np.quantile(delta_draws, 0.025))
    return {
        "n_genes": int(n_genes),
        "n_individuals": int(n_individuals),
        "macro_r2": point,
        "per_gene_r2": {arm: v for arm, v in per_gene.items()},
        "contrast": f"{left}-{right}",
        "macro_delta_r2": float(delta),
        "paired_bootstrap_ci95": [low, float(np.quantile(delta_draws, 0.975))],
        "pass": bool(delta > 0 and low > 0),
        "bootstrap_replicates": int(replicates),
        "bootstrap_seed": int(seed),
        "resampling_units": ["evaluation_individual", "gene"],
    }


def clustered_macro_bootstrap_shared(
    truth: np.ndarray,
    predictions: dict[str, np.ndarray],
    train_mean: float,
    contrast: tuple[str, str],
    seed: int,
    replicates: int,
    batch: int = 200,
) -> dict[str, object]:
    """Macro R2 over genes when every gene predicts the *same* outcome.

    In a cohort study the trait is one column, not one column per gene, so the
    truth vector must not be materialised once per gene: at UKB scale that
    alone would be gigabytes. ``truth`` is ``[n_individuals]`` and each
    ``predictions[arm]`` is ``[n_genes, n_individuals]``.
    """
    truth = np.asarray(truth, dtype=np.float64)
    n_individuals = len(truth)
    n_genes = next(iter(predictions.values())).shape[0]
    for arm, value in predictions.items():
        if value.shape != (n_genes, n_individuals):
            raise ValueError(f"prediction shape mismatch for arm {arm}")
    denominator = float(np.sum((truth - train_mean) ** 2))
    if denominator <= 0:
        raise ValueError("nonpositive R2 denominator")
    per_gene = {
        arm: 1.0 - np.sum((truth[None, :] - value) ** 2, axis=1) / denominator
        for arm, value in predictions.items()
    }
    point = {arm: float(v.mean()) for arm, v in per_gene.items()}

    rng = np.random.default_rng(seed)
    draws = {arm: np.zeros(replicates) for arm in predictions}
    done = 0
    while done < replicates:
        take = min(batch, replicates - done)
        genes = rng.integers(0, n_genes, size=(take, n_genes))
        people = rng.integers(0, n_individuals, size=(take, n_individuals))
        yb = truth[people]
        denom = np.sum((yb - train_mean) ** 2, axis=1)
        ok = denom > 0
        for arm, value in predictions.items():
            pb = value[genes[:, :, None], people[:, None, :]]
            r2 = 1.0 - np.sum((yb[:, None, :] - pb) ** 2, axis=2) / np.where(
                ok, denom, 1.0
            )[:, None]
            draws[arm][done : done + take] = np.where(ok, r2.mean(axis=1), np.nan)
        done += take

    left, right = contrast
    delta_draws = draws[left] - draws[right]
    delta = point[left] - point[right]
    low = float(np.nanquantile(delta_draws, 0.025))
    return {
        "n_genes": int(n_genes),
        "n_individuals": int(n_individuals),
        "macro_r2": point,
        "per_gene_r2": {arm: v for arm, v in per_gene.items()},
        "contrast": f"{left}-{right}",
        "macro_delta_r2": float(delta),
        "paired_bootstrap_ci95": [low, float(np.nanquantile(delta_draws, 0.975))],
        "pass": bool(delta > 0 and low > 0),
        "bootstrap_replicates": int(replicates),
        "bootstrap_seed": int(seed),
        "resampling_units": ["evaluation_individual", "gene"],
    }
