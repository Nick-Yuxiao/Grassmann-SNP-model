"""Ridge-specification diagnostic: is ``B < A`` an artefact of the shared penalty?

Every arm comparison in this programme so far has been fitted the same way:
standardize every column of ``[covariates | genotype features]`` and solve one
dual ridge with a single alpha over the whole design matrix
(``run_stage2b.fit_ridge_predictions``, ``tqg1_core.fit_ridge``).  After
standardization each column contributes equally to the Gram matrix, so in arm B
the 12 covariate columns are 4.5% of the kernel and the 256 dosage columns are
95.5%.  ``(K + alpha I)^-1`` shrinks every direction together, so no alpha
keeps the covariate fit while controlling the dosage block, and arm B cannot
recover arm A even though arm B's design matrix contains arm A's.

This module implements the same estimand under three fits and nothing else
changes between them:

``S1``  every column standardized, one alpha over ``[covariates | dosage]``.
        The programme's current specification, reproduced here.
``S2``  covariates unpenalized, dosage penalized, fitted jointly.  Partial
        ridge, solved through Frisch-Waugh-Lovell.
``S3``  covariates fitted by OLS, dosage ridged on the residual.  The standard
        two-stage eQTL recipe.

Arm A is the covariate-only model under the same treatment: penalized ridge in
``S1``, OLS in ``S2`` and ``S3``.  Arm B adds the dosage block.  ``S2`` and
``S3`` are identical for arm A, which is a free consistency check.

Nothing here reads an outcome for a row it was not handed.
"""

from __future__ import annotations

import numpy as np

SPECIFICATIONS = ("S1", "S2", "S3")


def standardize(x: np.ndarray, fit_idx: np.ndarray) -> np.ndarray:
    """Centre and scale every column on ``fit_idx`` rows, dropping constants."""
    fitted = x[fit_idx].astype(np.float64)
    mean = fitted.mean(axis=0)
    sd = fitted.std(axis=0)
    keep = sd >= 1e-8
    if not keep.any():
        raise ValueError("every feature is constant on the fitting rows")
    return (x.astype(np.float64)[:, keep] - mean[keep]) / sd[keep]


def _pick(best, mse: float, alpha: float, payload):
    """Lowest validation MSE wins; ties take the larger alpha, as the programme does."""
    if best is None or mse < best[0] - 1e-12 or (abs(mse - best[0]) <= 1e-12 and alpha > best[1]):
        return (mse, float(alpha), payload)
    return best


def _ols(design: np.ndarray, y: np.ndarray) -> np.ndarray:
    coefficients, *_ = np.linalg.lstsq(design, y, rcond=None)
    return coefficients


def _fit_s1(zc, zd, y, fit_idx, tune_idx, alphas):
    z = zc if zd is None else np.column_stack([zc, zd])
    z_fit = z[fit_idx]
    y_fit = y[fit_idx]
    y_mean = float(y_fit.mean())
    y_sd = float(y_fit.std(ddof=0))
    if y_sd <= 0:
        raise ValueError("fitting outcome has zero variance")
    y_std = (y_fit - y_mean) / y_sd
    gram = z_fit @ z_fit.T
    identity = np.eye(len(fit_idx))
    cross = z[tune_idx] @ z_fit.T
    best = None
    for alpha in alphas:
        dual = np.linalg.solve(gram + float(alpha) * identity, y_std)
        prediction = y_mean + y_sd * (cross @ dual)
        best = _pick(best, float(np.mean((y[tune_idx] - prediction) ** 2)), alpha, dual)
    _, alpha, dual = best

    def predict(rows: np.ndarray) -> np.ndarray:
        return y_mean + y_sd * ((z[rows] @ z_fit.T) @ dual)

    return predict, alpha


def _fit_s2(zc, zd, y, fit_idx, tune_idx, alphas):
    design = np.column_stack([np.ones(len(y)), zc])
    design_fit = design[fit_idx]
    y_fit = y[fit_idx]
    if zd is None:
        beta = _ols(design_fit, y_fit)
        return (lambda rows: design[rows] @ beta), 0.0

    def residualize(values: np.ndarray) -> np.ndarray:
        return values - design_fit @ _ols(design_fit, values)

    dosage_fit = zd[fit_idx]
    y_residual = residualize(y_fit)
    dosage_residual = residualize(dosage_fit)
    gram = dosage_residual @ dosage_residual.T
    identity = np.eye(len(fit_idx))
    best = None
    for alpha in alphas:
        dual = np.linalg.solve(gram + float(alpha) * identity, y_residual)
        beta_dosage = dosage_residual.T @ dual
        beta_covariate = _ols(design_fit, y_fit - dosage_fit @ beta_dosage)
        prediction = design[tune_idx] @ beta_covariate + zd[tune_idx] @ beta_dosage
        best = _pick(
            best,
            float(np.mean((y[tune_idx] - prediction) ** 2)),
            alpha,
            (beta_covariate, beta_dosage),
        )
    _, alpha, (beta_covariate, beta_dosage) = best

    def predict(rows: np.ndarray) -> np.ndarray:
        return design[rows] @ beta_covariate + zd[rows] @ beta_dosage

    return predict, alpha


def _fit_s3(zc, zd, y, fit_idx, tune_idx, alphas):
    design = np.column_stack([np.ones(len(y)), zc])
    design_fit = design[fit_idx]
    y_fit = y[fit_idx]
    beta_covariate = _ols(design_fit, y_fit)
    if zd is None:
        return (lambda rows: design[rows] @ beta_covariate), 0.0

    residual = y_fit - design_fit @ beta_covariate
    dosage_fit = zd[fit_idx]
    gram = dosage_fit @ dosage_fit.T
    identity = np.eye(len(fit_idx))
    best = None
    for alpha in alphas:
        dual = np.linalg.solve(gram + float(alpha) * identity, residual)
        beta_dosage = dosage_fit.T @ dual
        prediction = design[tune_idx] @ beta_covariate + zd[tune_idx] @ beta_dosage
        best = _pick(best, float(np.mean((y[tune_idx] - prediction) ** 2)), alpha, beta_dosage)
    _, alpha, beta_dosage = best

    def predict(rows: np.ndarray) -> np.ndarray:
        return design[rows] @ beta_covariate + zd[rows] @ beta_dosage

    return predict, alpha


_FITTERS = {"S1": _fit_s1, "S2": _fit_s2, "S3": _fit_s3}


def fit_arm(
    specification: str,
    covariates: np.ndarray,
    dosage: np.ndarray | None,
    y: np.ndarray,
    fit_idx: np.ndarray,
    tune_idx: np.ndarray,
    alphas: list[float],
):
    """Fit one arm under one specification.

    Columns are standardized on ``fit_idx`` only; alpha is chosen only by MSE on
    ``tune_idx``.  Pass a single-element ``alphas`` to refit at a fixed alpha.
    Returns ``(predict, alpha)`` where ``predict`` takes row indices.
    """
    if specification not in _FITTERS:
        raise ValueError(f"unknown specification {specification!r}")
    if not alphas:
        raise ValueError("alphas must not be empty")
    zc = standardize(covariates, fit_idx)
    zd = None if dosage is None else standardize(dosage, fit_idx)
    return _FITTERS[specification](zc, zd, y, fit_idx, tune_idx, [float(a) for a in alphas])


def group_kfold(groups: list[str], n_folds: int, seed: int) -> np.ndarray:
    """Assign each row a fold, keeping every group whole. Deterministic in ``seed``."""
    if n_folds < 2:
        raise ValueError("n_folds must be at least 2")
    unique = sorted(set(groups))
    if len(unique) < n_folds:
        raise ValueError(f"{len(unique)} groups cannot fill {n_folds} folds")
    order = np.random.default_rng(seed).permutation(len(unique))
    assignment = {unique[int(order[i])]: i % n_folds for i in range(len(unique))}
    return np.asarray([assignment[g] for g in groups], dtype=np.int64)


def cross_validated_predictions(
    specification: str,
    covariates: np.ndarray,
    dosage: np.ndarray | None,
    y: np.ndarray,
    folds: np.ndarray,
    groups: list[str],
    alphas: list[float],
    inner_seed: int,
) -> tuple[np.ndarray, list[float]]:
    """Out-of-fold predictions for every row.

    Alpha is chosen inside each outer training set on a group-disjoint inner
    split, then the arm is refitted on the whole outer training set at that
    alpha.  No outer test row takes part in standardization, alpha selection or
    fitting.
    """
    predictions = np.full(len(y), np.nan)
    chosen: list[float] = []
    for fold in sorted(set(folds.tolist())):
        test_idx = np.flatnonzero(folds == fold)
        train_idx = np.flatnonzero(folds != fold)
        train_groups = [groups[i] for i in train_idx]
        inner = group_kfold(train_groups, 5, inner_seed + int(fold))
        inner_fit = train_idx[inner != 0]
        inner_tune = train_idx[inner == 0]
        _, alpha = fit_arm(
            specification, covariates, dosage, y, inner_fit, inner_tune, alphas
        )
        predict, _ = fit_arm(
            specification, covariates, dosage, y, train_idx, train_idx, [alpha]
        )
        predictions[test_idx] = predict(test_idx)
        chosen.append(float(alpha))
    if not np.isfinite(predictions).all():
        raise ValueError("some rows never landed in a test fold")
    return predictions, chosen
