"""Representation arms and the adaptation baselines they have to beat.

The arms form a ladder of decreasing invariance:

    grassmann       vech(P)              invariant under GL(2)  <- the claim
    spectrum        vech(P) + (s1, s2)   invariant under O(2)
    second_moment   vech(E E.T)          invariant under O(2), richer than P
    bilinear_o2     tr(A_k E E.T)        invariant under O(2), matched capacity
    bilinear_raw    random quadratic     no invariance, matched capacity
    aligned         vec(E)               no invariance
    hybrid          grassmann + spectrum + aligned

and two adaptation baselines that get something the Grassmann arm does not: a
calibration set of pairs measured in *both* contexts.

    procrustes      vec(Q E), Q in O(r), fitted on calibration pairs
    gl_transport    vec(W E), W in R^{r x r} ridge, fitted on calibration pairs

Both are fitted on calibration pairs and applied to held-out test pairs.  Fitting
M on the test pair's own target data is not a baseline: least squares puts that
residual exactly at the chordal subspace distance, so it *is* the Grassmann arm
with the adaptation cost hidden.  ``geometry.gl2_alignment_residual`` states the
identity and the test suite asserts it.

Both transports act on the left, in trait space, because that is the only side a
map shared across pairs can act on.  Under the LD mechanism the nuisance is the
pair's own ``Lambda``, so it is right-acting and pair-specific, and no amount of
calibration data repairs it -- which is the real argument for an invariant
representation.  That argument is only honest if the sweep also contains the
context-wide nuisance where transport does win, so
``simulate.apply_distortion(shared=True)`` is part of the protocol.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "ARMS",
    "INVARIANCE",
    "batched_svd",
    "vech",
    "features",
    "fit_procrustes",
    "fit_gl_transport",
    "ridge_fit_predict",
    "transfer_r2",
]

ARMS = [
    "grassmann",
    "spectrum",
    "second_moment",
    "bilinear_o2",
    "bilinear_raw",
    "aligned",
    "procrustes",
    "gl_transport",
    "hybrid",
]

INVARIANCE = {
    "grassmann": "GL(2)",
    "spectrum": "O(2)",
    "second_moment": "O(2)",
    "bilinear_o2": "O(2)",
    "bilinear_raw": "none",
    "aligned": "none",
    "procrustes": "none (+calibration)",
    "gl_transport": "none (+calibration)",
    "hybrid": "mixed",
}


def batched_svd(E: np.ndarray):
    """``(n, r, 2) -> (U, s)`` with ``U`` the left singular vectors."""
    U, s, _ = np.linalg.svd(np.asarray(E, dtype=float), full_matrices=False)
    return U, s


def vech(M: np.ndarray) -> np.ndarray:
    """Lower-triangular half-vectorisation of a batch of symmetric matrices.

    Off-diagonal entries are scaled by ``sqrt(2)`` so the Euclidean norm of the
    result equals the Frobenius norm of the matrix; otherwise a ridge penalty
    would weight diagonal and off-diagonal information differently.
    """
    M = np.asarray(M, dtype=float)
    n_dim = M.shape[-1]
    rows, cols = np.tril_indices(n_dim)
    scale = np.where(rows == cols, 1.0, np.sqrt(2.0))
    return M[..., rows, cols] * scale


def _projector_and_spectrum(E: np.ndarray):
    U, s = batched_svd(E)
    P = np.einsum("nrk,nsk->nrs", U, U)
    return P, s


def features(arm: str, E: np.ndarray, rng: np.random.Generator, transport: np.ndarray | None = None,
             n_random: int | None = None) -> np.ndarray:
    """Feature matrix ``(n_pairs, n_features)`` for one arm."""
    E = np.asarray(E, dtype=float)
    n, r, k = E.shape

    if arm == "grassmann":
        P, _ = _projector_and_spectrum(E)
        return vech(P)
    if arm == "spectrum":
        P, s = _projector_and_spectrum(E)
        return np.hstack([vech(P), s])
    if arm == "second_moment":
        return vech(np.einsum("nrk,nsk->nrs", E, E))
    if arm in ("bilinear_o2", "bilinear_raw"):
        # Matched capacity: same feature count as the grassmann arm.
        m = n_random if n_random is not None else r * (r + 1) // 2
        if arm == "bilinear_o2":
            A = rng.standard_normal((m, r, r))
            G = np.einsum("nrk,nsk->nrs", E, E)
            return np.einsum("mrs,nrs->nm", A, G) / np.sqrt(r)
        V = E.reshape(n, r * k)
        a = rng.standard_normal((m, r * k))
        b = rng.standard_normal((m, r * k))
        return (V @ a.T) * (V @ b.T) / np.sqrt(r * k)
    if arm == "aligned":
        return E.reshape(n, r * k)
    if arm in ("procrustes", "gl_transport"):
        if transport is None:
            raise ValueError(f"arm {arm!r} needs a fitted transport matrix")
        return np.einsum("rs,nsk->nrk", transport, E).reshape(n, r * k)
    if arm == "hybrid":
        P, s = _projector_and_spectrum(E)
        return np.hstack([vech(P), s, E.reshape(n, r * k)])
    raise ValueError(f"unknown arm: {arm!r}")


def fit_procrustes(E_target: np.ndarray, E_source: np.ndarray) -> np.ndarray:
    """``argmin_{Q in O(r)} sum_p ||target_p - Q source_p||_F^2``."""
    C = np.einsum("nrk,nsk->rs", np.asarray(E_target, float), np.asarray(E_source, float))
    U, _, Vt = np.linalg.svd(C)
    return U @ Vt


def fit_gl_transport(E_target: np.ndarray, E_source: np.ndarray, alpha: float = 1e-3) -> np.ndarray:
    """Ridge-regularised ``argmin_W sum_p ||target_p - W source_p||_F^2``."""
    T = np.asarray(E_target, float)
    S = np.asarray(E_source, float)
    Cts = np.einsum("nrk,nsk->rs", T, S)
    Css = np.einsum("nrk,nsk->rs", S, S)
    scale = np.trace(Css) / Css.shape[0] if Css.shape[0] else 1.0
    return Cts @ np.linalg.inv(Css + alpha * scale * np.eye(Css.shape[0]))


def ridge_fit_predict(X_train, y_train, X_eval, alphas=(1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0, 1e3),
                      val_frac: float = 0.25, seed: int = 0):
    """Ridge with an internal validation split for ``alpha``; returns predictions.

    Columns are standardised on the training rows so one ``alpha`` means the same
    thing for every feature, and the intercept is not penalised.  Each arm picks
    its own ``alpha``: a shared one would hand an advantage to whichever arm
    happens to suit it.
    """
    X_train = np.asarray(X_train, float)
    y_train = np.asarray(y_train, float)
    X_eval = np.asarray(X_eval, float)

    rng = np.random.default_rng(seed)
    n = X_train.shape[0]
    perm = rng.permutation(n)
    n_val = max(1, int(round(val_frac * n)))
    val_idx, fit_idx = perm[:n_val], perm[n_val:]

    def solve(Xf, yf, alpha):
        mu, sd = Xf.mean(0), Xf.std(0)
        sd = np.where(sd > 0, sd, 1.0)
        Z = (Xf - mu) / sd
        ybar = yf.mean()
        A = Z.T @ Z + alpha * np.eye(Z.shape[1])
        w = np.linalg.solve(A, Z.T @ (yf - ybar))
        return lambda X: ((X - mu) / sd) @ w + ybar

    best_alpha, best_mse = None, np.inf
    for alpha in alphas:
        pred = solve(X_train[fit_idx], y_train[fit_idx], alpha)
        mse = float(np.mean((pred(X_train[val_idx]) - y_train[val_idx]) ** 2))
        if mse < best_mse:
            best_alpha, best_mse = alpha, mse

    final = solve(X_train, y_train, best_alpha)
    return final(X_eval), float(best_alpha)


def transfer_r2(y_eval, y_pred, y_train) -> float:
    """``1 - SSE / sum((y_eval - mean(y_train))^2)``.

    The training mean is the baseline, matching this project's existing
    convention, so a transfer that is worse than "predict the training mean"
    reads as negative rather than being flattered by the test set's own mean.
    """
    y_eval = np.asarray(y_eval, float)
    sse = float(np.sum((y_eval - np.asarray(y_pred, float)) ** 2))
    sst = float(np.sum((y_eval - np.mean(y_train)) ** 2))
    return 1.0 - sse / sst if sst > 0 else float("nan")
