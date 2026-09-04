"""Reference estimator for Gate 0A (used by the sigma pilot; Gate 0A will reuse it).

Scope of the sigma pilot: estimate ONLY the replicate-level SD of the paired
difference Delta_r = R2_genetic(A_test) - R2_genetic(B_pca_z). It never decides a
winner, the margin, the DGP set, k, arms, the success threshold, or which regime
is primary -- those are frozen upstream.

Key facts baked in here:
- Inputs are standardized per block with MAF-z using OUTER-TRAINING-FOLD stats only.
- Arms are fit on the training fold only; representations project all samples.
- Primary metric R2_genetic fits the downstream to the KNOWN genetic value g and
  evaluates R2(g, g_hat) on held-out folds. Because g is noise-free, this metric is
  heritability-independent, so the pilot needs no h2 grid.
- Interaction regimes use the matched bilinear cross-product head z_A (x) z_B.

Everything is numpy-only (no sklearn) so it runs under a minimal frozen python.
"""
from __future__ import annotations

import numpy as np


# --------------------------------------------------------------------------- #
# preprocessing
# --------------------------------------------------------------------------- #
def maf_z_stats(G_block_train):
    """MAF-z mean and sd from the training fold of one block (dosage 0/1/2)."""
    p = G_block_train.mean(axis=0) / 2.0
    var = 2.0 * p * (1.0 - p)
    sd = np.sqrt(np.where(var > 1e-12, var, 1.0))
    mean = 2.0 * p
    return mean, sd


def _sqdist(A, B):
    a2 = (A * A).sum(1)[:, None]
    b2 = (B * B).sum(1)[None, :]
    d2 = a2 + b2 - 2.0 * A @ B.T
    return np.maximum(d2, 0.0)


# --------------------------------------------------------------------------- #
# arms  (return dict: block_id -> representation array (N, k_b))
# --------------------------------------------------------------------------- #
def block_pca_z(G, blocks, k, tr, block_ids=None):
    reps = {}
    ids = block_ids if block_ids is not None else np.unique(blocks)
    for b in ids:
        idx = np.where(blocks == b)[0]
        Gb = G[:, idx]
        mean, sd = maf_z_stats(Gb[tr])
        Zb = (Gb - mean) / sd
        Ztr = Zb[tr]
        mu = Ztr.mean(0)
        _, _, Vt = np.linalg.svd(Ztr - mu, full_matrices=False)
        comp = Vt[: min(k, Vt.shape[0])].T
        reps[b] = (Zb - mu) @ comp
    return reps


def block_kpca_rbf(G, blocks, k, tr, block_ids=None):
    reps = {}
    ids = block_ids if block_ids is not None else np.unique(blocks)
    for b in ids:
        idx = np.where(blocks == b)[0]
        Gb = G[:, idx]
        mean, sd = maf_z_stats(Gb[tr])
        Zb = (Gb - mean) / sd
        Ztr = Zb[tr]
        n = Ztr.shape[0]
        d2tr = _sqdist(Ztr, Ztr)
        pos = d2tr[np.triu_indices(n, 1)]
        med = np.median(pos[pos > 0]) if np.any(pos > 0) else 1.0
        gamma = 1.0 / med
        Ktr = np.exp(-gamma * d2tr)
        one = np.full((n, n), 1.0 / n)
        Ktr_c = Ktr - one @ Ktr - Ktr @ one + one @ Ktr @ one
        vals, vecs = np.linalg.eigh((Ktr_c + Ktr_c.T) / 2.0)
        order = np.argsort(vals)[::-1]
        vals, vecs = vals[order], vecs[:, order]
        kk = min(k, n)
        vals_k, vecs_k = vals[:kk], vecs[:, :kk]
        alphas = np.zeros_like(vecs_k)
        good = vals_k > 1e-10
        alphas[:, good] = vecs_k[:, good] / np.sqrt(vals_k[good])
        Kall = np.exp(-gamma * _sqdist(Zb, Ztr))
        row = Ktr.mean(1)
        Kall_c = Kall - Kall.mean(1, keepdims=True) - row[None, :] + Ktr.mean()
        reps[b] = Kall_c @ alphas
    return reps


ARMS = {"pca_z": block_pca_z, "kpca_rbf": block_kpca_rbf}


# --------------------------------------------------------------------------- #
# downstream design + ridge
# --------------------------------------------------------------------------- #
def _cross(repA, repB):
    return (repA[:, :, None] * repB[:, None, :]).reshape(repA.shape[0], -1)


def design_matrix(reps, regime, sig_blocks):
    if regime in ("spectral_tail_adversarial", "major_LD_aligned"):
        return np.concatenate([reps[b] for b in sorted(reps)], axis=1)
    if regime in ("within_block_interaction",):
        b = sig_blocks[0]
        return _cross(reps[b], reps[b])
    if regime in ("between_block_interaction",):
        a, b = sig_blocks[0], sig_blocks[1]
        return _cross(reps[a], reps[b])
    raise ValueError(regime)


def _ridge_fit(X, y, lam):
    d = X.shape[1]
    A = X.T @ X + lam * np.eye(d)
    return np.linalg.solve(A, X.T @ y)


def ridge_cv_predict(Xtr, ytr, Xval, lambdas):
    # pick lambda by an inner 80/20 split of the training fold
    n = Xtr.shape[0]
    rs = np.random.default_rng(0)
    perm = rs.permutation(n)
    cut = int(0.8 * n)
    it, iv = perm[:cut], perm[cut:]
    best_lam, best_mse = lambdas[0], np.inf
    ytr_c = ytr - ytr.mean()
    for lam in lambdas:
        w = _ridge_fit(Xtr[it], ytr_c[it], lam)
        pred = Xtr[iv] @ w
        mse = float(np.mean((ytr_c[iv] - pred) ** 2))
        if mse < best_mse:
            best_mse, best_lam = mse, lam
    w = _ridge_fit(Xtr, ytr_c, best_lam)
    return Xval @ w + ytr.mean()


def r2(y, yhat):
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0


# --------------------------------------------------------------------------- #
# DGP: build the noise-free genetic value g for a regime (per replicate)
# --------------------------------------------------------------------------- #
def _block_dirs(G, blocks, b, tr, n_dir, which):
    idx = np.where(blocks == b)[0]
    Gb = G[:, idx]
    mean, sd = maf_z_stats(Gb[tr])
    Zb = (Gb - mean) / sd
    mu = Zb[tr].mean(0)
    _, _, Vt = np.linalg.svd(Zb[tr] - mu, full_matrices=False)
    m = Vt.shape[0]
    if which == "leading":
        dirs = Vt[:n_dir].T
    else:  # tail
        dirs = Vt[max(0, m - n_dir):].T
    scores = (Zb - mu) @ dirs  # (N, n_dir)
    return scores


def simulate_g(G, blocks, regime, sig_blocks, tr, rng, n_dir=3):
    if regime == "spectral_tail_adversarial":
        s = _block_dirs(G, blocks, sig_blocks[0], tr, n_dir, "tail")
        g = s @ rng.standard_normal(s.shape[1])
    elif regime == "major_LD_aligned":
        s = _block_dirs(G, blocks, sig_blocks[0], tr, n_dir, "leading")
        g = s @ rng.standard_normal(s.shape[1])
    elif regime == "within_block_interaction":
        s = _block_dirs(G, blocks, sig_blocks[0], tr, n_dir, "leading")
        s = (s - s[tr].mean(0)) / (s[tr].std(0) + 1e-12)
        g = (s[:, 0] * s[:, 1]) * rng.standard_normal()
    elif regime == "between_block_interaction":
        sa = _block_dirs(G, blocks, sig_blocks[0], tr, 1, "leading")[:, 0]
        sb = _block_dirs(G, blocks, sig_blocks[1], tr, 1, "leading")[:, 0]
        sa = (sa - sa[tr].mean()) / (sa[tr].std() + 1e-12)
        sb = (sb - sb[tr].mean()) / (sb[tr].std() + 1e-12)
        g = sa * sb * rng.standard_normal()
    else:
        raise ValueError(regime)
    g = g - g[tr].mean()
    s = g[tr].std()
    return g / s if s > 0 else g


# --------------------------------------------------------------------------- #
# one replicate: paired Delta of R2_genetic between two arms, via outer CV
# --------------------------------------------------------------------------- #
def replicate_delta(G, blocks, regime, sig_blocks, k, rng, arm_a="kpca_rbf",
                    arm_b="pca_z", folds=5, lambdas=(1e-2, 1e-1, 1.0, 10.0),
                    eval_block_ids=None):
    N = G.shape[0]
    perm = rng.permutation(N)
    fold_id = np.array_split(perm, folds)
    r2a, r2b = [], []
    for f in range(folds):
        val = fold_id[f]
        tr = np.concatenate([fold_id[j] for j in range(folds) if j != f])
        g = simulate_g(G, blocks, regime, sig_blocks, tr, rng)
        for arm, store in ((arm_a, r2a), (arm_b, r2b)):
            reps = ARMS[arm](G, blocks, k, tr, block_ids=eval_block_ids)
            X = design_matrix(reps, regime, sig_blocks)
            pred = ridge_cv_predict(X[tr], g[tr], X[val], list(lambdas))
            store.append(r2(g[val], pred))
    return float(np.mean(r2a) - np.mean(r2b))


def sigma_for_cell(G, blocks, regime, sig_blocks, k, R_pilot, seed,
                   folds=5, eval_block_ids=None):
    rng = np.random.default_rng(seed)
    deltas = [replicate_delta(G, blocks, regime, sig_blocks, k, rng, folds=folds,
                              eval_block_ids=eval_block_ids) for _ in range(R_pilot)]
    deltas = np.asarray(deltas)
    return {
        "mean_delta": float(deltas.mean()),
        "sigma_hat": float(deltas.std(ddof=1)),
        "R_pilot": int(R_pilot),
        "deltas": deltas.tolist(),
    }
