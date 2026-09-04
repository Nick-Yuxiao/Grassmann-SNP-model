"""Reference estimator for Gate 0A (used by the sigma pilot; Gate 0A will reuse it).

Scope of the sigma pilot: estimate ONLY the replicate-level SD of the paired
difference Delta_r = R2_genetic(A_test) - R2_genetic(B_pca_z). It never decides a
winner, the margin, the DGP set, k, arms, the success threshold, or which regime
is primary -- those are frozen upstream.

Contract-faithful design (fixes over the first draft):
- Runs on the TRAINING population only; the 249 donor-validation rows are never
  touched by the pilot (the caller passes G already restricted to train_idx).
- The outer 5-fold split is drawn ONCE and shared byte-identically across all
  regimes, budgets, arms, and simulation replicates.
- A simulation replicate is ONE fixed genetic value g_r: g_r is drawn once and the
  same g_r is predicted across all five outer folds (g does not change per fold).
- Penalty selection is a real nested inner 5-fold CV.
- Each regime's g_r = named causal component + diffuse polygenic background.
- Arm representations do NOT depend on g, so they are computed once per
  (arm, fold, block) at k_max and cached; k=8 is the length-8 prefix of k=16.

numpy-only (no sklearn) so it runs under a minimal frozen python.
"""
from __future__ import annotations

import numpy as np

ADDITIVE = {"spectral_tail_adversarial", "major_LD_aligned"}
WITHIN_INT = {"within_block_interaction"}
BETWEEN_INT = {"between_block_interaction"}


# --------------------------------------------------------------------------- #
# preprocessing
# --------------------------------------------------------------------------- #
def maf_z_stats(G_block_train):
    p = G_block_train.mean(axis=0) / 2.0
    var = 2.0 * p * (1.0 - p)
    sd = np.sqrt(np.where(var > 1e-12, var, 1.0))
    return 2.0 * p, sd


def _sqdist(A, B):
    a2 = (A * A).sum(1)[:, None]
    b2 = (B * B).sum(1)[None, :]
    return np.maximum(a2 + b2 - 2.0 * A @ B.T, 0.0)


def _standardize(v, ref=None):
    ref = v if ref is None else ref
    s = ref.std()
    return (v - ref.mean()) / (s if s > 0 else 1.0)


# --------------------------------------------------------------------------- #
# folds (drawn once, shared everywhere)
# --------------------------------------------------------------------------- #
def build_folds(n, folds, seed):
    rng = np.random.default_rng(seed)
    parts = np.array_split(rng.permutation(n), folds)
    out = []
    for f in range(folds):
        val = parts[f]
        tr = np.concatenate([parts[j] for j in range(folds) if j != f])
        out.append((tr, val))
    return out


# --------------------------------------------------------------------------- #
# arms: fit on fold-train, project all training rows (return (N_train, k_max))
# --------------------------------------------------------------------------- #
def _pca_z_project(Gb, tr, k):
    mean, sd = maf_z_stats(Gb[tr])
    Zb = (Gb - mean) / sd
    mu = Zb[tr].mean(0)
    _, _, Vt = np.linalg.svd(Zb[tr] - mu, full_matrices=False)
    comp = Vt[: min(k, Vt.shape[0])].T
    return (Zb - mu) @ comp


def _kpca_rbf_project(Gb, tr, k):
    mean, sd = maf_z_stats(Gb[tr])
    Zb = (Gb - mean) / sd
    Ztr = Zb[tr]
    n = Ztr.shape[0]
    d2tr = _sqdist(Ztr, Ztr)
    off = d2tr[np.triu_indices(n, 1)]
    med = np.median(off[off > 0]) if np.any(off > 0) else 1.0
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
    return Kall_c @ alphas


_ARM_FN = {"pca_z": _pca_z_project, "kpca_rbf": _kpca_rbf_project}


def precompute_reps(G, blocks, folds, k_max, block_ids, arms=("pca_z", "kpca_rbf")):
    """cache[(arm, fold_index, block)] = representation (N_train, <=k_max).

    Independent of g and regime, so computed once and reused across all
    replicates, regimes, and budgets. k=8 is the prefix cache[...][:, :8]."""
    cache = {}
    for arm in arms:
        fn = _ARM_FN[arm]
        for fi, (tr, _) in enumerate(folds):
            for b in block_ids:
                idx = np.where(blocks == b)[0]
                cache[(arm, fi, b)] = fn(G[:, idx], tr, k_max)
    return cache


# --------------------------------------------------------------------------- #
# downstream design + nested-CV ridge
# --------------------------------------------------------------------------- #
def _cross(a, b):
    return (a[:, :, None] * b[:, None, :]).reshape(a.shape[0], -1)


def design(cache, arm, fi, regime, sig_blocks, eval_ids, k):
    if regime in ADDITIVE:
        return np.concatenate([cache[(arm, fi, b)][:, :k] for b in sorted(eval_ids)], axis=1)
    if regime in WITHIN_INT:
        r = cache[(arm, fi, sig_blocks[0])][:, :k]
        return _cross(r, r)
    if regime in BETWEEN_INT:
        a = cache[(arm, fi, sig_blocks[0])][:, :k]
        b = cache[(arm, fi, sig_blocks[1])][:, :k]
        return _cross(a, b)
    raise ValueError(regime)


def _ridge_fit(X, y, lam):
    d = X.shape[1]
    return np.linalg.solve(X.T @ X + lam * np.eye(d), X.T @ y)


def ridge_nested_predict(Xtr, ytr, Xval, lambdas, inner_folds, seed):
    n = len(ytr)
    parts = np.array_split(np.random.default_rng(seed).permutation(n), inner_folds)
    ym = ytr.mean()
    ytr_c = ytr - ym
    mse = np.zeros(len(lambdas))
    for f in range(inner_folds):
        iv = parts[f]
        it = np.concatenate([parts[j] for j in range(inner_folds) if j != f])
        for li, lam in enumerate(lambdas):
            w = _ridge_fit(Xtr[it], ytr_c[it], lam)
            mse[li] += float(np.sum((ytr_c[iv] - Xtr[iv] @ w) ** 2))
    lam = lambdas[int(np.argmin(mse))]
    w = _ridge_fit(Xtr, ytr_c, lam)
    return Xval @ w + ym


def r2(y, yhat):
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    return 1.0 - float(np.sum((y - yhat) ** 2)) / ss_tot if ss_tot > 0 else 0.0


# --------------------------------------------------------------------------- #
# DGP: one fixed genetic value g per replicate (named component + polygenic bg)
# --------------------------------------------------------------------------- #
def _dirs(G, blocks, b, n_dir, which):
    idx = np.where(blocks == b)[0]
    Gb = G[:, idx]
    mean, sd = maf_z_stats(Gb)          # population (train) stats -- g is the truth
    Zb = (Gb - mean) / sd
    mu = Zb.mean(0)
    _, _, Vt = np.linalg.svd(Zb - mu, full_matrices=False)
    m = Vt.shape[0]
    dirs = Vt[:n_dir].T if which == "leading" else Vt[max(0, m - n_dir):].T
    return (Zb - mu) @ dirs


def _polygenic(G, blocks, eval_ids, rng):
    cols = []
    for b in sorted(eval_ids):
        idx = np.where(blocks == b)[0]
        Gb = G[:, idx]
        mean, sd = maf_z_stats(Gb)
        cols.append((Gb - mean) / sd)
    Z = np.concatenate(cols, axis=1)
    beta = rng.standard_normal(Z.shape[1])
    return Z @ beta


def simulate_g(G, blocks, regime, sig_blocks, eval_ids, rng, n_dir, poly_frac):
    if regime == "spectral_tail_adversarial":
        s = _dirs(G, blocks, sig_blocks[0], n_dir, "tail")
        named = s @ rng.standard_normal(s.shape[1])
    elif regime == "major_LD_aligned":
        s = _dirs(G, blocks, sig_blocks[0], n_dir, "leading")
        named = s @ rng.standard_normal(s.shape[1])
    elif regime == "within_block_interaction":
        s = _dirs(G, blocks, sig_blocks[0], 2, "leading")
        s0, s1 = _standardize(s[:, 0]), _standardize(s[:, 1])
        named = s0 * s1 * rng.standard_normal()
    elif regime == "between_block_interaction":
        sa = _standardize(_dirs(G, blocks, sig_blocks[0], 1, "leading")[:, 0])
        sb = _standardize(_dirs(G, blocks, sig_blocks[1], 1, "leading")[:, 0])
        named = sa * sb * rng.standard_normal()
    else:
        raise ValueError(regime)
    named = _standardize(named)
    poly = _standardize(_polygenic(G, blocks, eval_ids, rng))
    g = np.sqrt(1.0 - poly_frac) * named + np.sqrt(poly_frac) * poly
    return _standardize(g)


# --------------------------------------------------------------------------- #
# per-cell sigma
# --------------------------------------------------------------------------- #
def sigma_for_cell(cache, G, blocks, folds, regime, sig_blocks, eval_ids, k,
                   R_pilot, seed, lambdas=(1e-2, 1e-1, 1.0, 10.0),
                   inner_folds=5, n_dir=3, poly_frac=0.5, inner_seed=12345):
    rng = np.random.default_rng(seed)
    deltas = []
    for _ in range(R_pilot):
        g = simulate_g(G, blocks, regime, sig_blocks, eval_ids, rng, n_dir, poly_frac)
        arm_r2 = {}
        for arm in ("kpca_rbf", "pca_z"):
            fr = []
            for fi, (tr, val) in enumerate(folds):
                X = design(cache, arm, fi, regime, sig_blocks, eval_ids, k)
                pred = ridge_nested_predict(X[tr], g[tr], X[val], list(lambdas),
                                            inner_folds, inner_seed + fi)
                fr.append(r2(g[val], pred))
            arm_r2[arm] = float(np.mean(fr))
        deltas.append(arm_r2["kpca_rbf"] - arm_r2["pca_z"])
    deltas = np.asarray(deltas)
    return {
        "mean_delta_NON_EVIDENCE_DIAGNOSTIC_ONLY": float(deltas.mean()),
        "sigma_hat": float(deltas.std(ddof=1)),
        "R_pilot": int(R_pilot),
        "deltas": deltas.tolist(),
    }
