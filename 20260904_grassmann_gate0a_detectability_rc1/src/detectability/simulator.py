"""Detectability planning simulator for Gate 0A.

Purpose
-------
Fix the number of simulation replicates R needed to reliably detect the FIXED
primary margin Delta R2_genetic = 0.005, using the SAME decision procedure Gate
0A will use: the inferential unit is the simulation replicate, and a regime is
declared to have positive headroom iff the percentile-bootstrap 95% CI of the
paired difference Delta_r = R2(A) - R2(B), taken over R replicates, has a LOWER
bound greater than 0 (a one-sided positive-headroom decision).

This simulator uses NO real genotype and makes NO biological claim. It models
the replicate-level paired difference as Delta_r ~ Normal(Delta_true, sigma^2),
where sigma is the replicate-to-replicate SD of the paired R2_genetic difference
(trait-architecture + estimation noise). Because sigma is not known until a real
panel is bound, it is a PLANNING AXIS: the simulator tabulates power/FPR across a
sigma grid, and a small real-panel pilot later pins sigma and reads off R.

Key identity (why this is cheap and exact)
------------------------------------------
For X_r = Delta + sigma * eps_r with eps_r ~ Normal(0,1), the percentile
bootstrap CI lower bound is  L(Delta, sigma) = Delta + sigma * L0,  where L0 is
the CI lower bound computed on the centered unit sample eps_r. The sample sd is
invariant to the shift Delta. Hence the decision  L(Delta,sigma) > 0  is exactly
 L0 > -Delta / sigma. We therefore simulate the unit lower bound L0 once per R
and obtain every (Delta, sigma) cell in closed form.
"""
from __future__ import annotations

import numpy as np


def unit_ci_lower_bounds(R, M, B, ci=0.95, rng=None, mem_budget=2e7):
    """Return M Monte-Carlo samples of the percentile-bootstrap CI lower bound
    for the mean of R i.i.d. Normal(0,1) draws.

    Parameters
    ----------
    R : int    number of replicates (the quantity Gate 0A must choose)
    M : int    outer Monte-Carlo draws
    B : int    bootstrap resamples per outer draw
    ci : float two-sided CI level (0.95 -> lower percentile at 2.5%)
    rng : np.random.Generator or int or None
    mem_budget : approx max number of ints held per chunk (controls chunking)
    """
    rng = np.random.default_rng(rng)
    alpha = (1.0 - ci) / 2.0
    out = np.empty(M, dtype=np.float64)
    chunk = max(1, int(mem_budget // (B * R)))
    i = 0
    while i < M:
        m = min(chunk, M - i)
        eps = rng.standard_normal((m, R))
        idx = rng.integers(0, R, size=(m, B, R))
        rows = np.arange(m)[:, None, None]
        boot_means = eps[rows, idx].mean(axis=2)  # (m, B)
        out[i:i + m] = np.quantile(boot_means, alpha, axis=1)
        i += m
    return out


def power_table(delta_grid, sigma_grid, R_grid, M, B, ci=0.95, seed=0):
    """Compute P(declare positive headroom) for every (delta, sigma, R) cell.

    Returns
    -------
    table : dict[(delta, sigma, R)] -> power (float)
    L0 : dict[R] -> np.ndarray of unit CI lower bounds (M,)
    """
    rng = np.random.default_rng(seed)
    L0 = {int(R): unit_ci_lower_bounds(int(R), M, B, ci, rng) for R in R_grid}
    table = {}
    for sigma in sigma_grid:
        for R in R_grid:
            low = L0[int(R)]
            for delta in delta_grid:
                thr = -float(delta) / float(sigma)
                table[(float(delta), float(sigma), int(R))] = float(np.mean(low > thr))
    return table, L0


def min_replicates(table, sigma, R_grid, target_delta, target_power):
    """Smallest R in R_grid reaching >= target_power at target_delta for a given
    sigma; returns None if no gridded R suffices."""
    for R in sorted(int(r) for r in R_grid):
        if table[(float(target_delta), float(sigma), R)] >= target_power:
            return R
    return None
