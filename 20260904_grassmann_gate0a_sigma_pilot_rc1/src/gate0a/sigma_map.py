"""Map an estimated replicate SD sigma_hat to the required replicate count R.

Mirrors the detectability gate's decision procedure exactly (percentile bootstrap
CI over R replicates, one-sided lower bound > 0) via the shift/scale identity
power(Delta, sigma, R) = mean(L0_R > -Delta/sigma), where L0_R is the unit
bootstrap CI lower bound for the mean of R iid N(0,1) draws. The margin 0.005 and
the power target 0.90 are fixed upstream and never chosen here.
"""
from __future__ import annotations

import numpy as np

PRIMARY_MARGIN = 0.005
TARGET_POWER = 0.90


def _unit_lower_bounds(R, M, B, ci, rng):
    alpha = (1.0 - ci) / 2.0
    out = np.empty(M)
    chunk = max(1, int(2e7 // (B * R)))
    i = 0
    while i < M:
        m = min(chunk, M - i)
        eps = rng.standard_normal((m, R))
        idx = rng.integers(0, R, size=(m, B, R))
        rows = np.arange(m)[:, None, None]
        out[i:i + m] = np.quantile(eps[rows, idx].mean(axis=2), alpha, axis=1)
        i += m
    return out


def power_at(sigma, R, delta=PRIMARY_MARGIN, M=3000, B=1500, ci=0.95, seed=0):
    rng = np.random.default_rng(seed)
    low = _unit_lower_bounds(int(R), M, B, ci, rng)
    return float(np.mean(low > -float(delta) / float(sigma)))


def required_replicates(sigma, R_grid=(20, 50, 100, 200, 400, 800),
                        delta=PRIMARY_MARGIN, target_power=TARGET_POWER,
                        M=3000, B=1500, seed=0):
    """Smallest R in R_grid reaching >= target_power at the fixed margin; None if
    none suffice (then R must be increased, never the margin lowered)."""
    for R in sorted(int(r) for r in R_grid):
        if power_at(sigma, R, delta, M, B, seed=seed) >= target_power:
            return R
    return None
