"""Level 0 estimands: standardised excess instability, and excess joint conservation.

Primary, as specified:

    I_metric = ( mean_p D_context(p) - mean_p D_null(p) ) / sd_p( D_null(p) )

Denominator is the spread of the null *across pairs*, not the standard error of
its mean, so ``I`` is an effect size and does not grow with the number of pairs.
The routing claim is ``I_sub << I_coord``.

Secondary, metric-free: each pair's context distance is placed in its own null
draws as a tail probability.  Probabilities are comparable between metrics
without any normalisation choice at all, so if the two summaries disagree the
normalisation of the distances is doing the talking and neither should be
believed.  Report both; pre-register the primary; never pick after the fact.

Third, and the one that can end the programme:  excess joint conservation.
``I_sub << I_coord`` on its own is compatible with a world where the Grassmannian
is decorative -- two individually stable vectors span a stable plane with no
joint structure whatsoever.  ``excess_joint_conservation`` asks whether the
observed plane is more stable than the single-SNP marginals already imply.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np

from .geometry import subspace_distance, coord_distance, vector_angle
from .nulls import deflation_surrogate

__all__ = [
    "InstabilityResult",
    "standardised_excess",
    "tail_quantiles",
    "paired_bootstrap_ci",
    "excess_joint_conservation",
    "orientation_diagnostic",
    "observed_distances",
    "permutation_reference",
    "conservation_fraction",
    "level0_verdict",
]

# Pre-registered absolute floor for the conservation fraction. A pair whose
# plane has gone most of the way to an unrelated pair's plane is not conserved,
# however badly its coordinates did. Fixed here, before any real data is seen.
CONSERVATION_FLOOR = 0.20


@dataclass
class InstabilityResult:
    metric: str
    mean_context: float
    mean_null: float
    sd_null: float
    index: float
    mean_tail_quantile: float

    def to_dict(self) -> dict:
        return asdict(self)


def observed_distances(EA, EB, metric: str = "chordal"):
    """Per-pair ``(subspace, coordinate)`` distances between two contexts."""
    EA = np.asarray(EA, dtype=float)
    EB = np.asarray(EB, dtype=float)
    sub = np.array([subspace_distance(EA[p], EB[p], metric=metric) for p in range(EA.shape[0])])
    coord = np.array([coord_distance(EA[p], EB[p]) for p in range(EA.shape[0])])
    return sub, coord


def permutation_reference(EA, EB, rng, metric: str = "chordal", n_draws: int = 20):
    """Distances between *mismatched* pairs: the "nothing survives" reference.

    Without it the comparison has no upper end. Two contexts where the locus
    simply does not transfer at all will still show the subspace moving less
    than the coordinates -- the subspace metric saturates earlier -- and that
    reads as conservation when it is the opposite. This supplies the scale on
    which "as unstable as an unrelated pair" is the zero.

    A derangement, not a plain shuffle: a pair must never be matched to itself.
    """
    EA = np.asarray(EA, dtype=float)
    EB = np.asarray(EB, dtype=float)
    n = EA.shape[0]
    if n < 2:
        raise ValueError("permutation reference needs at least two pairs")
    sub, coord = [], []
    for _ in range(n_draws):
        perm = rng.permutation(n)
        clash = perm == np.arange(n)
        if np.any(clash):  # rotate the offenders off their own index
            perm[clash] = np.roll(perm[clash], 1) if clash.sum() > 1 else (perm[clash] + 1) % n
        sub.append([subspace_distance(EA[p], EB[perm[p]], metric=metric) for p in range(n)])
        coord.append([coord_distance(EA[p], EB[perm[p]]) for p in range(n)])
    return np.array(sub).T, np.array(coord).T


def conservation_fraction(d_context, d_null, d_perm) -> np.ndarray:
    """Per-pair ``1 - (D_ctx - D_null) / (D_perm - D_null)``.

    1.0 means the pair moved no more than estimation noise forces; 0.0 means it
    moved as far as an unrelated pair. Both ends are measured in the metric's
    own units and then divided out, so the subspace and coordinate versions land
    on one absolute scale and can be compared directly -- which the raw
    distances, living on different spaces, never can.
    """
    d_context = np.asarray(d_context, dtype=float)
    null = np.nanmean(np.asarray(d_null, dtype=float), axis=1)
    perm = np.nanmean(np.asarray(d_perm, dtype=float), axis=1)
    span = perm - null
    with np.errstate(invalid="ignore", divide="ignore"):
        frac = 1.0 - (d_context - null) / span
    return np.where(np.abs(span) > 1e-12, frac, np.nan)


def standardised_excess(d_context: np.ndarray, d_null: np.ndarray, metric: str) -> InstabilityResult:
    """``I`` for one metric.  ``d_null`` is ``(n_pairs, n_draws)``."""
    d_context = np.asarray(d_context, dtype=float)
    d_null = np.asarray(d_null, dtype=float)
    per_pair_null = np.nanmean(d_null, axis=1)
    sd_null = float(np.nanstd(per_pair_null, ddof=1))
    mean_ctx = float(np.nanmean(d_context))
    mean_null = float(np.nanmean(per_pair_null))
    index = (mean_ctx - mean_null) / sd_null if sd_null > 0 else float("nan")
    return InstabilityResult(
        metric=metric,
        mean_context=mean_ctx,
        mean_null=mean_null,
        sd_null=sd_null,
        index=float(index),
        mean_tail_quantile=float(np.nanmean(tail_quantiles(d_context, d_null))),
    )


def tail_quantiles(d_context: np.ndarray, d_null: np.ndarray) -> np.ndarray:
    """Per-pair ``P(D_null >= D_context)``, add-one smoothed.

    Small values mean the context moved the pair further than noise alone can.
    Being a probability, this is directly comparable across metrics.
    """
    d_context = np.asarray(d_context, dtype=float)[:, None]
    d_null = np.asarray(d_null, dtype=float)
    n = d_null.shape[1]
    return (1.0 + np.nansum(d_null >= d_context, axis=1)) / (n + 1.0)


def paired_bootstrap_ci(
    values_a: np.ndarray,
    values_b: np.ndarray,
    n_boot: int,
    rng: np.random.Generator,
    alpha: float = 0.05,
):
    """Percentile CI for ``mean(a) - mean(b)``, resampling *pairs* jointly.

    Pairs are the inference unit; the two arms are resampled on the same index
    draw so the comparison stays paired.
    """
    a = np.asarray(values_a, dtype=float)
    b = np.asarray(values_b, dtype=float)
    n = a.size
    idx = rng.integers(0, n, size=(n_boot, n))
    diffs = np.nanmean(a[idx], axis=1) - np.nanmean(b[idx], axis=1)
    lo, hi = np.nanpercentile(diffs, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return float(np.nanmean(a) - np.nanmean(b)), float(lo), float(hi)


def excess_joint_conservation(
    EA,
    EB,
    n_draws: int,
    rng: np.random.Generator,
    metric: str = "chordal",
    n_boot: int = 2000,
) -> dict:
    """Is the pair's plane more stable than its single-SNP marginals imply?

    Positive ``excess`` means the observed plane moved *less* than surrogates
    that reproduce both columns' angular displacement and length but scramble
    the direction of that displacement.  The lower CI bound clearing zero is the
    condition for the pair-level, and therefore Grassmann, framing to be worth
    anything at all.
    """
    EA = np.asarray(EA, dtype=float)
    EB = np.asarray(EB, dtype=float)
    n_pairs = EA.shape[0]

    observed = np.empty(n_pairs)
    surrogate_mean = np.empty(n_pairs)
    quantile = np.empty(n_pairs)
    for p in range(n_pairs):
        observed[p] = subspace_distance(EA[p], EB[p], metric=metric)
        draws = deflation_surrogate(EA[p], EB[p], n_draws, rng, metric=metric)
        surrogate_mean[p] = float(np.nanmean(draws))
        quantile[p] = (1.0 + np.nansum(draws <= observed[p])) / (n_draws + 1.0)

    excess, lo, hi = paired_bootstrap_ci(surrogate_mean, observed, n_boot, rng)
    return {
        "metric": metric,
        "mean_observed": float(np.nanmean(observed)),
        "mean_surrogate": float(np.nanmean(surrogate_mean)),
        "excess": excess,
        "ci_low": lo,
        "ci_high": hi,
        "mean_lower_tail_quantile": float(np.nanmean(quantile)),
        "n_pairs": int(n_pairs),
        "n_draws": int(n_draws),
        "n_boot": int(n_boot),
        "verdict": "EXCESS_JOINT_CONSERVATION" if lo > 0 else "NO_EXCESS_JOINT_CONSERVATION",
    }


def orientation_diagnostic(EA, EB, threshold: float = 0.5) -> dict:
    """Flag columns whose effect vector nearly reverses between contexts.

    A flipped allele code is a sign flip, a sign flip is in ``GL(2)``, so it
    leaves ``span`` untouched while it maximises the coordinate distance.  That
    is a bookkeeping error which reads exactly like the hypothesis, and on a
    panel joined by position alone -- as this project's is -- it is a live risk,
    not a hypothetical.  A non-trivial ``flip_fraction`` blocks the run.
    """
    EA = np.asarray(EA, dtype=float)
    EB = np.asarray(EB, dtype=float)
    angles = np.array(
        [
            vector_angle(EA[p, :, c], EB[p, :, c])
            for p in range(EA.shape[0])
            for c in range(EA.shape[2])
        ]
    )
    flips = np.nansum(angles > np.pi * (1.0 - threshold / 2.0))
    n = int(np.sum(~np.isnan(angles)))
    frac = float(flips) / n if n else float("nan")
    return {
        "n_columns": n,
        "n_near_reversed": int(flips),
        "flip_fraction": frac,
        "median_angle_rad": float(np.nanmedian(angles)),
        "blocking": bool(frac > 0.02),
        "note": "flip_fraction > 0.02 blocks the run: fix allele orientation first",
    }


def level0_verdict(
    c_sub: float,
    c_coord: float,
    gap_ci_low: float,
    joint: dict,
    floor: float = CONSERVATION_FLOOR,
) -> dict:
    """Three independent conditions, and each failure means something different.

    ``absolute``  the plane is conserved at all: ``C_sub`` clears the floor.
                  Without this, a locus that transfers nothing still passes,
                  because the subspace metric saturates before the coordinate
                  one does.
    ``relative``  the plane is conserved *more than the coordinates are*. If the
                  coordinates survive too, there is no nuisance to quotient out
                  and the invariant representation is unmotivated -- that is a
                  different finding from the hypothesis failing.
    ``joint``     the plane's stability exceeds what the single-SNP marginals
                  already imply. Otherwise the Grassmannian is decorative.
    """
    absolute = c_sub > floor
    relative = gap_ci_low > 0
    jointly = joint["verdict"] == "EXCESS_JOINT_CONSERVATION"

    if not absolute:
        verdict = "NO_SUBSPACE_CONSERVATION_STOP"
    elif not relative:
        verdict = "NO_CONTEXT_DISTORTION_ROUTE_UNMOTIVATED"
    elif not jointly:
        verdict = "CONSERVED_BUT_DECORATIVE_STOP"
    else:
        verdict = "PROCEED_TO_LEVEL_1"

    return {
        "conservation_sub": float(c_sub),
        "conservation_coord": float(c_coord),
        "conservation_floor": float(floor),
        "gap_ci_low": float(gap_ci_low),
        "subspace_conserved_absolutely": bool(absolute),
        "subspace_more_conserved_than_coordinates": bool(relative),
        "excess_joint_conservation": bool(jointly),
        "verdict": verdict,
    }
