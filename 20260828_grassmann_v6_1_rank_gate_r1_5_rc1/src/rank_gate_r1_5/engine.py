from __future__ import annotations

from dataclasses import dataclass
import hashlib
import numpy as np

from grassmann_v6_1.core import conditional_matrices, fit_joint_ridge, geometry_scores
from grassmann_v6_1.null_engines import fit_common_subspace_null


@dataclass(frozen=True)
class RankGatedResult:
    observed: np.ndarray
    observed_raw: np.ndarray
    observed_gaps: np.ndarray
    observed_eligible: np.ndarray
    resampled: np.ndarray
    resampled_gaps: np.ndarray
    resampled_eligible: np.ndarray
    resampled_max: np.ndarray
    candidate_p_values: np.ndarray
    family_p_value: float
    multiplier_fingerprints: tuple[str, ...]


def rank_gated_direction(y, g, x, covariates, *, rank, ridge_lambda, minimum_gap):
    fit = fit_joint_ridge(y, g, x, covariates, ridge_lambda=ridge_lambda, whiten=True)
    geometry = geometry_scores(conditional_matrices(fit.B, fit.Gamma), rank)
    raw = float(geometry["direction_score"])
    gap = float(min(geometry["relative_rank_gaps"]))
    eligible = bool(gap >= minimum_gap)
    return (raw if eligible else 0.0), raw, gap, eligible


def _multipliers(subject_ids, rng):
    order = np.argsort(subject_ids)
    canonical = rng.choice((-1.0, 1.0), size=len(subject_ids))
    weights = np.empty(len(subject_ids), dtype=float)
    weights[order] = canonical
    fingerprint = hashlib.sha256(canonical.astype(np.int8).tobytes()).hexdigest()[:16]
    return weights, fingerprint


def run_rank_gated_maxT(family, *, resamples, seed, rank, ridge_lambda, minimum_gap):
    if minimum_gap < 0 or minimum_gap >= 1:
        raise ValueError("minimum_gap must be in [0,1)")
    if resamples < 19:
        raise ValueError("at least 19 resamples are required")
    k = len(family.regions)
    nulls = [fit_common_subspace_null(family.y, family.g, x, family.covariates,
             ridge_lambda=ridge_lambda, rank=rank) for x in family.regions]
    observed_parts = [rank_gated_direction(family.y, family.g, x, family.covariates,
                      rank=rank, ridge_lambda=ridge_lambda, minimum_gap=minimum_gap)
                      for x in family.regions]
    observed = np.array([x[0] for x in observed_parts])
    observed_raw = np.array([x[1] for x in observed_parts])
    observed_gaps = np.array([x[2] for x in observed_parts])
    observed_eligible = np.array([x[3] for x in observed_parts], dtype=bool)
    statistics = np.empty((k, resamples)); gaps = np.empty((k, resamples)); eligible = np.empty((k, resamples), dtype=bool)
    rng = np.random.default_rng(seed); fingerprints=[]
    for b in range(resamples):
        weights, fingerprint = _multipliers(family.subject_ids, rng); fingerprints.append(fingerprint)
        for j, (x, null) in enumerate(zip(family.regions, nulls)):
            y_star = null.fitted + null.residuals * weights[:, None]
            value, _, gap, ok = rank_gated_direction(y_star, family.g, x, family.covariates,
                rank=rank, ridge_lambda=ridge_lambda, minimum_gap=minimum_gap)
            statistics[j,b] = value; gaps[j,b] = gap; eligible[j,b] = ok
    resampled_max = statistics.max(axis=0); observed_max=float(observed.max())
    candidate_p=np.array([(1+np.count_nonzero(statistics[j]>=observed[j]))/(resamples+1) for j in range(k)])
    family_p=(1+np.count_nonzero(resampled_max>=observed_max))/(resamples+1)
    return RankGatedResult(observed, observed_raw, observed_gaps, observed_eligible,
        statistics, gaps, eligible, resampled_max, candidate_p, float(family_p), tuple(fingerprints))

