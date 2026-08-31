from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from grassmann_v6_1.core import conditional_matrices, fit_joint_ridge, geometry_scores, top_subspace


@dataclass(frozen=True)
class ScenarioData:
    y: np.ndarray
    g: np.ndarray
    x: np.ndarray
    covariates: np.ndarray
    B: np.ndarray
    Gamma: np.ndarray
    seed: int


@dataclass(frozen=True)
class FeasibilityFit:
    genotype_counts: tuple[int, int, int]
    minimum_group_count: int
    fitted_minimum_gap: float
    fitted_direction_score: float
    fitted_max_principal_angle_deg: float
    true_minimum_gap: float
    true_direction_score: float
    true_max_principal_angle_deg: float
    group_count_eligible: bool
    gap_eligible: bool
    jointly_estimable: bool
    gated_direction_score: float


def principal_angles_deg(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    singular_values = np.linalg.svd(left.T @ right, compute_uv=False)
    radians = np.arccos(np.clip(singular_values, -1.0, 1.0))
    return np.degrees(radians)


def truth_matrices(
    *,
    population_relative_gap: float,
    true_max_principal_angle_deg: float,
    effect_scale: float,
    region_features: int = 6,
    traits: int = 5,
) -> tuple[np.ndarray, np.ndarray]:
    if region_features < 4 or traits < 3:
        raise ValueError("truth construction requires at least four region features and three traits")
    if not 0 <= population_relative_gap < 1:
        raise ValueError("population_relative_gap must be in [0,1)")
    if not 0 <= true_max_principal_angle_deg <= 90:
        raise ValueError("principal angle must be in [0,90]")
    if effect_scale <= 0:
        raise ValueError("effect_scale must be positive")

    angle = math.radians(true_max_principal_angle_deg)
    u1_0 = np.eye(region_features)[:, 0]
    u1_2 = math.cos(angle) * np.eye(region_features)[:, 0] + math.sin(angle) * np.eye(region_features)[:, 2]
    u2 = np.eye(region_features)[:, 1]
    u3 = np.eye(region_features)[:, 3]
    v1 = np.eye(traits)[:, 0]
    v2 = np.eye(traits)[:, 1]
    v3 = np.eye(traits)[:, 2]

    singular_1 = effect_scale
    singular_2 = 0.70 * effect_scale
    singular_3 = singular_2 * (1.0 - population_relative_gap)
    m0 = (
        singular_1 * np.outer(u1_0, v1)
        + singular_2 * np.outer(u2, v2)
        + singular_3 * np.outer(u3, v3)
    )
    m2 = (
        singular_1 * np.outer(u1_2, v1)
        + singular_2 * np.outer(u2, v2)
        + singular_3 * np.outer(u3, v3)
    )
    return m0, (m2 - m0) / 2.0


def _ar1_covariance(width: int, rho: float) -> np.ndarray:
    indices = np.arange(width)
    return rho ** np.abs(indices[:, None] - indices[None, :])


def generate_scenario(
    *,
    seed: int,
    n: int,
    maf: float,
    population_relative_gap: float,
    true_max_principal_angle_deg: float,
    effect_scale: float,
    residual_sd: float,
    conditional_ld_rhos_by_dosage: tuple[float, float, float],
    residual_scales_by_dosage: tuple[float, float, float],
    region_features: int = 6,
    traits: int = 5,
    covariate_count: int = 4,
) -> ScenarioData:
    if not 0 < maf <= 0.5:
        raise ValueError("MAF must be in (0,0.5]")
    rng = np.random.default_rng(seed)
    g = rng.binomial(2, maf, size=n).astype(float)
    covariates = rng.normal(size=(n, covariate_count))
    x = np.empty((n, region_features), dtype=float)
    for dosage in (0, 1, 2):
        mask = g == dosage
        count = int(mask.sum())
        if count:
            covariance = _ar1_covariance(region_features, conditional_ld_rhos_by_dosage[dosage])
            x[mask] = rng.multivariate_normal(np.zeros(region_features), covariance, size=count)

    b, gamma = truth_matrices(
        population_relative_gap=population_relative_gap,
        true_max_principal_angle_deg=true_max_principal_angle_deg,
        effect_scale=effect_scale,
        region_features=region_features,
        traits=traits,
    )
    covariate_effect = rng.normal(scale=0.15, size=(covariate_count, traits))
    target_effect = np.linspace(0.10, 0.30, traits)
    mean = covariates @ covariate_effect + g[:, None] * target_effect
    mean += x @ b + (g[:, None] * x) @ gamma
    residual_scale = residual_sd * np.take(np.asarray(residual_scales_by_dosage), g.astype(int))
    y = mean + rng.normal(size=(n, traits)) * residual_scale[:, None]
    return ScenarioData(y, g, x, covariates, b, gamma, seed)


def fit_scenario(
    data: ScenarioData,
    *,
    rank: int,
    ridge_lambda: float,
    minimum_group_count: int,
    minimum_fitted_rank_gap: float,
) -> FeasibilityFit:
    fit = fit_joint_ridge(
        data.y,
        data.g,
        data.x,
        data.covariates,
        ridge_lambda=ridge_lambda,
        whiten=True,
    )
    fitted_matrices = conditional_matrices(fit.B, fit.Gamma)
    fitted_geometry = geometry_scores(fitted_matrices, rank)
    fitted_bases = [top_subspace(matrix, rank)[0] for matrix in fitted_matrices]
    fitted_angle = float(principal_angles_deg(fitted_bases[0], fitted_bases[2]).max())

    true_matrices = conditional_matrices(data.B, data.Gamma)
    true_geometry = geometry_scores(true_matrices, rank)
    true_bases = [top_subspace(matrix, rank)[0] for matrix in true_matrices]
    true_angle = float(principal_angles_deg(true_bases[0], true_bases[2]).max())

    counts = tuple(int(np.count_nonzero(data.g == dosage)) for dosage in (0, 1, 2))
    fitted_gap = float(min(fitted_geometry["relative_rank_gaps"]))
    group_ok = min(counts) >= minimum_group_count
    gap_ok = fitted_gap >= minimum_fitted_rank_gap
    estimable = bool(group_ok and gap_ok)
    return FeasibilityFit(
        genotype_counts=counts,
        minimum_group_count=min(counts),
        fitted_minimum_gap=fitted_gap,
        fitted_direction_score=float(fitted_geometry["direction_score"]),
        fitted_max_principal_angle_deg=fitted_angle,
        true_minimum_gap=float(min(true_geometry["relative_rank_gaps"])),
        true_direction_score=float(true_geometry["direction_score"]),
        true_max_principal_angle_deg=true_angle,
        group_count_eligible=bool(group_ok),
        gap_eligible=bool(gap_ok),
        jointly_estimable=estimable,
        gated_direction_score=float(fitted_geometry["direction_score"] if estimable else 0.0),
    )


def wilson_interval(successes: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if total <= 0 or not 0 <= successes <= total:
        raise ValueError("invalid binomial counts")
    proportion = successes / total
    denominator = 1.0 + z * z / total
    center = (proportion + z * z / (2.0 * total)) / denominator
    half = z * math.sqrt(proportion * (1.0 - proportion) / total + z * z / (4.0 * total * total)) / denominator
    return max(0.0, center - half), min(1.0, center + half)
