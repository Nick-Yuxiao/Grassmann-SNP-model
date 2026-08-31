"""Exploratory, non-evidentiary design-feasibility utilities."""

from .engine import (
    FeasibilityFit,
    generate_scenario,
    fit_scenario,
    principal_angles_deg,
    truth_matrices,
    wilson_interval,
)

__all__ = [
    "FeasibilityFit",
    "generate_scenario",
    "fit_scenario",
    "principal_angles_deg",
    "truth_matrices",
    "wilson_interval",
]
