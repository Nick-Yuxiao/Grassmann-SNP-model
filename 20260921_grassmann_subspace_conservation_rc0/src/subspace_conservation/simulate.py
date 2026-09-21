"""Synthetic worlds.

Nothing here is evidence about genetics.  These generators exist so that the
Level 0 estimators can be calibrated against known truth, and so that the Level 2
phase surface has an axis definition before real genotype is ever touched.

Four worlds, and the fourth is the one that matters:

``S`` subspace-conserved  -- ``E_B = E_A M``, ``M in GL(2)``.  The hypothesis.
``C`` fully conserved     -- ``E_B = E_A``.  Coordinates survive too.
``N`` unrelated           -- independent draw.  Nothing survives.
``D`` decorative          -- each column rotated by *the same angle it moved in
      world S*, independently.  Single-SNP marginals identical to ``S``; joint
      structure destroyed.  An estimator that cannot separate ``D`` from ``S``
      cannot support a pair-level claim, and the calibration run exists mainly
      to show whether this one does.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .geometry import random_gl2, random_orthogonal, vector_angle
from .nulls import rotate_toward_random

__all__ = [
    "TruthSpec",
    "draw_truth",
    "apply_distortion",
    "apply_trait_distortion",
    "make_decorative_matched",
    "observe",
    "phenotype",
    "DISTORTION_LEVELS",
]

# Ordered by increasing distortion.  ``none`` and ``orth`` are not filler:
# right-O(2) leaves E E.T exactly invariant, so the spectrum and random-bilinear
# arms are invariant there too and the Grassmann arm has nothing to win.  Only
# from ``gl_*`` onward is the column space the maximal invariant.  A sweep that
# omits ``orth`` cannot tell those two claims apart.
DISTORTION_LEVELS = ["none", "orth", "gl_1.5", "gl_2", "gl_4", "gl_8", "gl_16"]


@dataclass
class TruthSpec:
    """Spectrum of the true pair effect matrix, in units of the noise SD."""

    r: int = 8
    sigma1: float = 8.0
    sigma2: float = 4.0
    noise: float = 1.0


def draw_truth(n_pairs: int, spec: TruthSpec, rng: np.random.Generator) -> np.ndarray:
    """``(n_pairs, r, 2)`` true effect matrices with the prescribed spectrum.

    Fixing the spectrum rather than drawing it is deliberate: ``sigma2 / noise``
    is the identifiability axis of the phase surface, so it has to be set, not
    observed.
    """
    out = np.empty((n_pairs, spec.r, 2))
    S = np.diag([spec.sigma1, spec.sigma2])
    for p in range(n_pairs):
        A = rng.standard_normal((spec.r, 2))
        U, _ = np.linalg.qr(A)
        out[p] = U @ S @ random_orthogonal(2, rng).T
    return out


def apply_distortion(E_true: np.ndarray, level: str, rng: np.random.Generator, shared: bool = False):
    """Re-express each pair in the second context's coordinates.

    ``shared=False`` draws one ``M`` per pair.  That is what the LD mechanism
    implies -- ``M`` comes from that pair's own tagging matrix ``Lambda`` -- and
    it is why a global left-side transport map cannot repair it.  ``shared=True``
    draws a single ``M`` for every pair, standing in for a context-wide nuisance
    such as a trait-scale difference, where transport *should* win.  Running only
    the first would rig the comparison; both are in the protocol.
    """
    E_true = np.asarray(E_true, dtype=float)
    n = E_true.shape[0]

    def draw_M():
        if level == "none":
            return np.eye(2)
        if level == "orth":
            return random_orthogonal(2, rng)
        if level.startswith("gl_"):
            return random_gl2(float(level.split("_", 1)[1]), rng)
        raise ValueError(f"unknown distortion level: {level!r}")

    if shared:
        M = draw_M()
        return E_true @ M
    return np.stack([E_true[p] @ draw_M() for p in range(n)])


def apply_trait_distortion(E_true: np.ndarray, strength: float, rng: np.random.Generator):
    """Left-acting, context-wide nuisance in trait space: ``E -> W E``.

    This is the anti-rigging control, and it is the mirror image of everything
    else in this module.  A right-acting nuisance mixes the two SNPs and leaves
    the column space alone, which is what the Grassmann arm is built for.  A
    left-acting one re-expresses the *traits* -- different measurement scales,
    different assay, a differently normalised panel between cohorts -- and it
    moves the column space like any other arm's input.  Grassmann cannot absorb
    it; a shared transport fitted on calibration pairs can.

    Without this slice in the sweep the comparison is stacked: every nuisance in
    the grid would be one the invariant representation is immune to by
    construction.  ``strength = 0`` is the identity.
    """
    E_true = np.asarray(E_true, dtype=float)
    r = E_true.shape[1]
    W = np.eye(r) + strength * rng.standard_normal((r, r)) / np.sqrt(r)
    return np.einsum("rs,nsk->nrk", W, E_true), W


def make_decorative_matched(E_A: np.ndarray, E_B: np.ndarray, rng: np.random.Generator):
    """World ``D``: reproduce ``E_B``'s single-SNP marginals, scramble the joint part."""
    E_A = np.asarray(E_A, dtype=float)
    E_B = np.asarray(E_B, dtype=float)
    out = np.empty_like(E_B)
    for p in range(E_A.shape[0]):
        for c in range(E_A.shape[2]):
            angle = vector_angle(E_A[p, :, c], E_B[p, :, c])
            norm = float(np.linalg.norm(E_B[p, :, c]))
            out[p, :, c] = rotate_toward_random(E_A[p, :, c], angle, norm, rng)
    return out


def observe(E_true: np.ndarray, noise: float, rng: np.random.Generator):
    """Add estimation noise; return ``(E_hat, SE)`` as a real study would report."""
    E_true = np.asarray(E_true, dtype=float)
    E_hat = E_true + rng.standard_normal(E_true.shape) * noise
    return E_hat, np.full_like(E_true, noise)


def phenotype(E_true: np.ndarray, gamma: float, reference_plane: np.ndarray, rng: np.random.Generator,
              noise_sd: float = 0.3):
    """``y = gamma * f(span) + (1 - gamma) * g(coords) + eps``, standardised.

    ``f`` depends on the pair's plane only -- its principal angles to a fixed
    disease-relevant plane.  ``g`` reads two raw coordinates, so it is destroyed
    by any column mixing.  ``gamma`` interpolates instead of switching, because
    the interesting quantity is where the crossover sits, not which of two
    rigged corners wins.

    ``y`` is defined once from the context-A truth and carried unchanged to
    context B: same biology, different coordinates.
    """
    from .geometry import subspace_distance

    E_true = np.asarray(E_true, dtype=float)
    n = E_true.shape[0]
    f = np.array([np.exp(-3.0 * subspace_distance(E_true[p], reference_plane) ** 2) for p in range(n)])
    g = E_true[:, 0, 0] - E_true[:, 1, 1]

    def z(v):
        s = np.std(v)
        return (v - np.mean(v)) / s if s > 0 else np.zeros_like(v)

    signal = gamma * z(f) + (1.0 - gamma) * z(g)
    return z(signal) + rng.standard_normal(n) * noise_sd
