"""Null distributions for context instability.

The whole Level 0 estimand is a difference against a null, because the raw
numbers are not comparable:  ``subspace_distance`` and ``coord_distance`` are
different metrics on different spaces and reading one against the other tells
you nothing.  What *is* comparable is how far each sits above the instability
that pure estimation noise already produces in that same metric.

Two nulls, because two kinds of input:

``parametric``  needs only ``(beta, se)`` summary statistics.  This is the one
                the real Level 0 run should use: no individual-level row ever
                leaves the server, which is the same data-governance boundary
                the Concept Gate package already operates under.

``split_half``  needs individual-level data and is provided for the calibration
                study, where the truth is known.  Its one trap is sample size:
                two halves of one context are each estimated on n/2, so the null
                is noisier than a full-n cross-context comparison and the excess
                comes out biased downward.  Contexts must be sub-sampled to the
                same n/2.  ``matched_n`` records whether that was done.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .geometry import subspace_distance, coord_distance, vector_angle

__all__ = ["NullDraws", "parametric_null", "deflation_surrogate", "rotate_toward_random"]


@dataclass
class NullDraws:
    """Per-pair null draws of both metrics under noise only."""

    sub: np.ndarray  # (n_pairs, n_draws)
    coord: np.ndarray  # (n_pairs, n_draws)
    meta: dict = field(default_factory=dict)


def _precision_weighted_mean(E_list, S_list) -> np.ndarray:
    """Pool contexts into one common truth, weighting by inverse variance."""
    W = np.stack([1.0 / np.maximum(S, 1e-12) ** 2 for S in S_list])
    X = np.stack([np.asarray(E, dtype=float) for E in E_list])
    return np.sum(W * X, axis=0) / np.sum(W, axis=0)


def parametric_null(
    E_by_context,
    S_by_context,
    n_draws: int,
    rng: np.random.Generator,
    metric: str = "chordal",
) -> NullDraws:
    """Null draws under "the contexts share one true E, only noise differs".

    For each pair, the shared truth is the precision-weighted mean across
    contexts, and each draw re-adds each context's own Gaussian noise with that
    context's own standard errors.  Feeding each context its *own* ``se`` keeps
    the null honest when the contexts have very different sample sizes: the
    null then carries the same asymmetry the observed comparison does.

    Parameters
    ----------
    E_by_context, S_by_context
        Two sequences of ``(n_pairs, r, 2)`` arrays: the effect estimates and
        their standard errors, context A first, context B second.
    """
    if len(E_by_context) != 2 or len(S_by_context) != 2:
        raise ValueError("parametric_null compares exactly two contexts")
    EA, EB = (np.asarray(x, dtype=float) for x in E_by_context)
    SA, SB = (np.asarray(x, dtype=float) for x in S_by_context)
    if not (EA.shape == EB.shape == SA.shape == SB.shape):
        raise ValueError("effect and standard-error arrays must share a shape")

    n_pairs = EA.shape[0]
    sub = np.empty((n_pairs, n_draws))
    coord = np.empty((n_pairs, n_draws))

    for p in range(n_pairs):
        truth = _precision_weighted_mean([EA[p], EB[p]], [SA[p], SB[p]])
        for d in range(n_draws):
            a = truth + rng.standard_normal(truth.shape) * SA[p]
            b = truth + rng.standard_normal(truth.shape) * SB[p]
            sub[p, d] = subspace_distance(a, b, metric=metric)
            coord[p, d] = coord_distance(a, b)

    return NullDraws(
        sub=sub,
        coord=coord,
        meta={"kind": "parametric", "n_draws": n_draws, "metric": metric, "matched_n": True},
    )


def rotate_toward_random(e, angle: float, target_norm: float, rng: np.random.Generator):
    """Rotate ``e`` by ``angle`` into a uniformly random orthogonal direction.

    Norm is set to ``target_norm``, so the surrogate matches the observed
    context-B column in both length and angular displacement and differs only in
    *where* it went.
    """
    e = np.asarray(e, dtype=float)
    n = np.linalg.norm(e)
    if n == 0.0:
        return e.copy()
    u = e / n
    w = rng.standard_normal(e.shape)
    w -= u * (u @ w)
    nw = np.linalg.norm(w)
    if nw == 0.0:
        return e.copy()
    w /= nw
    return target_norm * (np.cos(angle) * u + np.sin(angle) * w)


def deflation_surrogate(
    E_A,
    E_B,
    n_draws: int,
    rng: np.random.Generator,
    metric: str = "chordal",
) -> np.ndarray:
    """Subspace distances under "each SNP moved on its own, jointly by chance".

    Each column of ``E_A`` is rotated by *the angle that column actually moved*
    between the two contexts, but into an independently drawn random direction,
    and its length is set to the length that column actually has in context B.
    So the surrogate reproduces every single-SNP marginal -- per-SNP direction
    change and per-SNP magnitude -- and destroys only the joint structure.

    The observed pair sitting below this distribution is the one thing that
    makes the pair, and therefore the Grassmannian, worth talking about.  Two
    individually stable, non-collinear effect vectors span a stable plane for
    free; that is not joint conservation, it is arithmetic.
    """
    E_A = np.asarray(E_A, dtype=float)
    E_B = np.asarray(E_B, dtype=float)
    angles = [vector_angle(E_A[:, c], E_B[:, c]) for c in range(E_A.shape[1])]
    norms = [float(np.linalg.norm(E_B[:, c])) for c in range(E_B.shape[1])]

    out = np.empty(n_draws)
    for d in range(n_draws):
        cols = [
            rotate_toward_random(E_A[:, c], angles[c], norms[c], rng)
            for c in range(E_A.shape[1])
        ]
        out[d] = subspace_distance(E_A, np.column_stack(cols), metric=metric)
    return out
