"""Subspace geometry for pair effect matrices.

A pair (i, j) measured over r contexts-of-measurement (traits) gives
``E = [e_i, e_j] in R^{r x 2}``.  The object this package is about is the
column space ``span(E)``, a point of the Grassmannian ``Gr(2, r)``.

The nuisance group is ``GL(2)`` acting on the right, not ``O(2)``:

    E = B @ Lambda.T            B: r x K causal effect vectors
                                Lambda: 2 x K LD loadings of the typed SNPs

Changing ancestry changes ``Lambda``.  With ``K = 2`` and ``Lambda`` invertible
this is exactly ``E -> E M`` with ``M = (Lambda^{-1} Lambda').T in GL(2)``, and
``span(E M) = span(E)``.  The distinction matters: the second moment
``E E.T`` is *exactly* invariant under right-``O(2)`` (``E R R.T E.T = E E.T``)
and is strictly richer than the projector, so a design that only rotates gives
a spectrum arm a free win over the Grassmann arm.  Only under ``GL(2)`` is the
column space the maximal invariant.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "orthonormalize",
    "principal_angles",
    "subspace_distance",
    "projector",
    "coord_distance",
    "vector_angle",
    "random_gl2",
    "random_orthogonal",
    "gl2_alignment_residual",
]

# Below this cosine the arccos branch is well conditioned; above it we must use
# the arcsin branch or small angles are destroyed (Bjorck-Golub / Knyazev).
_COS_BRANCH = 1.0 / np.sqrt(2.0)


def orthonormalize(E, rcond: float = 1e-10):
    """Return ``(Q, s)``: an orthonormal basis of ``span(E)`` and its spectrum.

    ``Q`` has one column per numerically non-degenerate singular direction, so a
    rank-deficient ``E`` yields fewer columns than ``E`` has.  Callers that need
    a fixed ``k`` must check ``Q.shape[1]`` themselves rather than assume it.
    """
    E = np.asarray(E, dtype=float)
    U, s, _ = np.linalg.svd(E, full_matrices=False)
    if s.size == 0 or s[0] == 0.0:
        return U[:, :0], s
    keep = s > rcond * s[0]
    return U[:, keep], s


def principal_angles(E_A, E_B, rcond: float = 1e-10) -> np.ndarray:
    """Principal angles between ``span(E_A)`` and ``span(E_B)``, ascending.

    Uses the sine branch for small angles: ``cos`` loses all precision near 1,
    which is exactly the regime this package spends its time in (two contexts
    that agree).  Returns ``min(rank_A, rank_B)`` angles in ``[0, pi/2]``.
    """
    QA, _ = orthonormalize(E_A, rcond)
    QB, _ = orthonormalize(E_B, rcond)
    if QA.shape[1] == 0 or QB.shape[1] == 0:
        return np.zeros(0)

    Y, cos_theta, Zt = np.linalg.svd(QA.T @ QB, full_matrices=False)
    cos_theta = np.clip(cos_theta, 0.0, 1.0)

    # Principal vectors on the B side; their components orthogonal to span(A)
    # have norm sin(theta_i) exactly, which stays accurate as theta -> 0.
    V = QB @ Zt.T
    residual = V - QA @ (QA.T @ V)
    sin_theta = np.clip(np.linalg.norm(residual, axis=0), 0.0, 1.0)

    theta = np.where(cos_theta > _COS_BRANCH, np.arcsin(sin_theta), np.arccos(cos_theta))
    return np.sort(theta)


def subspace_distance(E_A, E_B, metric: str = "chordal", rcond: float = 1e-10) -> float:
    """Distance between the two column spaces, normalised to ``[0, 1]``.

    ``chordal``  -- ``sqrt(sum sin^2 theta / k)``; equals
                    ``||P_A - P_B||_F / sqrt(2k)``.
    ``geodesic`` -- ``sqrt(sum theta^2) / (sqrt(k) * pi/2)``.

    Both are normalised so that the two metrics, and subspaces of different
    dimension, land on a common ``[0, 1]`` scale.  That normalisation is a
    convenience only: the comparison against the coordinate metric is made
    through the null calibration in :mod:`estimators`, never by reading these
    numbers against each other directly.
    """
    theta = principal_angles(E_A, E_B, rcond)
    k = theta.size
    if k == 0:
        return float("nan")
    if metric == "chordal":
        return float(np.sqrt(np.sum(np.sin(theta) ** 2) / k))
    if metric == "geodesic":
        return float(np.sqrt(np.sum(theta**2)) / (np.sqrt(k) * (np.pi / 2)))
    raise ValueError(f"unknown metric: {metric!r}")


def projector(E, rcond: float = 1e-10) -> np.ndarray:
    """Orthogonal projector onto ``span(E)``; invariant under ``E -> E M``."""
    Q, _ = orthonormalize(E, rcond)
    return Q @ Q.T


def coord_distance(E_A, E_B) -> float:
    """Aligned-coordinate distance, normalised to ``[0, 1]``.

    Each matrix is divided by its Frobenius norm first, so a pure change of
    overall effect scale between contexts -- different phenotype units, a
    different amount of attenuation -- is not charged to the coordinate arm.
    Anything beyond that, including the column mixing that ``span`` absorbs, is.

    PRECONDITION: the two columns are the same two SNPs in the same order, with
    the same allele coding, in both contexts.  A flipped allele is a sign flip,
    a sign flip lies in ``GL(2)``, so it leaves the subspace untouched while it
    inflates this number -- that would manufacture a Grassmann advantage out of
    a bookkeeping error.  ``estimators.orientation_diagnostic`` checks it.
    """
    A = np.asarray(E_A, dtype=float)
    B = np.asarray(E_B, dtype=float)
    na, nb = np.linalg.norm(A), np.linalg.norm(B)
    if na == 0.0 or nb == 0.0:
        return float("nan")
    return float(np.linalg.norm(A / na - B / nb) / 2.0)


def vector_angle(u, v) -> float:
    """Angle in ``[0, pi]`` between two vectors, stable near 0 and pi."""
    u = np.asarray(u, dtype=float)
    v = np.asarray(v, dtype=float)
    nu, nv = np.linalg.norm(u), np.linalg.norm(v)
    if nu == 0.0 or nv == 0.0:
        return float("nan")
    uh, vh = u / nu, v / nv
    # atan2 of the orthogonal and parallel parts: no cancellation at either end.
    return float(np.arctan2(np.linalg.norm(uh - vh * (vh @ uh)), uh @ vh))


def random_orthogonal(k: int, rng: np.random.Generator) -> np.ndarray:
    """Haar-distributed ``O(k)`` matrix."""
    A = rng.standard_normal((k, k))
    Q, R = np.linalg.qr(A)
    return Q * np.sign(np.diag(R))


def random_gl2(cond: float, rng: np.random.Generator) -> np.ndarray:
    """Random ``GL(2)`` matrix with exactly the requested condition number.

    ``M = U diag(1, 1/cond) V.T`` with Haar ``U, V``.  ``cond = 1`` degenerates
    to a pure rotation/reflection, which is the ``O(2)`` sub-case -- keep that
    level in any sweep, it is what separates "the spectrum arm is also
    invariant" from "only the subspace survives".
    """
    if cond < 1.0:
        raise ValueError("cond must be >= 1")
    U = random_orthogonal(2, rng)
    V = random_orthogonal(2, rng)
    return U @ np.diag([1.0, 1.0 / cond]) @ V.T


def gl2_alignment_residual(E_A, E_B) -> float:
    """Residual of the best per-pair ``GL(2)`` alignment of ``E_A`` onto ``E_B``.

    ``min_M ||Bq - Aq M||_F / sqrt(k)`` on orthonormalised inputs.  This is not
    a separate method: least squares puts the residual at
    ``||(I - P_A) Q_B||_F``, which *is* the chordal subspace distance.  So the
    Grassmann arm is exactly the aligned arm with a per-pair ``GL(2)``
    adaptation already applied, at zero adaptation cost -- and any baseline that
    fits ``M`` on the test pair's own target data is not a baseline, it is the
    Grassmann arm plus leakage.  :mod:`tests` asserts the identity.
    """
    QA, _ = orthonormalize(E_A)
    QB, _ = orthonormalize(E_B)
    if QA.shape[1] == 0 or QB.shape[1] == 0:
        return float("nan")
    M, *_ = np.linalg.lstsq(QA, QB, rcond=None)
    return float(np.linalg.norm(QB - QA @ M) / np.sqrt(QB.shape[1]))
