"""Invariance, numerical-stability and identity tests for the geometry layer."""

import unittest
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from subspace_conservation.geometry import (  # noqa: E402
    coord_distance,
    gl2_alignment_residual,
    orthonormalize,
    principal_angles,
    projector,
    random_gl2,
    random_orthogonal,
    subspace_distance,
    vector_angle,
)


class TestInvariance(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(20260921)

    def test_subspace_distance_is_gl2_invariant(self):
        for cond in (1.0, 2.0, 10.0, 1e3):
            E = self.rng.standard_normal((9, 2))
            M = random_gl2(cond, self.rng)
            self.assertLess(subspace_distance(E, E @ M), 1e-10, f"cond={cond}")

    def test_projector_is_gl2_invariant(self):
        E = self.rng.standard_normal((7, 2))
        M = random_gl2(25.0, self.rng)
        np.testing.assert_allclose(projector(E), projector(E @ M), atol=1e-10)

    def test_second_moment_is_o2_invariant_but_not_gl2(self):
        # The correction that reshapes the whole design: under right-O(2) the
        # second moment survives untouched and is strictly richer than the
        # projector, so an O(2)-only world gives the spectrum arm a free win.
        E = self.rng.standard_normal((6, 2))
        R = random_orthogonal(2, self.rng)
        np.testing.assert_allclose(E @ R @ (E @ R).T, E @ E.T, atol=1e-10)
        M = random_gl2(8.0, self.rng)
        self.assertGreater(np.linalg.norm(E @ M @ (E @ M).T - E @ E.T), 1e-3)

    def test_coord_distance_is_not_invariant(self):
        E = self.rng.standard_normal((6, 2))
        self.assertGreater(coord_distance(E, E @ random_gl2(4.0, self.rng)), 1e-2)

    def test_column_sign_flip_moves_coordinates_not_subspace(self):
        # A flipped allele is a sign flip, a sign flip is in GL(2). This is the
        # bookkeeping error that would read exactly like the hypothesis.
        E = self.rng.standard_normal((8, 2))
        F = E * np.array([1.0, -1.0])
        self.assertLess(subspace_distance(E, F), 1e-10)
        self.assertGreater(coord_distance(E, F), 0.1)


class TestNumerics(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(7)

    def test_small_angles_survive_the_sine_branch(self):
        # arccos of a cosine this close to 1 returns garbage; the sine branch
        # must recover the angle to several digits.
        for eps in (1e-4, 1e-6, 1e-8):
            Q, _ = np.linalg.qr(self.rng.standard_normal((10, 2)))
            perturb = self.rng.standard_normal((10, 2))
            perturb -= Q @ (Q.T @ perturb)
            perturb /= np.linalg.norm(perturb, axis=0)
            theta = principal_angles(Q, Q + eps * perturb)
            self.assertTrue(np.all(theta < 10 * eps), f"eps={eps}, theta={theta}")
            self.assertTrue(np.all(theta > 0.1 * eps), f"eps={eps}, theta={theta}")

    def test_chordal_equals_projector_frobenius_norm(self):
        A = self.rng.standard_normal((9, 2))
        B = self.rng.standard_normal((9, 2))
        d = subspace_distance(A, B, metric="chordal")
        expect = np.linalg.norm(projector(A) - projector(B)) / np.sqrt(2 * 2)
        self.assertAlmostEqual(d, expect, places=10)

    def test_distances_are_bounded(self):
        for _ in range(50):
            A = self.rng.standard_normal((5, 2))
            B = self.rng.standard_normal((5, 2))
            for metric in ("chordal", "geodesic"):
                d = subspace_distance(A, B, metric=metric)
                self.assertGreaterEqual(d, -1e-12)
                self.assertLessEqual(d, 1.0 + 1e-12)
            self.assertLessEqual(coord_distance(A, B), 1.0 + 1e-12)

    def test_orthonormalize_drops_degenerate_directions(self):
        v = self.rng.standard_normal((6, 1))
        E = np.hstack([v, v * 2.0])  # rank 1
        Q, _ = orthonormalize(E)
        self.assertEqual(Q.shape[1], 1)

    def test_vector_angle_stable_at_both_ends(self):
        u = self.rng.standard_normal(8)
        self.assertLess(vector_angle(u, u * 3.0), 1e-12)
        self.assertAlmostEqual(vector_angle(u, -u), np.pi, places=10)


class TestGL2AlignmentIdentity(unittest.TestCase):
    def test_per_pair_gl2_residual_is_the_chordal_distance(self):
        # Grassmann == aligned with a per-pair GL(2) adaptation already applied.
        # So a baseline that fits M on the test pair's own target data is not a
        # baseline; it is this arm plus leakage.
        rng = np.random.default_rng(11)
        for _ in range(20):
            A = rng.standard_normal((10, 2))
            B = rng.standard_normal((10, 2))
            self.assertAlmostEqual(
                gl2_alignment_residual(A, B),
                subspace_distance(A, B, metric="chordal"),
                places=10,
            )


if __name__ == "__main__":
    unittest.main()
