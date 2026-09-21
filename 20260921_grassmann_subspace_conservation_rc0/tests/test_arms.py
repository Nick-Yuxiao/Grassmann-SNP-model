"""Feature maps carry the invariance they claim, and the ridge head is fair."""

import unittest
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from subspace_conservation.arms import (  # noqa: E402
    ARMS,
    INVARIANCE,
    features,
    fit_gl_transport,
    fit_procrustes,
    ridge_fit_predict,
    transfer_r2,
    vech,
)
from subspace_conservation.geometry import random_gl2, random_orthogonal  # noqa: E402
from subspace_conservation.geometry import subspace_distance  # noqa: E402
from subspace_conservation.simulate import (  # noqa: E402
    DISTORTION_LEVELS,
    apply_distortion,
    apply_trait_distortion,
)


def _distort(E, M):
    return np.einsum("nrk,kl->nrl", E, M)


class TestArmInvariance(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(20260921)
        self.E = self.rng.standard_normal((30, 8, 2))

    def test_declared_invariance_matches_behaviour(self):
        R = random_orthogonal(2, self.rng)
        M = random_gl2(8.0, self.rng)
        for arm, claim in INVARIANCE.items():
            if arm in ("procrustes", "gl_transport", "hybrid"):
                continue
            f = features(arm, self.E, np.random.default_rng(0))
            f_o2 = features(arm, _distort(self.E, R), np.random.default_rng(0))
            f_gl = features(arm, _distort(self.E, M), np.random.default_rng(0))
            o2_inv = np.allclose(f, f_o2, atol=1e-8)
            gl_inv = np.allclose(f, f_gl, atol=1e-8)
            if claim == "GL(2)":
                self.assertTrue(o2_inv and gl_inv, arm)
            elif claim == "O(2)":
                self.assertTrue(o2_inv, arm)
                self.assertFalse(gl_inv, f"{arm} must lose GL(2)")
            else:
                self.assertFalse(gl_inv, arm)

    def test_bilinear_control_matches_grassmann_capacity(self):
        g = features("grassmann", self.E, np.random.default_rng(0))
        for arm in ("bilinear_o2", "bilinear_raw"):
            self.assertEqual(
                features(arm, self.E, np.random.default_rng(0)).shape[1], g.shape[1], arm
            )

    def test_every_declared_arm_is_buildable(self):
        T = np.eye(8)
        for arm in ARMS:
            f = features(arm, self.E, np.random.default_rng(0), transport=T)
            self.assertEqual(f.shape[0], self.E.shape[0], arm)
            self.assertTrue(np.all(np.isfinite(f)), arm)

    def test_vech_preserves_the_frobenius_norm(self):
        S = self.rng.standard_normal((5, 6, 6))
        S = S + np.transpose(S, (0, 2, 1))
        np.testing.assert_allclose(
            np.linalg.norm(vech(S), axis=1), np.linalg.norm(S, axis=(1, 2)), rtol=1e-10
        )


class TestTransports(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(5)

    def test_procrustes_recovers_a_shared_left_rotation(self):
        E = self.rng.standard_normal((300, 6, 2))
        Q_true = random_orthogonal(6, self.rng)
        E_ctx = np.einsum("rs,nsk->nrk", Q_true, E)
        Q_hat = fit_procrustes(E, E_ctx)  # map context back to source frame
        np.testing.assert_allclose(Q_hat, Q_true.T, atol=1e-8)

    def test_gl_transport_recovers_a_shared_left_map(self):
        E = self.rng.standard_normal((400, 6, 2))
        W_true = self.rng.standard_normal((6, 6))
        E_ctx = np.einsum("rs,nsk->nrk", W_true, E)
        W_hat = fit_gl_transport(E, E_ctx, alpha=1e-8)
        np.testing.assert_allclose(W_hat, np.linalg.inv(W_true), atol=1e-4)

    def test_pair_specific_distortion_defeats_any_shared_transport(self):
        # The point of the whole design: under the LD mechanism M is the pair's
        # own Lambda, so no map shared across pairs can undo it.
        E = self.rng.standard_normal((400, 6, 2))
        E_ctx = apply_distortion(E, "gl_8", self.rng, shared=False)
        W = fit_gl_transport(E, E_ctx, alpha=1e-3)
        resid_raw = np.linalg.norm(E - E_ctx)
        resid_fit = np.linalg.norm(E - np.einsum("rs,nsk->nrk", W, E_ctx))
        self.assertGreater(resid_fit, 0.5 * resid_raw)

    def test_shared_right_distortion_also_defeats_a_left_transport(self):
        # Sharing M across pairs does not make a right-acting nuisance
        # left-repairable: the two sides never meet. Only a right-side 2x2
        # transport could undo this one, and the subspace is immune to it
        # anyway, so it is not the control that keeps the comparison honest.
        E = self.rng.standard_normal((400, 6, 2))
        E_ctx = apply_distortion(E, "gl_4", self.rng, shared=True)
        W = fit_gl_transport(E, E_ctx, alpha=1e-3)
        resid_raw = np.linalg.norm(E - E_ctx)
        resid_fit = np.linalg.norm(E - np.einsum("rs,nsk->nrk", W, E_ctx))
        self.assertGreater(resid_fit, 0.3 * resid_raw)

    def test_left_acting_trait_nuisance_is_repairable_and_moves_the_subspace(self):
        # The anti-rigging control, both halves. If only the first held, the
        # transport arm would be facing a nuisance it cannot ever fix; if only
        # the second, the Grassmann arm would never be exposed to a nuisance it
        # cannot absorb. The sweep needs a regime where the ordering reverses.
        E = self.rng.standard_normal((400, 6, 2))
        E_ctx, _ = apply_trait_distortion(E, 0.6, self.rng)
        W = fit_gl_transport(E, E_ctx, alpha=1e-8)
        resid_raw = np.linalg.norm(E - E_ctx)
        resid_fit = np.linalg.norm(E - np.einsum("rs,nsk->nrk", W, E_ctx))
        self.assertLess(resid_fit, 0.05 * resid_raw)

        moved = np.mean([subspace_distance(E[p], E_ctx[p]) for p in range(100)])
        self.assertGreater(moved, 0.05)


class TestRidgeHead(unittest.TestCase):
    def test_recovers_a_linear_signal(self):
        rng = np.random.default_rng(3)
        X = rng.standard_normal((500, 10))
        w = rng.standard_normal(10)
        y = X @ w + 0.1 * rng.standard_normal(500)
        Xe = rng.standard_normal((200, 10))
        ye = Xe @ w
        pred, alpha = ridge_fit_predict(X, y, Xe)
        self.assertGreater(transfer_r2(ye, pred, y), 0.95)
        self.assertIn(alpha, (1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0, 1e3))

    def test_r2_is_negative_when_prediction_is_worse_than_the_train_mean(self):
        y_train = np.zeros(50)
        self.assertLess(transfer_r2(np.ones(50), np.full(50, 5.0), y_train), 0.0)

    def test_distortion_levels_are_ordered_and_complete(self):
        self.assertEqual(DISTORTION_LEVELS[0], "none")
        self.assertEqual(DISTORTION_LEVELS[1], "orth")
        self.assertTrue(all(lv.startswith("gl_") for lv in DISTORTION_LEVELS[2:]))


if __name__ == "__main__":
    unittest.main()
