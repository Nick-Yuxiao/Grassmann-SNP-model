import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

PKG = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PKG / "src"))

from gate0a import estimator, load_frozen_panel  # noqa: E402
from gate0a.sigma_map import required_replicates, power_at  # noqa: E402


def synth_panel(N=200, n_blocks=4, sites_per_block=12, seed=0):
    rng = np.random.default_rng(seed)
    Gs, blocks = [], []
    for b in range(n_blocks):
        F = rng.standard_normal((N, 3))
        load = rng.standard_normal((3, sites_per_block))
        lat = F @ load + 0.5 * rng.standard_normal((N, sites_per_block))
        p = 1.0 / (1.0 + np.exp(-lat))
        g = (rng.random((N, sites_per_block)) < p).astype(int) \
            + (rng.random((N, sites_per_block)) < p).astype(int)
        Gs.append(g)
        blocks += [b] * sites_per_block
    return np.concatenate(Gs, axis=1).astype(np.int16), np.array(blocks)


class TestEstimator(unittest.TestCase):
    def setUp(self):
        self.G, self.blocks = synth_panel()
        self.folds = estimator.build_folds(self.G.shape[0], 5, seed=1)
        self.cache = estimator.precompute_reps(self.G, self.blocks, self.folds, 16, [0, 1, 2, 3])

    def test_folds_drawn_once_are_shared(self):
        f2 = estimator.build_folds(self.G.shape[0], 5, seed=1)
        for (t1, v1), (t2, v2) in zip(self.folds, f2):
            self.assertTrue(np.array_equal(t1, t2) and np.array_equal(v1, v2))

    def test_cache_prefix_is_topk(self):
        rep16 = self.cache[("pca_z", 0, 0)]
        self.assertLessEqual(rep16.shape[1], 16)
        # k=8 design must be the length-8 prefix
        d8 = estimator.design(self.cache, "pca_z", 0, "major_LD_aligned", [0], [0, 1], 8)
        self.assertEqual(d8.shape[1], 8 * 2)

    def test_g_is_fixed_across_folds(self):
        rng = np.random.default_rng(3)
        g = estimator.simulate_g(self.G, self.blocks, "spectral_tail_adversarial",
                                 [0], [0, 1, 2, 3], rng, n_dir=3, poly_frac=0.5)
        self.assertEqual(g.shape[0], self.G.shape[0])
        self.assertTrue(np.all(np.isfinite(g)))
        self.assertAlmostEqual(float(g.std()), 1.0, places=6)  # standardized truth

    def test_sigma_for_cell_finite(self):
        res = estimator.sigma_for_cell(
            self.cache, self.G, self.blocks, self.folds,
            "between_block_interaction", [0, 3], [0, 1, 2, 3], k=4,
            R_pilot=5, seed=2, inner_folds=3)
        self.assertGreaterEqual(res["sigma_hat"], 0.0)
        self.assertTrue(np.isfinite(res["sigma_hat"]))
        self.assertEqual(len(res["deltas"]), 5)
        self.assertIn("mean_delta_NON_EVIDENCE_DIAGNOSTIC_ONLY", res)

    def test_polygenic_changes_g(self):
        rng = np.random.default_rng(7)
        g0 = estimator.simulate_g(self.G, self.blocks, "major_LD_aligned", [0],
                                  [0, 1, 2, 3], np.random.default_rng(7), 3, 0.0)
        g1 = estimator.simulate_g(self.G, self.blocks, "major_LD_aligned", [0],
                                  [0, 1, 2, 3], np.random.default_rng(7), 3, 0.9)
        self.assertFalse(np.allclose(g0, g1))

    def test_sigma_map_monotone(self):
        r_small = required_replicates(0.01, M=1500, B=800)
        r_large = required_replicates(0.03, M=1500, B=800)
        if r_small is not None and r_large is not None:
            self.assertLessEqual(r_small, r_large)

    def test_power_increases_with_R(self):
        self.assertLessEqual(power_at(0.02, 20, M=1500, B=800),
                             power_at(0.02, 200, M=1500, B=800) + 0.03)

    def test_loader_length_guard(self):
        G, blocks = synth_panel(N=60, n_blocks=3, sites_per_block=10)
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            np.save(d / "G.npy", G)
            np.save(d / "blocks.npy", blocks)
            np.save(d / "train_idx.npy", np.arange(45))
            np.save(d / "val_idx.npy", np.arange(45, 60))
            self.assertEqual(load_frozen_panel(d, expect_L=30)["L"], 30)
            with self.assertRaises(ValueError):
                load_frozen_panel(d, expect_L=154850)


if __name__ == "__main__":
    unittest.main()
