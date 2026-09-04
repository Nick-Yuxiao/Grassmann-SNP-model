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
        # LD via a few shared latent factors per block
        F = rng.standard_normal((N, 3))
        load = rng.standard_normal((3, sites_per_block))
        lat = F @ load + 0.5 * rng.standard_normal((N, sites_per_block))
        p = 1.0 / (1.0 + np.exp(-lat))          # per-cell prob
        g = (rng.random((N, sites_per_block)) < p).astype(int) \
            + (rng.random((N, sites_per_block)) < p).astype(int)  # 0/1/2
        Gs.append(g)
        blocks += [b] * sites_per_block
    return np.concatenate(Gs, axis=1).astype(np.int16), np.array(blocks)


class TestEstimator(unittest.TestCase):
    def test_arms_shapes(self):
        G, blocks = synth_panel()
        tr = np.arange(0, 150)
        for arm in ("pca_z", "kpca_rbf"):
            reps = estimator.ARMS[arm](G, blocks, k=4, tr=tr)
            self.assertEqual(len(reps), len(np.unique(blocks)))
            for b, r in reps.items():
                self.assertEqual(r.shape[0], G.shape[0])
                self.assertLessEqual(r.shape[1], 4)
                self.assertTrue(np.all(np.isfinite(r)))

    def test_replicate_delta_finite(self):
        G, blocks = synth_panel()
        rng = np.random.default_rng(1)
        d = estimator.replicate_delta(G, blocks, "spectral_tail_adversarial",
                                      [0], k=4, rng=rng, folds=3)
        self.assertTrue(np.isfinite(d))

    def test_sigma_for_cell(self):
        G, blocks = synth_panel()
        res = estimator.sigma_for_cell(G, blocks, "between_block_interaction",
                                       [0, 3], k=4, R_pilot=5, seed=2, folds=3)
        self.assertGreaterEqual(res["sigma_hat"], 0.0)
        self.assertTrue(np.isfinite(res["sigma_hat"]))
        self.assertEqual(len(res["deltas"]), 5)

    def test_sigma_map_monotone(self):
        # smaller sigma needs no more replicates than larger sigma
        r_small = required_replicates(0.01, M=1500, B=800)
        r_large = required_replicates(0.03, M=1500, B=800)
        if r_small is not None and r_large is not None:
            self.assertLessEqual(r_small, r_large)

    def test_power_increases_with_R(self):
        p20 = power_at(0.02, 20, M=1500, B=800)
        p200 = power_at(0.02, 200, M=1500, B=800)
        self.assertLessEqual(p20, p200 + 0.03)

    def test_loader_roundtrip_and_length_guard(self):
        G, blocks = synth_panel(N=60, n_blocks=3, sites_per_block=10)
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            np.save(d / "G.npy", G)
            np.save(d / "blocks.npy", blocks)
            np.save(d / "train_idx.npy", np.arange(45))
            np.save(d / "val_idx.npy", np.arange(45, 60))
            panel = load_frozen_panel(d, expect_L=30)   # 3*10
            self.assertEqual(panel["L"], 30)
            with self.assertRaises(ValueError):
                load_frozen_panel(d, expect_L=154850)   # contract mismatch


if __name__ == "__main__":
    unittest.main()
