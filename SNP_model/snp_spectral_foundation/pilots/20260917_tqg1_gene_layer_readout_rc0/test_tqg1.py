"""Asset-free tests for the TQ-G1 statistical and encoding contract.

    python -m unittest discover -s . -p "test_tqg1.py" -v

These run without the chr18 GEUVADIS assets, so the frozen estimator can be
checked before any real outcome is opened.
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

PILOT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(PILOT_DIR))
sys.path.insert(0, str(PILOT_DIR.parents[1] / "src"))

try:  # The statistical contract is testable without a torch install.
    import torch
    from snp_spectral_foundation.config import ModelConfig
    from snp_spectral_foundation.model import LocalGenotypeEncoder, encode_genotypes

    HAS_TORCH = True
except ModuleNotFoundError:  # pragma: no cover - exercised only on torch-free hosts
    HAS_TORCH = False
from tqg1_core import (  # noqa: E402
    af_missingness_dosage,
    decile_strata,
    deterministic_gene_order,
    fit_ridge,
    gene_pool_features,
    matched_null_spearman,
    paired_gene_bootstrap_difference,
    r2_against_train_mean,
    rankdata_average,
    spearman,
    tss_distance_weights,
    two_way_clustered_bootstrap,
    uniform_thin,
)

ALPHAS = [1e-6, 1e-4, 1e-2, 1.0, 100.0, 10000.0]


class TestSelection(unittest.TestCase):
    def test_uniform_thin_spans_and_counts(self) -> None:
        indices = np.arange(1000, 2000)
        picked = uniform_thin(indices, 64)
        self.assertEqual(len(picked), 64)
        self.assertTrue(np.all(np.diff(picked) > 0))
        self.assertEqual(picked[0], 1000)
        self.assertEqual(picked[-1], 1999)
        self.assertTrue(set(picked).issubset(set(indices.tolist())))

    def test_uniform_thin_refuses_upsampling(self) -> None:
        with self.assertRaises(ValueError):
            uniform_thin(np.arange(10), 20)

    def test_gene_order_is_deterministic_and_input_order_free(self) -> None:
        genes = [f"ENSG{i:08d}" for i in range(50)]
        first = deterministic_gene_order(genes, "tqg1-rc0")
        second = deterministic_gene_order(list(reversed(genes)), "tqg1-rc0")
        self.assertEqual(first, second)
        self.assertNotEqual(first, deterministic_gene_order(genes, "other-salt"))
        self.assertEqual(sorted(first), sorted(genes))


class TestGenotypeSummaries(unittest.TestCase):
    def test_af_missingness_and_imputation(self) -> None:
        genotype = np.array([[0, 2, -1], [2, 2, 1], [1, 2, 1], [1, 0, 1]], dtype=np.int8)
        train = np.array([0, 1, 2], dtype=np.int64)
        af, missing, dosage = af_missingness_dosage(genotype, train)
        self.assertAlmostEqual(af[0], 0.5)
        self.assertAlmostEqual(af[1], 1.0)
        self.assertAlmostEqual(af[2], 0.5)  # only two observed train carriers
        self.assertAlmostEqual(missing[2], 1.0 / 3.0)
        # The missing entry is filled with 2*AF, never with MAF.
        self.assertAlmostEqual(float(dosage[0, 2]), 2.0 * af[2], places=5)
        # Held-out rows are imputed with the train-derived value, not their own.
        self.assertAlmostEqual(float(dosage[3, 1]), 0.0)


class TestPooling(unittest.TestCase):
    def test_distance_kernel_is_normalised_and_monotone(self) -> None:
        positions = np.array([1000, 2000, 3000, 10000], dtype=np.int64)
        weights = tss_distance_weights(positions, 1000, 1000.0)
        self.assertAlmostEqual(float(weights.sum()), 1.0)
        self.assertTrue(np.all(np.diff(weights) < 0))

    def test_gene_pool_shape_and_individual_independence(self) -> None:
        rng = np.random.default_rng(0)
        hidden = rng.normal(size=(7, 4, 8, 5)).astype(np.float32)
        positions = np.sort(rng.integers(0, 1_000_000, size=32)).astype(np.int64)
        pooled = gene_pool_features(hidden, positions, 500_000, 100_000.0)
        self.assertEqual(pooled.shape, (7, 5 + 5 + 4 * 5))
        # Pooling one individual at a time must agree with the batched call.
        single = gene_pool_features(hidden[2:3], positions, 500_000, 100_000.0)
        np.testing.assert_allclose(single[0], pooled[2], rtol=1e-5, atol=1e-6)

    def test_gene_pool_rejects_mismatched_positions(self) -> None:
        hidden = np.zeros((2, 4, 8, 3), dtype=np.float32)
        with self.assertRaises(ValueError):
            gene_pool_features(hidden, np.arange(10, dtype=np.int64), 5, 1000.0)


class TestRidge(unittest.TestCase):
    def setUp(self) -> None:
        rng = np.random.default_rng(7)
        self.n = 90
        self.x = rng.normal(size=(self.n, 12))
        beta = np.array([2.0, -1.5] + [0.0] * 10)
        self.y = self.x @ beta + 0.1 * rng.normal(size=self.n)
        self.train = np.arange(0, 60, dtype=np.int64)
        self.val = np.arange(60, 75, dtype=np.int64)
        self.test = np.arange(75, 90, dtype=np.int64)

    def test_recovers_linear_signal(self) -> None:
        predict, info = fit_ridge(
            self.x, self.y[self.train], self.train, self.val, ALPHAS, self.y[self.val]
        )
        prediction = predict(self.x[self.test])
        r2 = r2_against_train_mean(
            self.y[self.test], prediction, float(self.y[self.train].mean())
        )
        self.assertGreater(r2, 0.9)
        self.assertIn(info["alpha"], ALPHAS)
        self.assertEqual(info["nonconstant_features"], 12)

    def test_constant_columns_are_dropped_not_fatal(self) -> None:
        x = np.column_stack([self.x, np.ones(self.n)])
        predict, info = fit_ridge(
            x, self.y[self.train], self.train, self.val, ALPHAS, self.y[self.val]
        )
        self.assertEqual(info["input_features"], 13)
        self.assertEqual(info["nonconstant_features"], 12)
        self.assertEqual(predict(x[self.test]).shape, (15,))

    def test_zeroing_a_feature_block_isolates_its_contribution(self) -> None:
        """The perturbation contrast must cancel every unperturbed column."""
        predict, _ = fit_ridge(
            self.x, self.y[self.train], self.train, self.val, ALPHAS, self.y[self.val]
        )
        perturbed = self.x.copy()
        perturbed[:, :2] = self.x[self.train][:, :2].mean(axis=0)
        delta = predict(self.x[self.test]) - predict(perturbed[self.test])
        # Columns 2..11 carry no signal, so the isolated delta must track y.
        self.assertGreater(abs(float(np.corrcoef(delta, self.y[self.test])[0, 1])), 0.9)

    def test_rejects_constant_outcome(self) -> None:
        with self.assertRaises(ValueError):
            fit_ridge(
                self.x, np.ones(len(self.train)), self.train, self.val, ALPHAS, self.y[self.val]
            )


class TestBootstrap(unittest.TestCase):
    def _panel(self, noise_b: float, noise_c: float, seed: int = 3):
        rng = np.random.default_rng(seed)
        genes, people = 40, 30
        truth = rng.normal(size=(genes, people))
        train_means = np.zeros(genes)
        predictions = {
            "B": truth * 0.5 + noise_b * rng.normal(size=(genes, people)),
            "C_gene": truth * 0.5 + noise_c * rng.normal(size=(genes, people)),
        }
        return truth, predictions, train_means

    def test_better_arm_passes(self) -> None:
        truth, predictions, means = self._panel(noise_b=0.9, noise_c=0.2)
        out = two_way_clustered_bootstrap(
            truth, predictions, means, ("C_gene", "B"), 11, 400
        )
        self.assertGreater(out["macro_delta_r2"], 0)
        self.assertTrue(out["pass"])
        self.assertEqual(out["resampling_units"], ["evaluation_individual", "gene"])

    def test_identical_arms_give_zero_delta_and_a_covering_interval(self) -> None:
        truth, predictions, means = self._panel(noise_b=0.5, noise_c=0.5)
        predictions["C_gene"] = predictions["B"].copy()
        out = two_way_clustered_bootstrap(
            truth, predictions, means, ("C_gene", "B"), 12, 400
        )
        self.assertAlmostEqual(out["macro_delta_r2"], 0.0, places=12)
        self.assertFalse(out["pass"])
        low, high = out["paired_bootstrap_ci95"]
        self.assertLessEqual(low, 0.0)
        self.assertGreaterEqual(high, 0.0)

    def test_gene_resampling_widens_the_interval(self) -> None:
        """Dropping the gene cluster would understate uncertainty."""
        truth, predictions, means = self._panel(noise_b=0.9, noise_c=0.85)
        two_way = two_way_clustered_bootstrap(
            truth, predictions, means, ("C_gene", "B"), 13, 800
        )
        low, high = two_way["paired_bootstrap_ci95"]
        one_gene = two_way_clustered_bootstrap(
            truth[:1], {k: v[:1] for k, v in predictions.items()}, means[:1], ("C_gene", "B"), 13, 800
        )
        self.assertGreater(high - low, 0.0)
        self.assertEqual(one_gene["n_genes"], 1)

    def test_shape_mismatch_is_rejected(self) -> None:
        truth, predictions, means = self._panel(0.5, 0.5)
        predictions["C_gene"] = predictions["C_gene"][:, :5]
        with self.assertRaises(ValueError):
            two_way_clustered_bootstrap(truth, predictions, means, ("C_gene", "B"), 1, 10)


class TestRankStatistics(unittest.TestCase):
    def test_average_ties(self) -> None:
        np.testing.assert_allclose(
            rankdata_average(np.array([10.0, 20.0, 20.0, 30.0])), [1.0, 2.5, 2.5, 4.0]
        )

    def test_spearman_endpoints(self) -> None:
        x = np.arange(20.0)
        self.assertAlmostEqual(spearman(x, x), 1.0, places=10)
        self.assertAlmostEqual(spearman(x, -x), -1.0, places=10)
        self.assertAlmostEqual(spearman(x, np.zeros(20)), 0.0, places=10)

    def test_decile_strata_are_bounded(self) -> None:
        strata = decile_strata(np.random.default_rng(1).normal(size=200))
        self.assertGreaterEqual(strata.min(), 0)
        self.assertLessEqual(strata.max(), 9)


class TestMatchedNull(unittest.TestCase):
    def test_null_score_is_not_significant(self) -> None:
        rng = np.random.default_rng(5)
        target = rng.normal(size=120)
        score = rng.normal(size=120)
        strata = decile_strata(rng.normal(size=120))
        out = matched_null_spearman(score, target, strata, 21, 500)
        self.assertGreater(out["p_two_sided"], 0.05)

    def test_strong_signal_is_significant(self) -> None:
        rng = np.random.default_rng(6)
        target = rng.normal(size=120)
        score = target + 0.05 * rng.normal(size=120)
        strata = np.zeros(120, dtype=np.int64)
        out = matched_null_spearman(score, target, strata, 22, 500)
        self.assertGreater(out["observed_spearman"], 0.9)
        self.assertLess(out["p_one_sided_greater"], 0.01)

    def test_paired_difference_detects_a_better_score(self) -> None:
        rng = np.random.default_rng(8)
        target = rng.normal(size=150)
        good = target + 0.2 * rng.normal(size=150)
        poor = target + 2.0 * rng.normal(size=150)
        out = paired_gene_bootstrap_difference(good, poor, target, 31, 400)
        self.assertTrue(out["pass"])
        reversed_out = paired_gene_bootstrap_difference(poor, good, target, 31, 400)
        self.assertFalse(reversed_out["pass"])


@unittest.skipUnless(HAS_TORCH, "torch is required for the encoder contract tests")
class TestEncoderPerturbation(unittest.TestCase):
    """The gene perturbation must remove exactly the cis genotype information."""

    def setUp(self) -> None:
        torch.manual_seed(0)
        self.cfg = ModelConfig(
            max_blocks=4,
            snps_per_block=8,
            d_model=16,
            n_heads=4,
            local_layers=1,
            local_ff_dim=32,
            spectral_rank=4,
            dropout=0.0,
            local_position_mode="relative_continuous",
            use_local_block_embedding=False,
        )
        self.model = LocalGenotypeEncoder(self.cfg).eval()
        rng = np.random.default_rng(2)
        self.genotype = torch.from_numpy(
            rng.integers(0, 3, size=(6, 4, 8)).astype(np.int64)
        )
        self.af = torch.from_numpy(rng.uniform(0.05, 0.95, size=(4, 8)).astype(np.float32))
        positions = np.sort(rng.integers(0, 2_000_000, size=32)).astype(np.int64)
        self.positions = torch.from_numpy(positions.reshape(4, 8))

    def test_full_mask_erases_individual_differences(self) -> None:
        mask = torch.ones_like(self.genotype, dtype=torch.bool)
        with torch.no_grad():
            hidden, _ = self.model(
                encode_genotypes(self.genotype, mask),
                self.af,
                genomic_positions=self.positions,
            )
        spread = hidden.std(dim=0).max().item()
        self.assertLess(spread, 1e-6, "masked encoding must not depend on the individual")

    def test_unmasked_encoding_does_depend_on_the_individual(self) -> None:
        with torch.no_grad():
            hidden, _ = self.model(
                encode_genotypes(self.genotype), self.af, genomic_positions=self.positions
            )
        self.assertGreater(hidden.std(dim=0).max().item(), 1e-3)

    def test_encoder_has_no_locus_specific_parameters(self) -> None:
        """A shared encoder must carry no block lookup table."""
        self.assertIsNone(self.model.block_embedding)
        self.assertIsNone(self.model.position_embedding)
        names = [n for n, _ in self.model.named_parameters()]
        self.assertFalse([n for n in names if "block_embedding" in n])

    def test_block_ids_are_prohibited(self) -> None:
        with self.assertRaises(ValueError):
            self.model(
                encode_genotypes(self.genotype),
                self.af,
                torch.zeros(4, dtype=torch.long),
                genomic_positions=self.positions,
            )

    def test_float_positions_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            self.model(
                encode_genotypes(self.genotype),
                self.af,
                genomic_positions=self.positions.float(),
            )


if __name__ == "__main__":
    unittest.main(verbosity=2)
