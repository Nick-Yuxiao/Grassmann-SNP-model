"""Asset-free tests for TQ-B1. numpy only, no cohort, no GPU.

    python -m unittest test_tqb1 -v
"""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from smoke_bar import write_bed  # noqa: E402
from tqb1_core import (  # noqa: E402
    PrimalRidge,
    af_missingness_dosage,
    clustered_macro_bootstrap_shared,
    decile_strata,
    deterministic_order,
    gene_bootstrap_spearman,
    matched_null_spearman,
    perturbation_score,
    r2_against_train_mean,
    rankdata_average,
    read_bed_variants,
    spearman,
    uniform_thin,
)
from tqb1_io import assign_splits  # noqa: E402

ALPHAS = [1e-4, 1e-2, 1.0, 100.0, 1e4, 1e6]


class TestBedReader(unittest.TestCase):
    def test_round_trip_for_awkward_sample_counts(self) -> None:
        rng = np.random.default_rng(1)
        for n in (1, 2, 3, 4, 5, 7, 8, 9, 13, 256):
            genotype = rng.integers(-1, 3, size=(n, 11)).astype(np.int8)
            path = Path(tempfile.mkdtemp()) / "t.bed"
            write_bed(path, genotype)
            back = read_bed_variants(path, n, np.arange(11), 11)
            np.testing.assert_array_equal(genotype, back, err_msg=f"n={n}")

    def test_random_access_matches_full_read(self) -> None:
        rng = np.random.default_rng(2)
        genotype = rng.integers(0, 3, size=(37, 25)).astype(np.int8)
        path = Path(tempfile.mkdtemp()) / "t.bed"
        write_bed(path, genotype)
        picked = np.array([24, 3, 17], dtype=np.int64)
        np.testing.assert_array_equal(
            read_bed_variants(path, 37, picked, 25), genotype[:, picked]
        )

    def test_rejects_bad_magic_and_out_of_range(self) -> None:
        path = Path(tempfile.mkdtemp()) / "bad.bed"
        path.write_bytes(b"\x00\x00\x00" + b"\x00" * 16)
        with self.assertRaises(ValueError):
            read_bed_variants(path, 4, np.array([0]), 1)
        good = Path(tempfile.mkdtemp()) / "t.bed"
        write_bed(good, np.zeros((4, 2), dtype=np.int8))
        with self.assertRaises(ValueError):
            read_bed_variants(good, 4, np.array([5]), 2)


class TestPenaltyMask(unittest.TestCase):
    """The defect this package exists to avoid: a shared penalty on covariates."""

    def setUp(self) -> None:
        rng = np.random.default_rng(11)
        self.n = 600
        self.cov = rng.normal(size=(self.n, 4))
        self.noise_genotype = rng.normal(size=(self.n, 60))
        self.y = self.cov @ np.array([1.5, -1.0, 0.7, 0.3]) + 0.4 * rng.normal(size=self.n)
        self.train = np.arange(400)
        self.test = np.arange(400, 600)

    def test_unpenalised_covariates_degenerate_onto_the_covariate_only_fit(self) -> None:
        design = np.column_stack([self.cov, self.noise_genotype])
        mask = np.concatenate([np.zeros(4), np.ones(60)])
        big = PrimalRidge(design[self.train], self.y[self.train], mask)
        covariate_only = PrimalRidge(
            self.cov[self.train], self.y[self.train], np.zeros(4)
        )
        np.testing.assert_allclose(
            big.predict(design[self.test], 1e12),
            covariate_only.predict(self.cov[self.test], 1.0),
            rtol=1e-5,
            atol=1e-5,
        )

    def test_a_shared_penalty_makes_the_genetic_arm_lose_on_null_features(self) -> None:
        """Documents the failure mode, so a regression cannot pass silently."""
        design = np.column_stack([self.cov, self.noise_genotype])
        shared = PrimalRidge(design[self.train], self.y[self.train])
        shared_alpha, _ = shared.select_alpha(design[self.test], self.y[self.test], ALPHAS)
        masked = PrimalRidge(
            design[self.train], self.y[self.train], np.concatenate([np.zeros(4), np.ones(60)])
        )
        masked_alpha, _ = masked.select_alpha(design[self.test], self.y[self.test], ALPHAS)
        train_mean = float(self.y[self.train].mean())
        r2_shared = r2_against_train_mean(
            self.y[self.test], shared.predict(design[self.test], shared_alpha), train_mean
        )
        r2_masked = r2_against_train_mean(
            self.y[self.test], masked.predict(design[self.test], masked_alpha), train_mean
        )
        self.assertGreater(r2_masked, r2_shared)

    def test_penalty_mask_shape_is_checked(self) -> None:
        with self.assertRaises(ValueError):
            PrimalRidge(self.cov[self.train], self.y[self.train], np.zeros(3))


class TestRidge(unittest.TestCase):
    def test_recovers_a_linear_signal(self) -> None:
        rng = np.random.default_rng(3)
        x = rng.normal(size=(500, 8))
        y = x @ np.array([2.0, -1.0, 0, 0, 0, 0, 0, 0]) + 0.1 * rng.normal(size=500)
        train, test = np.arange(350), np.arange(350, 500)
        model = PrimalRidge(x[train], y[train])
        alpha, _ = model.select_alpha(x[test], y[test], ALPHAS)
        r2 = r2_against_train_mean(y[test], model.predict(x[test], alpha), float(y[train].mean()))
        self.assertGreater(r2, 0.9)

    def test_cost_is_independent_of_sample_count(self) -> None:
        """Primal form: the normal equations are p-by-p, never n-by-n."""
        rng = np.random.default_rng(4)
        x = rng.normal(size=(20000, 12))
        y = rng.normal(size=20000)
        model = PrimalRidge(x, y)
        self.assertEqual(model.xtx.shape, (12, 12))
        self.assertEqual(model.n, 20000)

    def test_tie_takes_the_larger_alpha(self) -> None:
        rng = np.random.default_rng(5)
        x = rng.normal(size=(120, 3))
        y = rng.normal(size=120)
        model = PrimalRidge(x[:80], y[:80], np.zeros(3))  # unpenalised: every alpha ties
        alpha, _ = model.select_alpha(x[80:], y[80:], ALPHAS)
        self.assertEqual(alpha, max(ALPHAS))

    def test_perturbation_cancels_unperturbed_columns(self) -> None:
        rng = np.random.default_rng(6)
        cov = rng.normal(size=(300, 3))
        dose = rng.normal(size=(300, 10))
        y = cov @ np.array([1.0, 0.5, -0.5]) + dose[:, 0] * 0.8 + 0.2 * rng.normal(size=300)
        design = np.column_stack([cov, dose])
        train, test = np.arange(200), np.arange(200, 300)
        model = PrimalRidge(design[train], y[train], np.concatenate([np.zeros(3), np.ones(10)]))
        perturbed = design.copy()
        perturbed[:, 3:] = dose[train].mean(axis=0)
        delta, score = perturbation_score(model, design[test], perturbed[test], 1.0)
        self.assertGreater(score, 0)
        # Covariates are identical in both designs, so they cannot leak in.
        self.assertLess(abs(float(np.corrcoef(delta, cov[test, 0])[0, 1])), 0.5)
        self.assertGreater(abs(float(np.corrcoef(delta, dose[test, 0])[0, 1])), 0.9)


class TestSelectionAndSplits(unittest.TestCase):
    def test_uniform_thin(self) -> None:
        picked = uniform_thin(np.arange(500, 900), 32)
        self.assertEqual(len(picked), 32)
        self.assertTrue(np.all(np.diff(picked) > 0))
        self.assertEqual(picked[0], 500)
        self.assertEqual(picked[-1], 899)

    def test_deterministic_order_is_input_order_free(self) -> None:
        keys = [f"G{i}" for i in range(40)]
        self.assertEqual(
            deterministic_order(keys, "s"), deterministic_order(list(reversed(keys)), "s")
        )

    def test_groups_never_straddle_a_split(self) -> None:
        ids = [f"S{i}" for i in range(400)]
        groups = [f"FAM{i // 4}" for i in range(400)]  # quads
        labels = assign_splits(
            ids, groups, [""] * 400, {"train": 0.7, "validation": 0.15, "evaluation": 0.15}, 7
        )
        by_group: dict[str, set[str]] = {}
        for g, l in zip(groups, labels):
            by_group.setdefault(g, set()).add(l)
        self.assertTrue(all(len(v) == 1 for v in by_group.values()))
        self.assertEqual(set(labels), {"train", "validation", "evaluation"})

    def test_supplied_splits_are_honoured_and_validated(self) -> None:
        ids = ["a", "b"]
        labels = assign_splits(ids, ids, ["train", "evaluation"], {"train": 1.0}, 1)
        self.assertEqual(list(labels), ["train", "evaluation"])
        with self.assertRaises(ValueError):
            assign_splits(ids, ids, ["train", "typo"], {"train": 1.0}, 1)

    def test_af_imputation_uses_train_frequency(self) -> None:
        genotype = np.array([[0, 2], [2, 2], [1, -1], [-1, 0]], dtype=np.int8)
        train = np.array([0, 1, 2])
        af, missing, dosage = af_missingness_dosage(genotype, train)
        self.assertAlmostEqual(af[0], 0.5)
        self.assertAlmostEqual(missing[1], 1.0 / 3.0)
        self.assertAlmostEqual(float(dosage[2, 1]), 2.0 * af[1], places=5)


class TestStatistics(unittest.TestCase):
    def _panel(self, noise_b: float, seed: int = 8):
        rng = np.random.default_rng(seed)
        genes, people = 30, 200
        truth = rng.normal(size=people)
        predictions = {
            "A": 0.1 * truth + 0.9 * rng.normal(size=(genes, people)),
            "B": 0.6 * truth + noise_b * rng.normal(size=(genes, people)),
        }
        return truth, predictions

    def test_shared_outcome_bootstrap_detects_a_better_arm(self) -> None:
        truth, predictions = self._panel(0.2)
        out = clustered_macro_bootstrap_shared(
            truth, predictions, 0.0, ("B", "A"), 1, 300
        )
        self.assertTrue(out["pass"])
        self.assertEqual(out["n_individuals"], 200)

    def test_identical_arms_cannot_pass(self) -> None:
        truth, predictions = self._panel(0.5)
        predictions["B"] = predictions["A"].copy()
        out = clustered_macro_bootstrap_shared(truth, predictions, 0.0, ("B", "A"), 1, 300)
        self.assertAlmostEqual(out["macro_delta_r2"], 0.0, places=12)
        self.assertFalse(out["pass"])

    def test_rank_helpers(self) -> None:
        np.testing.assert_allclose(rankdata_average(np.array([5.0, 9.0, 9.0])), [1.0, 2.5, 2.5])
        x = np.arange(15.0)
        self.assertAlmostEqual(spearman(x, x), 1.0, places=10)
        self.assertAlmostEqual(spearman(x, -x), -1.0, places=10)
        strata = decile_strata(np.random.default_rng(1).normal(size=100))
        self.assertGreaterEqual(strata.min(), 0)
        self.assertLessEqual(strata.max(), 9)

    def test_matched_null_separates_signal_from_noise(self) -> None:
        rng = np.random.default_rng(9)
        target = rng.normal(size=150)
        strata = np.zeros(150, dtype=np.int64)
        self.assertLess(
            matched_null_spearman(target + 0.05 * rng.normal(size=150), target, strata, 1, 400)[
                "p_one_sided_greater"
            ],
            0.01,
        )
        self.assertGreater(
            matched_null_spearman(rng.normal(size=150), target, strata, 2, 400)["p_two_sided"],
            0.05,
        )

    def test_gene_bootstrap_interval_brackets_the_estimate(self) -> None:
        rng = np.random.default_rng(10)
        target = rng.normal(size=120)
        score = target + 0.5 * rng.normal(size=120)
        out = gene_bootstrap_spearman(score, target, 3, 400)
        low, high = out["bootstrap_ci95"]
        self.assertLessEqual(low, out["observed_spearman"])
        self.assertGreaterEqual(high, out["observed_spearman"])
        self.assertGreater(low, 0.0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
