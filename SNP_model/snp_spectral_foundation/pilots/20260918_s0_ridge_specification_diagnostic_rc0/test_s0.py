"""Contract tests for the S0 ridge-specification diagnostic. No assets required."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

PILOT_DIR = Path(__file__).resolve().parent
TQG1_DIR = PILOT_DIR.parent / "20260917_tqg1_gene_layer_readout_rc0"
for path in (str(PILOT_DIR), str(TQG1_DIR)):
    if path not in sys.path:
        sys.path.insert(0, path)

from tqg1_core import fit_ridge  # noqa: E402

from s0_core import (  # noqa: E402
    SPECIFICATIONS,
    cross_validated_predictions,
    fit_arm,
    group_kfold,
    standardize,
)

ALPHAS = [1e-6, 1e-4, 1e-2, 1.0, 10.0, 100.0, 1e4, 1e6]


def synthetic(n=150, n_cov=12, n_dos=256, signal="covariate", seed=11):
    """Covariates carry the signal; dosage is pure noise unless asked otherwise."""
    rng = np.random.default_rng(seed)
    covariates = rng.normal(size=(n, n_cov))
    dosage = rng.normal(size=(n, n_dos))
    if signal == "covariate":
        y = covariates @ rng.normal(size=n_cov) + 0.25 * rng.normal(size=n)
    else:
        y = dosage[:, :4] @ rng.normal(size=4) * 3.0 + 0.25 * rng.normal(size=n)
    families = [f"F{i // 3:03d}" for i in range(n)]
    return covariates, dosage, y, families


class TestStandardize(unittest.TestCase):
    def test_constant_columns_are_dropped(self):
        x = np.column_stack([np.ones(20), np.arange(20.0)])
        z = standardize(x, np.arange(20))
        self.assertEqual(z.shape[1], 1)

    def test_scaling_uses_only_the_fitting_rows(self):
        x = np.arange(40.0).reshape(20, 2)
        fit_idx = np.arange(10)
        z = standardize(x, fit_idx)
        self.assertAlmostEqual(float(z[fit_idx].mean()), 0.0, places=10)
        self.assertAlmostEqual(float(z[fit_idx, 0].std()), 1.0, places=10)
        self.assertGreater(float(z[10:].mean()), 0.0)

    def test_all_constant_is_rejected(self):
        with self.assertRaises(ValueError):
            standardize(np.ones((10, 3)), np.arange(10))


class TestS1Fidelity(unittest.TestCase):
    """S1 must reproduce the fit the programme has always used."""

    def test_s1_matches_tqg1_fit_ridge(self):
        covariates, dosage, y, _ = synthetic()
        fit_idx = np.arange(100)
        tune_idx = np.arange(100, 130)
        test_idx = np.arange(130, 150)
        design = np.column_stack([covariates, dosage])
        reference, info = fit_ridge(
            design, y[fit_idx], fit_idx, tune_idx, ALPHAS, y[tune_idx]
        )
        predict, alpha = fit_arm(
            "S1", covariates, dosage, y, fit_idx, tune_idx, ALPHAS
        )
        self.assertEqual(alpha, info["alpha"])
        np.testing.assert_allclose(predict(test_idx), reference(design[test_idx]), rtol=1e-9)

    def test_s1_arm_a_matches_tqg1_fit_ridge(self):
        covariates, _, y, _ = synthetic()
        fit_idx, tune_idx, test_idx = np.arange(100), np.arange(100, 130), np.arange(130, 150)
        reference, _ = fit_ridge(covariates, y[fit_idx], fit_idx, tune_idx, ALPHAS, y[tune_idx])
        predict, _ = fit_arm("S1", covariates, None, y, fit_idx, tune_idx, ALPHAS)
        np.testing.assert_allclose(predict(test_idx), reference(covariates[test_idx]), rtol=1e-9)


class TestUnpenalizedSpecifications(unittest.TestCase):
    def test_s2_and_s3_agree_on_arm_a(self):
        covariates, _, y, _ = synthetic()
        fit_idx, tune_idx, test_idx = np.arange(100), np.arange(100, 130), np.arange(130, 150)
        p2, _ = fit_arm("S2", covariates, None, y, fit_idx, tune_idx, ALPHAS)
        p3, _ = fit_arm("S3", covariates, None, y, fit_idx, tune_idx, ALPHAS)
        np.testing.assert_allclose(p2(test_idx), p3(test_idx), rtol=1e-10)

    def test_s2_arm_a_is_ordinary_least_squares(self):
        covariates, _, y, _ = synthetic()
        fit_idx = np.arange(120)
        predict, alpha = fit_arm("S2", covariates, None, y, fit_idx, fit_idx, ALPHAS)
        design = np.column_stack([np.ones(len(fit_idx)), standardize(covariates, fit_idx)[fit_idx]])
        expected = design @ np.linalg.lstsq(design, y[fit_idx], rcond=None)[0]
        self.assertEqual(alpha, 0.0)
        np.testing.assert_allclose(predict(fit_idx), expected, rtol=1e-9)

    def test_tiny_alpha_recovers_the_joint_least_squares_fit(self):
        covariates, dosage, y, _ = synthetic(n=60, n_cov=4, n_dos=6)
        fit_idx = np.arange(60)
        predict, _ = fit_arm("S2", covariates, dosage, y, fit_idx, fit_idx, [1e-10])
        design = np.column_stack(
            [
                np.ones(60),
                standardize(covariates, fit_idx)[fit_idx],
                standardize(dosage, fit_idx)[fit_idx],
            ]
        )
        expected = design @ np.linalg.lstsq(design, y, rcond=None)[0]
        # alpha is 1e-10, not 0, and the dual solve carries its own conditioning,
        # so agreement is to ~1e-5 on fitted values of order 1. 1e-3 still pins
        # the partial ridge to the joint least-squares fit.
        np.testing.assert_allclose(predict(fit_idx), expected, atol=1e-3)

    def test_every_specification_is_reachable(self):
        covariates, dosage, y, _ = synthetic()
        for specification in SPECIFICATIONS:
            predict, _ = fit_arm(
                specification, covariates, dosage, y, np.arange(100), np.arange(100, 130), ALPHAS
            )
            self.assertEqual(predict(np.arange(130, 150)).shape, (20,))

    def test_unknown_specification_is_rejected(self):
        covariates, dosage, y, _ = synthetic()
        with self.assertRaises(ValueError):
            fit_arm("S9", covariates, dosage, y, np.arange(100), np.arange(100, 130), ALPHAS)


class TestTheDilutionMechanism(unittest.TestCase):
    """The claim under test: a shared penalty lets noise columns destroy arm A's fit."""

    def _cv_r2(self, specification, covariates, dosage, y, folds, families):
        prediction, _ = cross_validated_predictions(
            specification, covariates, dosage, y, folds, families, ALPHAS, 7
        )
        return 1.0 - np.sum((y - prediction) ** 2) / np.sum((y - y.mean()) ** 2)

    def test_shared_penalty_loses_the_covariate_fit_and_partial_ridge_keeps_it(self):
        covariates, dosage, y, families = synthetic(signal="covariate")
        folds = group_kfold(families, 5, 3)
        damage = {}
        for specification in SPECIFICATIONS:
            only_a = self._cv_r2(specification, covariates, None, y, folds, families)
            with_b = self._cv_r2(specification, covariates, dosage, y, folds, families)
            damage[specification] = only_a - with_b
        self.assertGreater(damage["S1"], 0.05)
        self.assertLess(damage["S2"], damage["S1"] / 2)
        self.assertLess(damage["S3"], damage["S1"] / 2)

    def test_real_dosage_signal_is_still_found_by_every_specification(self):
        covariates, dosage, y, families = synthetic(n_dos=32, signal="dosage")
        folds = group_kfold(families, 5, 3)
        for specification in SPECIFICATIONS:
            only_a = self._cv_r2(specification, covariates, None, y, folds, families)
            with_b = self._cv_r2(specification, covariates, dosage, y, folds, families)
            self.assertGreater(with_b, only_a, f"{specification} lost real dosage signal")


class TestFolds(unittest.TestCase):
    def test_families_are_never_split(self):
        families = [f"F{i // 3:03d}" for i in range(150)]
        folds = group_kfold(families, 5, 3)
        for family in set(families):
            rows = [i for i, f in enumerate(families) if f == family]
            self.assertEqual(len({int(folds[i]) for i in rows}), 1)

    def test_assignment_is_deterministic(self):
        families = [f"F{i // 3:03d}" for i in range(150)]
        np.testing.assert_array_equal(group_kfold(families, 5, 3), group_kfold(families, 5, 3))

    def test_every_fold_is_populated(self):
        families = [f"F{i // 3:03d}" for i in range(150)]
        folds = group_kfold(families, 5, 3)
        self.assertEqual(sorted(set(folds.tolist())), [0, 1, 2, 3, 4])

    def test_too_few_groups_is_rejected(self):
        with self.assertRaises(ValueError):
            group_kfold(["A", "A", "B"], 5, 3)


class TestHeldOutHygiene(unittest.TestCase):
    def test_a_test_rows_outcome_never_reaches_its_own_prediction(self):
        covariates, dosage, y, families = synthetic()
        folds = group_kfold(families, 5, 3)
        baseline, _ = cross_validated_predictions(
            "S2", covariates, dosage, y, folds, families, ALPHAS, 7
        )
        corrupted = y.copy()
        held_out = np.flatnonzero(folds == 0)
        corrupted[held_out] += 1000.0
        after, _ = cross_validated_predictions(
            "S2", covariates, dosage, corrupted, folds, families, ALPHAS, 7
        )
        np.testing.assert_allclose(baseline[held_out], after[held_out], rtol=1e-9)

    def test_every_row_is_predicted_exactly_once_out_of_fold(self):
        covariates, dosage, y, families = synthetic()
        folds = group_kfold(families, 5, 3)
        prediction, alphas = cross_validated_predictions(
            "S1", covariates, dosage, y, folds, families, ALPHAS, 7
        )
        self.assertTrue(np.isfinite(prediction).all())
        self.assertEqual(len(alphas), 5)


class TestRoleGuard(unittest.TestCase):
    def test_forbidden_roles_abort_the_run(self):
        import run_s0

        config = {
            "roles_used": ["dev_train", "dev_validation"],
            "forbidden_roles": ["task_gate", "bridge_test"],
        }
        clean = [
            {"development_role": "dev_train", "primary_role": "development", "family": "F1"},
            {"development_role": "dev_validation", "primary_role": "development", "family": "F2"},
        ]
        dev_idx, families = run_s0.resolve_development_rows(clean, config)
        self.assertEqual(len(dev_idx), 2)
        self.assertEqual(families, ["F1", "F2"])

        leaking = clean + [
            {"development_role": "dev_train", "primary_role": "bridge_test", "family": "F3"}
        ]
        with self.assertRaises(SystemExit):
            run_s0.resolve_development_rows(leaking, config)


if __name__ == "__main__":
    unittest.main()
