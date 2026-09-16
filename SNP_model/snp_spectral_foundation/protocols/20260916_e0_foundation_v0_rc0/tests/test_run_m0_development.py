from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "run_m0_development.py"
SPEC = importlib.util.spec_from_file_location("run_m0_development", SCRIPT)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot load run_m0_development.py")
M0 = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = M0
SPEC.loader.exec_module(M0)


class M0DevelopmentTests(unittest.TestCase):
    def test_gt_to_dosage_rejects_missing_haploid_and_non_biallelic(self) -> None:
        self.assertEqual(M0.gt_to_dosage((0, 1)), 1)
        self.assertEqual(M0.gt_to_dosage((1, 1)), 2)
        self.assertIsNone(M0.gt_to_dosage((0, None)))
        self.assertIsNone(M0.gt_to_dosage((1,)))
        self.assertIsNone(M0.gt_to_dosage((0, 2)))

    def test_empirical_and_hwe_probabilities_are_normalized(self) -> None:
        empirical = M0.empirical_probabilities(np.asarray((8.0, 4.0, 2.0)), 0.5)
        hwe = M0.hwe_probabilities(0.25)
        self.assertAlmostEqual(float(empirical.sum()), 1.0)
        self.assertAlmostEqual(float(hwe.sum()), 1.0)
        np.testing.assert_allclose(hwe, np.asarray((0.5625, 0.375, 0.0625)))

    def test_mask_rng_seed_is_stable_and_keyed(self) -> None:
        first = M0.stable_rng_seed(7, "sample", "block")
        self.assertEqual(first, M0.stable_rng_seed(7, "sample", "block"))
        self.assertNotEqual(first, M0.stable_rng_seed(8, "sample", "block"))
        self.assertNotEqual(first, M0.stable_rng_seed(7, "other", "block"))


if __name__ == "__main__":
    unittest.main()
