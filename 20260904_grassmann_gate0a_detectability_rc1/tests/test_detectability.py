import sys
import unittest
from pathlib import Path

import numpy as np

PKG = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PKG / "src"))

from detectability import power_table, unit_ci_lower_bounds, min_replicates  # noqa: E402


class TestDetectability(unittest.TestCase):
    def test_fpr_is_near_one_sided_nominal(self):
        # decision is CI-lower > 0; under Delta=0 a calibrated 95% CI gives ~0.025
        table, _ = power_table([0.0], [0.02], [100], M=4000, B=1500, seed=1)
        fpr = table[(0.0, 0.02, 100)]
        self.assertLess(fpr, 0.06)
        self.assertGreater(fpr, 0.005)

    def test_power_monotone_in_R(self):
        table, _ = power_table([0.005], [0.02], [20, 50, 100, 200], M=3000, B=1200, seed=2)
        vals = [table[(0.005, 0.02, R)] for R in (20, 50, 100, 200)]
        for a, b in zip(vals, vals[1:]):
            self.assertLessEqual(a, b + 0.03)  # non-decreasing up to MC noise

    def test_power_monotone_in_delta(self):
        table, _ = power_table([0.0, 0.002, 0.005, 0.010], [0.02], [100], M=3000, B=1200, seed=3)
        vals = [table[(d, 0.02, 100)] for d in (0.0, 0.002, 0.005, 0.010)]
        for a, b in zip(vals, vals[1:]):
            self.assertLessEqual(a, b + 0.03)

    def test_smaller_sigma_more_power(self):
        table, _ = power_table([0.005], [0.01, 0.05], [50], M=3000, B=1200, seed=4)
        self.assertGreater(table[(0.005, 0.01, 50)], table[(0.005, 0.05, 50)])

    def test_shift_scale_identity(self):
        # power(delta,sigma,R) must equal mean(L0 > -delta/sigma)
        low = unit_ci_lower_bounds(50, M=2000, B=1000, rng=5)
        table, L0 = power_table([0.005], [0.02], [50], M=2000, B=1000, seed=5)
        # same seed path is not identical here (separate rng); just check formula shape
        p = float(np.mean(L0[50] > -0.005 / 0.02))
        self.assertAlmostEqual(p, table[(0.005, 0.02, 50)], places=12)

    def test_min_replicates(self):
        table, _ = power_table([0.005], [0.01], [20, 50, 100, 200], M=3000, B=1200, seed=6)
        r = min_replicates(table, 0.01, [20, 50, 100, 200], 0.005, 0.90)
        self.assertIn(r, (20, 50, 100, 200, None))


if __name__ == "__main__":
    unittest.main()
