"""Do the Level 0 estimators actually separate the four worlds?"""

import unittest
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from subspace_conservation.estimators import (  # noqa: E402
    excess_joint_conservation,
    observed_distances,
    orientation_diagnostic,
    paired_bootstrap_ci,
    standardised_excess,
    tail_quantiles,
    level0_verdict,
)
from subspace_conservation.nulls import parametric_null  # noqa: E402
from subspace_conservation.simulate import (  # noqa: E402
    TruthSpec,
    apply_distortion,
    draw_truth,
    make_decorative_matched,
    observe,
)


def _world(kind, n_pairs, rng, spec=None, cond="gl_8"):
    spec = spec or TruthSpec(r=8, sigma1=8.0, sigma2=4.0, noise=1.0)
    A = draw_truth(n_pairs, spec, rng)
    if kind == "S":
        B = apply_distortion(A, cond, rng)
    elif kind == "C":
        B = A.copy()
    elif kind == "N":
        B = draw_truth(n_pairs, spec, rng)
    elif kind == "D":
        B = make_decorative_matched(A, apply_distortion(A, cond, rng), rng)
    else:
        raise ValueError(kind)
    EA, SA = observe(A, spec.noise, rng)
    EB, SB = observe(B, spec.noise, rng)
    return EA, SA, EB, SB


class TestStandardisedExcess(unittest.TestCase):
    def test_subspace_beats_coordinates_in_the_conserved_world(self):
        rng = np.random.default_rng(1)
        EA, SA, EB, SB = _world("S", 200, rng)
        null = parametric_null([EA, EB], [SA, SB], 40, rng)
        sub, coord = observed_distances(EA, EB)
        i_sub = standardised_excess(sub, null.sub, "subspace").index
        i_coord = standardised_excess(coord, null.coord, "coordinate").index
        self.assertLess(i_sub, i_coord)

    def test_nothing_is_conserved_in_the_unrelated_world(self):
        rng = np.random.default_rng(2)
        EA, SA, EB, SB = _world("N", 200, rng)
        null = parametric_null([EA, EB], [SA, SB], 40, rng)
        sub, coord = observed_distances(EA, EB)
        self.assertGreater(standardised_excess(sub, null.sub, "subspace").index, 2.0)
        self.assertGreater(standardised_excess(coord, null.coord, "coordinate").index, 2.0)

    def test_fully_conserved_world_sits_near_its_null(self):
        rng = np.random.default_rng(3)
        EA, SA, EB, SB = _world("C", 200, rng)
        null = parametric_null([EA, EB], [SA, SB], 40, rng)
        sub, coord = observed_distances(EA, EB)
        self.assertLess(abs(standardised_excess(sub, null.sub, "subspace").index), 1.0)
        self.assertLess(abs(standardised_excess(coord, null.coord, "coordinate").index), 1.0)

    def test_index_is_an_effect_size_not_a_significance(self):
        # Doubling the pairs must not inflate it; that is the whole reason the
        # denominator is the across-pair spread and not a standard error.
        rng = np.random.default_rng(4)
        vals = []
        for n in (150, 600):
            EA, SA, EB, SB = _world("S", n, np.random.default_rng(99))
            null = parametric_null([EA, EB], [SA, SB], 30, rng)
            sub, _ = observed_distances(EA, EB)
            vals.append(standardised_excess(sub, null.sub, "subspace").index)
        self.assertLess(abs(vals[0] - vals[1]), 1.5)


class TestDeflation(unittest.TestCase):
    """The test that decides whether the pair framing earns its keep."""

    def test_detects_excess_in_the_conserved_world(self):
        rng = np.random.default_rng(5)
        EA, _, EB, _ = _world("S", 150, rng)
        out = excess_joint_conservation(EA, EB, 60, rng, n_boot=800)
        self.assertEqual(out["verdict"], "EXCESS_JOINT_CONSERVATION")

    def test_finds_no_excess_in_the_decorative_world(self):
        # World D reproduces every single-SNP marginal of world S and keeps no
        # joint structure. An estimator that calls this positive cannot support
        # a pair-level claim at all.
        rng = np.random.default_rng(6)
        EA, _, EB, _ = _world("D", 150, rng)
        out = excess_joint_conservation(EA, EB, 60, rng, n_boot=800)
        self.assertEqual(out["verdict"], "NO_EXCESS_JOINT_CONSERVATION")


class TestPlumbing(unittest.TestCase):
    def test_tail_quantiles_are_probabilities(self):
        rng = np.random.default_rng(8)
        d_ctx = rng.random(40)
        d_null = rng.random((40, 25))
        q = tail_quantiles(d_ctx, d_null)
        self.assertTrue(np.all((q > 0) & (q <= 1.0)))

    def test_paired_bootstrap_recovers_a_known_shift(self):
        rng = np.random.default_rng(9)
        a = rng.standard_normal(400)
        b = a - 0.5
        diff, lo, hi = paired_bootstrap_ci(a, b, 2000, rng)
        self.assertAlmostEqual(diff, 0.5, places=6)
        self.assertLess(lo, 0.5)
        self.assertGreater(hi, 0.5)

    def test_orientation_diagnostic_flags_flipped_alleles(self):
        rng = np.random.default_rng(10)
        EA = rng.standard_normal((100, 6, 2))
        EB = EA.copy()
        EB[:30, :, 0] *= -1.0
        out = orientation_diagnostic(EA, EB)
        self.assertTrue(out["blocking"])
        self.assertGreater(out["flip_fraction"], 0.1)

    def test_orientation_diagnostic_passes_clean_input(self):
        rng = np.random.default_rng(12)
        EA = rng.standard_normal((100, 6, 2))
        EB = EA + 0.05 * rng.standard_normal((100, 6, 2))
        self.assertFalse(orientation_diagnostic(EA, EB)["blocking"])

    def test_verdict_separates_all_four_failure_modes(self):
        # Arguments are conservation fractions in [0, 1]: 1 = moved no more than
        # noise forces, 0 = moved as far as an unrelated pair.
        joint_yes = {"verdict": "EXCESS_JOINT_CONSERVATION"}
        joint_no = {"verdict": "NO_EXCESS_JOINT_CONSERVATION"}

        # plane conserved, coordinates not, and the pair carries joint structure
        self.assertEqual(level0_verdict(0.9, 0.1, 0.2, joint_yes)["verdict"], "PROCEED_TO_LEVEL_1")
        # same, but the plane is no more stable than its two columns already imply
        self.assertEqual(
            level0_verdict(0.9, 0.1, 0.2, joint_no)["verdict"], "CONSERVED_BUT_DECORATIVE_STOP"
        )
        # coordinates survive too: nothing to quotient out, route unmotivated
        self.assertEqual(
            level0_verdict(0.9, 0.9, -0.01, joint_yes)["verdict"],
            "NO_CONTEXT_DISTORTION_ROUTE_UNMOTIVATED",
        )
        # nothing survives. Must not read as conservation just because the
        # subspace metric saturates before the coordinate one does.
        self.assertEqual(
            level0_verdict(0.02, 0.01, 0.2, joint_yes)["verdict"],
            "NO_SUBSPACE_CONSERVATION_STOP",
        )


if __name__ == "__main__":
    unittest.main()
