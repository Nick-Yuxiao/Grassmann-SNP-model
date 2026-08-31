from __future__ import annotations

import hashlib
import json
import sys
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
CONFIG = json.loads((ROOT / "config" / "FEASIBILITY_CONFIG.json").read_text(encoding="utf-8"))
META = json.loads((ROOT / "PARENT_EVIDENCE.json").read_text(encoding="utf-8"))
RANK_PARENT = ROOT.parent / META["implementation_reference"]
R1_PARENT = RANK_PARENT.parent / json.loads(
    (RANK_PARENT / "PARENT_EVIDENCE.json").read_text(encoding="utf-8")
)["code_parent"]
GC_PARENT = R1_PARENT.parent / json.loads(
    (R1_PARENT / "PARENT_EVIDENCE.json").read_text(encoding="utf-8")
)["code_parent"]
sys.path[:0] = [str(ROOT / "src"), str(GC_PARENT / "src")]

from design_feasibility import generate_scenario, fit_scenario, principal_angles_deg, truth_matrices
from grassmann_v6_1.core import conditional_matrices, geometry_scores, top_subspace


class FeasibilityContracts(unittest.TestCase):
    def test_truth_construction_has_requested_gap_and_angle(self):
        b, gamma = truth_matrices(
            population_relative_gap=0.10,
            true_max_principal_angle_deg=20.0,
            effect_scale=CONFIG["effect_scale"],
        )
        matrices = conditional_matrices(b, gamma)
        geometry = geometry_scores(matrices, CONFIG["rank"])
        bases = [top_subspace(matrix, CONFIG["rank"])[0] for matrix in matrices]
        angle = principal_angles_deg(bases[0], bases[2]).max()
        self.assertAlmostEqual(min(geometry["relative_rank_gaps"]), 0.10, places=12)
        self.assertAlmostEqual(angle, 20.0, places=10)

    def test_same_seed_is_reproducible_and_angle_paired(self):
        kwargs = dict(
            seed=CONFIG["seed_base"],
            n=500,
            maf=0.20,
            population_relative_gap=0.10,
            effect_scale=CONFIG["effect_scale"],
            residual_sd=CONFIG["residual_sd"],
            conditional_ld_rhos_by_dosage=tuple(CONFIG["conditional_ld_rhos_by_dosage"]),
            residual_scales_by_dosage=tuple(CONFIG["residual_scales_by_dosage"]),
        )
        null = generate_scenario(true_max_principal_angle_deg=0.0, **kwargs)
        alternative = generate_scenario(true_max_principal_angle_deg=20.0, **kwargs)
        repeat = generate_scenario(true_max_principal_angle_deg=0.0, **kwargs)
        self.assertEqual(hashlib.sha256(null.g.tobytes()).hexdigest(), hashlib.sha256(repeat.g.tobytes()).hexdigest())
        self.assertTrue(np.array_equal(null.g, alternative.g))
        self.assertTrue(np.allclose(null.x, alternative.x))
        self.assertFalse(np.allclose(null.y, alternative.y))

    def test_low_maf_small_n_is_explicit_no_call(self):
        scenario = generate_scenario(
            seed=CONFIG["seed_base"] + 1,
            n=500,
            maf=0.05,
            population_relative_gap=0.20,
            true_max_principal_angle_deg=20.0,
            effect_scale=CONFIG["effect_scale"],
            residual_sd=CONFIG["residual_sd"],
            conditional_ld_rhos_by_dosage=tuple(CONFIG["conditional_ld_rhos_by_dosage"]),
            residual_scales_by_dosage=tuple(CONFIG["residual_scales_by_dosage"]),
        )
        fitted = fit_scenario(
            scenario,
            rank=CONFIG["rank"],
            ridge_lambda=CONFIG["ridge_lambda"],
            minimum_group_count=CONFIG["minimum_group_count"],
            minimum_fitted_rank_gap=CONFIG["minimum_fitted_rank_gap"],
        )
        self.assertFalse(fitted.group_count_eligible)
        self.assertFalse(fitted.jointly_estimable)
        self.assertEqual(fitted.gated_direction_score, 0.0)

    def test_evidence_firewall_is_machine_readable(self):
        firewall = json.loads((ROOT / "EVIDENCE_FIREWALL.json").read_text(encoding="utf-8"))
        self.assertFalse(firewall["formal_evidence_eligible"])
        self.assertFalse(firewall["promotion_permitted"])
        self.assertIn("t14_or_t16_power_ranking", firewall["forbidden_uses"])
        self.assertIn("gpu_or_v7_work", firewall["forbidden_uses"])


if __name__ == "__main__":
    unittest.main()
