from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
CONFIG = json.loads((ROOT / "config" / "BOUNDED_SMOKE_CONFIG.json").read_text(encoding="utf-8"))
META = json.loads((ROOT / "PARENT_EVIDENCE.json").read_text(encoding="utf-8"))
RANK_PARENT = ROOT.parent / META["code_parent"]
R1_PARENT = RANK_PARENT.parent / json.loads(
    (RANK_PARENT / "PARENT_EVIDENCE.json").read_text(encoding="utf-8")
)["code_parent"]
GC_PARENT = R1_PARENT.parent / json.loads(
    (R1_PARENT / "PARENT_EVIDENCE.json").read_text(encoding="utf-8")
)["code_parent"]
sys.path[:0] = [str(RANK_PARENT / "src"), str(R1_PARENT / "src"), str(GC_PARENT / "src")]

from rank_gate_r1_5 import run_rank_gated_maxT
from shared_family_r1 import generate_shared_family, independent_target_selection


class DesignContracts(unittest.TestCase):
    def test_frozen_counts_and_seed_namespaces(self):
        self.assertEqual(CONFIG["planned_counts"]["independent_families"], 21)
        self.assertEqual(CONFIG["planned_counts"]["family_level_resamples"], 819)
        seeds = {
            CONFIG["data_seed_base"],
            CONFIG["selection_seed_base"],
            CONFIG["inference_seed_base"],
            CONFIG["bootstrap_seed_base"],
        }
        self.assertEqual(len(seeds), 4)
        self.assertGreaterEqual(min(seeds), 628000000)

    def test_selected_target_is_subject_disjoint_and_applied(self):
        selected = independent_target_selection(
            selection_seed=CONFIG["selection_seed_base"],
            inference_seed=CONFIG["inference_seed_base"],
            n_inference=120,
        )
        self.assertFalse(set(selected.selection_subject_ids) & set(selected.inference_subject_ids))
        family = generate_shared_family(
            seed=CONFIG["inference_seed_base"] + 10_000,
            family_size=4,
            g_override=selected.inference_selected_g,
            subject_ids_override=selected.inference_subject_ids,
        )
        self.assertTrue(
            np.array_equal(family.g, selected.inference_target_panel[:, selected.selected_index])
        )

    def test_rank_gate_preserves_family_and_maps_every_resample(self):
        family = generate_shared_family(
            seed=CONFIG["data_seed_base"],
            n=120,
            family_size=CONFIG["family_size"],
            conditional_ld=True,
            heteroskedastic=True,
        )
        result = run_rank_gated_maxT(
            family,
            resamples=19,
            seed=CONFIG["bootstrap_seed_base"],
            rank=CONFIG["rank"],
            ridge_lambda=CONFIG["ridge_lambda"],
            minimum_gap=CONFIG["minimum_fitted_rank_gap"],
        )
        self.assertEqual(len(result.observed), CONFIG["family_size"])
        self.assertEqual(result.resampled.shape, (CONFIG["family_size"], 19))
        self.assertTrue(np.all(result.observed[~result.observed_eligible] == 0))
        self.assertTrue(np.all(result.candidate_p_values[~result.observed_eligible] == 1))
        self.assertTrue(np.all(result.resampled[~result.resampled_eligible] == 0))
        self.assertEqual(len(result.multiplier_fingerprints), 19)

    def test_design_cannot_authorize_its_own_run(self):
        self.assertIn("bounded_smoke_rc2_run_without_new_approval", CONFIG["does_not_authorize"])
        source = (ROOT / "scripts" / "run_bounded_smoke.py").read_text(encoding="utf-8")
        self.assertIn("APPROVE_BOUNDED_SMOKE_RC2_RUN", source)
        self.assertIn("package_manifest_sha256", source)


if __name__ == "__main__":
    unittest.main()
