from __future__ import annotations

import importlib.util
import math
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
AUDIT_PATH = ROOT / "p0" / "audit_estimand_family_v7_3_0.py"
SPEC = importlib.util.spec_from_file_location("audit_v730", AUDIT_PATH)
AUDIT = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(AUDIT)


def summary(passed: bool, shape: str, drops: tuple[float, float, float]) -> dict:
    return {
        "all_absolute_changes_le_0p002": passed,
        "instability": False,
        "shape_class": shape,
        "changes": [
            {"start": start, "stop": stop, "nll_drop": drop}
            for (start, stop), drop in zip(
                ((34000, 36000), (36000, 38000), (38000, 40000)), drops
            )
        ],
    }


class TestV730(unittest.TestCase):
    def decisions(self) -> tuple[dict, dict]:
        at_30k = {}
        at_40k = {}
        for model in AUDIT.MODELS:
            ld = f"{model}|ld_block_0p90"
            longrange = f"{model}|within_chrom_longrange_0p90"
            at_30k[ld] = summary(True, "STABLE", (0.0007, 0.0006, 0.0005))
            at_40k[ld] = summary(True, "STABLE", (0.0004, 0.0003, 0.0004))
            at_30k[longrange] = summary(False, "NOT_STABLE", (0.003, 0.004, 0.005))
            at_40k[longrange] = summary(False, "NOT_STABLE", (0.005, 0.006, 0.007))
        gpc = "local_attn_gpc_8m_w256|within_chrom_longrange_0p90"
        at_30k[gpc] = summary(True, "STABLE", (0.0010, 0.0012, 0.0014))
        return {"terminal_4e_4": at_30k}, {"terminal_4e_4": at_40k}

    def test_family_fork_uses_out_of_sample_persistence(self) -> None:
        at_30k, at_40k = self.decisions()
        families, false_plateaus = AUDIT.family_eligibility(at_30k, at_40k)
        self.assertEqual(
            families["ld_block_0p90"]["capacity_status"],
            "ELIGIBLE_FOR_CONFIRMATORY_30K",
        )
        self.assertEqual(
            families["ld_block_0p90"]["basis"],
            "OUT_OF_SAMPLE_PERSISTENCE_30K_TO40K",
        )
        self.assertEqual(
            families["within_chrom_longrange_0p90"]["capacity_status"],
            "INCONCLUSIVE_BUDGET_NOT_ESTIMABLE",
        )
        self.assertEqual(len(false_plateaus), 1)
        self.assertIn("local_attn_gpc_8m_w256", false_plateaus[0])

    def test_positive_and_negative_efficiency_upgrades_are_forbidden(self) -> None:
        at_30k, at_40k = self.decisions()
        families, _ = AUDIT.family_eligibility(at_30k, at_40k)
        longrange = families["within_chrom_longrange_0p90"]
        self.assertFalse(longrange["positive_capacity_upgrade_permitted"])
        self.assertFalse(longrange["negative_capacity_upgrade_permitted"])
        self.assertEqual(longrange["registered_efficiency_budget_points"], [20000, 30000, 40000])

    def test_boundary_change_is_explicitly_not_pure_resume(self) -> None:
        extension = {}
        final = {}
        for index, key in enumerate(sorted(AUDIT.EXPECTED_CELLS)):
            extension[key] = {"rows": [(30000, 0.75 + index / 1000)]}
            final[key] = {"rows": [(30250, 0.749 + index / 1000)]}
        result = AUDIT.boundary_audit(extension, final)
        self.assertFalse(result["pure_resume_discontinuity_estimate"])
        self.assertEqual(result["confounded_with_optimizer_steps"], 250)
        self.assertAlmostEqual(result["maximum_absolute_change"], 0.001)

    def test_historical_variance_is_planning_proxy_not_data_seed_pilot(self) -> None:
        grouped = {}
        for key in sorted(AUDIT.EXPECTED_CELLS):
            grouped[key] = [
                {
                    "terminal_tail5": 0.700,
                    "mask_seed": 91001,
                    "init_seed": 81001,
                    "data_seed": None,
                },
                {
                    "terminal_tail5": 0.702,
                    "mask_seed": 91002,
                    "init_seed": 81002,
                    "data_seed": None,
                },
            ]
        proxy = AUDIT.variance_proxy_from_grouped(grouped)
        self.assertEqual(proxy["role"], "PLANNING_PROXY_ONLY")
        self.assertFalse(proxy["replaces_n5_data_seed_pilot"])
        self.assertFalse(proxy["sqrt2_is_unconditional_bound"])
        for cell in proxy["cells"].values():
            self.assertEqual(cell["observation_count"], 2)
            self.assertEqual(cell["independent_data_seed_count"], 0)
            self.assertTrue(math.isfinite(cell["terminal_tail5_sd"]))

    def test_post_hoc_tail_cannot_support_go(self) -> None:
        _, final = self.decisions()
        result = AUDIT.post_hoc_tail_audit(final)
        self.assertFalse(result["registered_in_v7_2_4"])
        self.assertFalse(result["go_direction_permitted"])
        self.assertEqual(
            result["h5_0p3_delta_min_gate_role"],
            "NOT_REGISTERED_NOT_EXECUTABLE_FOR_GO",
        )

    def test_protocol_contains_decay_and_replication_firewalls(self) -> None:
        protocol = (ROOT / "PROTOCOL_ADDENDUM.v7.3.0.md").read_text(encoding="utf-8")
        self.assertIn("Positive and negative results are both descriptive", protocol)
        self.assertIn("`INCONCLUSIVE_BUDGET_NOT_ESTIMABLE`", protocol)
        self.assertIn("independent repetition unit", protocol)
        self.assertIn("not valid for a WSD or cosine-decay tail", protocol)
        self.assertIn("It may not\nstart GPU work", protocol)

    def test_runner_is_cpu_only_and_verifies_source_manifests(self) -> None:
        runner = (ROOT / "p0" / "run_estimand_family_audit_v7_3_0.sh").read_text(
            encoding="utf-8"
        )
        self.assertNotIn("nvidia-smi", runner)
        self.assertNotIn("CUDA_VISIBLE_DEVICES", runner)
        self.assertIn("BUDGET_BRIDGE_MANIFEST.v7.2.1.sha256", runner)
        self.assertIn("BUDGET_EXTENSION_MANIFEST.v7.2.3.sha256", runner)
        self.assertIn("FINAL_BUDGET_MANIFEST.v7.2.4.sha256", runner)
        self.assertIn("C0_EXTENSION_MANIFEST.v7.1.13.sha256", runner)


if __name__ == "__main__":
    unittest.main()
