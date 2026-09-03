from __future__ import annotations

import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
P0 = ROOT / "p0"
sys.path.insert(0, str(P0))

import build_architecture_pilot_readiness_v7_7_3 as builder


def valid_source() -> dict:
    return {
        "status": "LONG_RANGE_TASK_VALIDITY_PASS",
        "protocol_version": "v7.7.2",
        "next_authorized_stage": "DRAFT_V7_7_3_LONG_RANGE_ARCHITECTURE_PILOT_NO_GPU",
        "gates": {
            "local_negative_pass": True,
            "oracle_positive_pass": True,
            "conventional_global_positive_pass": True,
            "target_shuffled_negative_pass": True,
            "all_pass": True,
        },
        "authorization": {
            "gpu_used": False,
            "grassmann_fitted": False,
            "architecture_decision_permitted": False,
        },
    }


class TestV773(unittest.TestCase):
    def test_source_pass_contract_accepted(self) -> None:
        builder.validate_source(valid_source())

    def test_source_rejects_non_pass(self) -> None:
        source = valid_source()
        source["status"] = "LONG_RANGE_TASK_INVALID"
        source["gates"]["all_pass"] = False
        with self.assertRaises(ValueError):
            builder.validate_source(source)

    def test_source_rejects_gpu_use(self) -> None:
        source = valid_source()
        source["authorization"]["gpu_used"] = True
        with self.assertRaises(ValueError):
            builder.validate_source(source)

    def test_source_rejects_grassmann_fitted(self) -> None:
        source = valid_source()
        source["authorization"]["grassmann_fitted"] = True
        with self.assertRaises(ValueError):
            builder.validate_source(source)

    def test_arm_map_is_complete_2x2_with_controls(self) -> None:
        arms = builder.load_arms(ROOT / "ARCHITECTURE_PILOT_ARM_MAP.v7.7.3.tsv")
        checks = builder.static_contract_checks(arms)
        self.assertTrue(checks["all_pass"])
        self.assertEqual(checks["decision_eligible_cell_count"], 4)
        self.assertTrue(checks["full_factorial_complete"])

    def test_router_present_arms_share_full_sequence_inputs(self) -> None:
        arms = builder.load_arms(ROOT / "ARCHITECTURE_PILOT_ARM_MAP.v7.7.3.tsv")
        checks = builder.static_contract_checks(arms)
        self.assertTrue(checks["router_present_inputs_identical"])
        self.assertTrue(checks["router_present_not_source_positions"])

    def test_no_arm_execution_authorized(self) -> None:
        arms = builder.load_arms(ROOT / "ARCHITECTURE_PILOT_ARM_MAP.v7.7.3.tsv")
        self.assertTrue(all(row["execution_authorized"] == "FALSE" for row in arms))

    def test_spec_remains_non_gpu_and_non_architectural(self) -> None:
        spec = (ROOT / "LONG_RANGE_ARCHITECTURE_PILOT_SPEC.v7.7.3.yaml").read_text(encoding="utf-8")
        self.assertIn("gpu_authorized: false", spec)
        self.assertIn("grassmann_training_authorized: false", spec)
        self.assertIn("architecture_decision_permitted: false", spec)
        self.assertIn("NLL_LR_A01_MINUS_NLL_LR_A11", spec)
        self.assertIn("IMPLEMENT_V7_7_4_CPU_BLINDED_VARIANCE_PILOT_NO_GPU", spec)


if __name__ == "__main__":
    unittest.main()
