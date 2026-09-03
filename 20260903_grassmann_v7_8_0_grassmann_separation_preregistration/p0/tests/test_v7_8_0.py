from __future__ import annotations

import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
P0 = ROOT / "p0"
sys.path.insert(0, str(P0))

import build_separation_prereg_readiness_v7_8_0 as builder


def valid_source() -> dict:
    return {
        "status": "LONG_RANGE_DIFFICULTY_REPLAN_CONTRACT_SIGNED_NO_LAUNCH",
        "protocol_version": "v7.7.7",
        "next_authorized_stage": "IMPLEMENT_V7_7_8_LONG_RANGE_DIFFICULTY_REPLAN_HARNESS_NO_LAUNCH",
        "authorization": {"gpu_authorized": False, "grassmann_fitted": False},
    }


class TestV780(unittest.TestCase):
    def test_source_replan_contract_accepted(self) -> None:
        builder.validate_source(valid_source())

    def test_source_rejects_wrong_status(self) -> None:
        source = valid_source()
        source["status"] = "SOMETHING_ELSE"
        with self.assertRaises(ValueError):
            builder.validate_source(source)

    def test_source_rejects_gpu_or_grassmann(self) -> None:
        source = valid_source()
        source["authorization"]["gpu_authorized"] = True
        with self.assertRaises(ValueError):
            builder.validate_source(source)

    def test_arm_map_has_fair_suite_grassmann_and_controls(self) -> None:
        arms = builder.load_arms(ROOT / "SEPARATION_ARM_MAP.v7.8.0.tsv")
        checks = builder.static_arm_checks(arms)
        self.assertTrue(checks["all_pass"])
        self.assertGreaterEqual(checks["conventional_family_count"], 3)
        self.assertTrue(checks["has_grassmann_arm"])
        self.assertTrue(checks["controls_present"])

    def test_no_arm_execution_authorized(self) -> None:
        arms = builder.load_arms(ROOT / "SEPARATION_ARM_MAP.v7.8.0.tsv")
        self.assertTrue(all(row["execution_authorized"] == "FALSE" for row in arms))

    def test_spec_is_gpu_gated_not_authorized(self) -> None:
        spec = (ROOT / "GRASSMANN_SEPARATION_PREREG_SPEC.v7.8.0.yaml").read_text(encoding="utf-8")
        self.assertIn("gpu_authorized: false", spec)
        self.assertIn("gpu_gated: true", spec)
        self.assertIn("lr1_required: false", spec)
        self.assertIn("IMPLEMENT_V7_8_1_CPU_GRASSMANN_SEPARATION_HARNESS_NO_LAUNCH", spec)

    def test_spec_guards_against_strawman_and_taskshopping(self) -> None:
        spec = (ROOT / "GRASSMANN_SEPARATION_PREREG_SPEC.v7.8.0.yaml").read_text(encoding="utf-8")
        self.assertIn("not_a_strawman: true", spec)
        self.assertIn("reverse_engineered_from_grassmann: false", spec)
        self.assertIn("conventional_convergence_compute_sufficiency_audit_required: true", spec)

    def test_protocol_states_separation_and_route_closure(self) -> None:
        protocol = (ROOT / "PROTOCOL_ADDENDUM.v7.8.0.md").read_text(encoding="utf-8")
        self.assertIn("separation", protocol.lower())
        self.assertIn("A1-R", protocol)
        self.assertIn("close the v7", protocol.lower())


if __name__ == "__main__":
    unittest.main()
