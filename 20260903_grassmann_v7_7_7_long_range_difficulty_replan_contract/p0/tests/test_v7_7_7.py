from __future__ import annotations

import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
P0 = ROOT / "p0"
sys.path.insert(0, str(P0))

import build_difficulty_replan_readiness_v7_7_7 as builder


def valid_source() -> dict:
    return {
        "status": "LONG_RANGE_TASK_DIFFICULTY_UNRESOLVED",
        "selected_k": None,
        "protocol_version": "v7.7.6",
        "next_authorized_stage": "STOP_NO_AUTOMATIC_TASK_EXPANSION_OR_GPU",
        "authorization": {"gpu_used": False, "grassmann_fitted": False},
    }


class TestV777(unittest.TestCase):
    def test_source_unresolved_accepted(self) -> None:
        builder.validate_source(valid_source())

    def test_source_rejects_selected_k(self) -> None:
        source = valid_source()
        source["selected_k"] = 6
        with self.assertRaises(ValueError):
            builder.validate_source(source)

    def test_source_rejects_non_unresolved_status(self) -> None:
        source = valid_source()
        source["status"] = "LONG_RANGE_TASK_DIFFICULTY_SELECTED"
        with self.assertRaises(ValueError):
            builder.validate_source(source)

    def test_source_rejects_gpu_or_grassmann(self) -> None:
        source = valid_source()
        source["authorization"]["gpu_used"] = True
        with self.assertRaises(ValueError):
            builder.validate_source(source)

    def test_search_grid_is_predeclared_and_nonexecutable(self) -> None:
        rows = builder.search_grid_rows()
        # 4 label families x 1 geometry x 3 distractor levels x 2 budgets = 24
        self.assertEqual(len(rows), 24)
        labels = {r[0] for r in rows}
        self.assertEqual(labels, {"parity", "majority_threshold", "weighted_threshold", "noisy_threshold"})
        self.assertTrue(all(r[1] == "random_positions_in_token_marker" for r in rows))

    def test_spec_forbids_handing_positions_and_gpu(self) -> None:
        spec = (ROOT / "LONG_RANGE_DIFFICULTY_REPLAN_SPEC.v7.7.7.yaml").read_text(encoding="utf-8")
        self.assertIn("baseline_handed_source_positions: false", spec)
        self.assertIn("gpu_authorized: false", spec)
        self.assertIn("grassmann_consulted_in_selection: false", spec)

    def test_spec_encodes_route_closure_branch(self) -> None:
        spec = (ROOT / "LONG_RANGE_DIFFICULTY_REPLAN_SPEC.v7.7.7.yaml").read_text(encoding="utf-8")
        self.assertIn("CLOSE_V7_GRASSMANN_PRIMARY_ROUTE_NO_TASK_EXPANSION", spec)
        self.assertIn("IMPLEMENT_V7_7_8_LONG_RANGE_DIFFICULTY_REPLAN_HARNESS_NO_LAUNCH", spec)

    def test_protocol_reports_bimodal_finding(self) -> None:
        protocol = (ROOT / "PROTOCOL_ADDENDUM.v7.7.7.md").read_text(encoding="utf-8")
        self.assertIn("bimodal", protocol)
        self.assertIn("LONG_RANGE_TASK_DIFFICULTY_UNRESOLVED", protocol)


if __name__ == "__main__":
    unittest.main()
