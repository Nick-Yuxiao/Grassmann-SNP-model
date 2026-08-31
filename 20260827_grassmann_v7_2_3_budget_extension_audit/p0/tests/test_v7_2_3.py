from __future__ import annotations

import csv
import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ASSESS_PATH = ROOT / "p0" / "assess_budget_extension_v7_2_3.py"
SPEC = importlib.util.spec_from_file_location("assess_v723", ASSESS_PATH)
ASSESS = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(ASSESS)


class TestV723(unittest.TestCase):
    def test_schedule_extends_every_source_once(self) -> None:
        with (ROOT / "BUDGET_EXTENSION_SCHEDULE.v7.2.3.tsv").open(
            encoding="utf-8-sig", newline=""
        ) as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 12)
        self.assertEqual(len({row["source_run_id"] for row in rows}), 12)
        self.assertEqual({row["source_step"] for row in rows}, {"20000"})
        self.assertEqual({row["target_step"] for row in rows}, {"30000"})

    def test_schedule_preserves_paired_gpu_blocks(self) -> None:
        with (ROOT / "BUDGET_EXTENSION_SCHEDULE.v7.2.3.tsv").open(
            encoding="utf-8-sig", newline=""
        ) as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        for gpu in {1, 3, 4, 5, 6, 7}:
            block = [row for row in rows if int(row["preferred_physical_gpu"]) == gpu]
            self.assertEqual(len(block), 2)
            self.assertEqual({row["learning_rate"] for row in block}, {"0.0001", "0.0004"})
            self.assertEqual(len({(row["model"], row["mask"]) for row in block}), 1)

    def test_isotonic_is_nonincreasing(self) -> None:
        fitted = ASSESS.isotonic_nonincreasing([1.0, 0.8, 0.9, 0.7])
        self.assertTrue(all(left >= right for left, right in zip(fitted, fitted[1:])))
        self.assertAlmostEqual(fitted[1], 0.85)
        self.assertAlmostEqual(fitted[2], 0.85)

    def test_acceleration_on_known_curves(self) -> None:
        control = [(step, 1.0 - step / 60000) for step in range(250, 30001, 250)]
        selected = [(step, 1.0 - step / 30000) for step in range(250, 30001, 250)]
        result = ASSESS.acceleration(control, selected)
        self.assertGreater(float(result["median_acceleration_ratio"]), 1.8)
        self.assertLess(float(result["median_acceleration_ratio"]), 2.2)
        self.assertGreater(float(result["implied_selected_step_reduction_fraction"]), 0.45)
        self.assertEqual(result["censoring"], "NONE_BY_CONSTRUCTION_COMMON_SUPPORT_ONLY")

    def test_stable_but_accelerating_is_flagged_without_redefining_primary_gate(self) -> None:
        values = []
        for step in range(250, 30001, 250):
            if step <= 24000:
                value = 1.0
            elif step <= 26000:
                value = 0.9995
            elif step <= 28000:
                value = 0.99875
            else:
                value = 0.99775
            values.append((step, value))
        train = [(step, 1.0 - step / 1000000) for step in range(20250, 30001, 250)]
        result = ASSESS.terminal_summary(values, train)
        self.assertTrue(result["all_absolute_changes_le_0p002"])
        self.assertEqual(result["shape_class"], "STABLE_BUT_ACCELERATING")
        self.assertTrue(result["shape_flag_is_non_primary"])

    def test_decision_branch_distinguishes_primary_failure_and_shape_review(self) -> None:
        self.assertEqual(ASSESS.decision_branch(False, False, True, ["cell"])[3], 7)
        self.assertEqual(ASSESS.decision_branch(False, False, False, ["cell"])[3], 6)
        self.assertEqual(ASSESS.decision_branch(False, False, True, [])[3], 0)
        self.assertEqual(ASSESS.decision_branch(False, True, True, [])[3], 4)

    def test_protocol_blocks_selective_extension_and_architecture_decision(self) -> None:
        text = (ROOT / "PROTOCOL_ADDENDUM.v7.2.3.md").read_text(encoding="utf-8")
        self.assertIn("Selective continuation", text)
        self.assertIn("GO/NO-GO decisions are forbidden", text)
        self.assertIn("No further automatic extension", text)
        self.assertIn("does not retroactively redefine", text)
        self.assertIn("not required to fail", text)
        self.assertIn("1-1/R", text)


if __name__ == "__main__":
    unittest.main()
