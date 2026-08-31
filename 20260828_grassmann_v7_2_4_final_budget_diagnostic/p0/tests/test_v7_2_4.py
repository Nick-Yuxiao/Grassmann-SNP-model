from __future__ import annotations

import csv
import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ASSESS_PATH = ROOT / "p0" / "assess_final_budget_v7_2_4.py"
SPEC = importlib.util.spec_from_file_location("assess_v724", ASSESS_PATH)
ASSESS = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(ASSESS)


class TestV724(unittest.TestCase):
    def schedule(self) -> list[dict[str, str]]:
        with (ROOT / "FINAL_BUDGET_SCHEDULE.v7.2.4.tsv").open(
            encoding="utf-8-sig", newline=""
        ) as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def test_schedule_is_complete_selected_lr_factorial(self) -> None:
        rows = self.schedule()
        self.assertEqual(len(rows), 6)
        self.assertEqual(len({row["source_run_id"] for row in rows}), 6)
        self.assertEqual(
            {(row["model"], row["mask"]) for row in rows},
            {(model, mask) for model in ASSESS.MODELS for mask in ASSESS.MASKS},
        )
        self.assertEqual({row["learning_rate"] for row in rows}, {"0.0004"})

    def test_schedule_is_exact_30k_to_40k_one_cell_per_gpu(self) -> None:
        rows = self.schedule()
        self.assertEqual({row["source_step"] for row in rows}, {"30000"})
        self.assertEqual({row["target_step"] for row in rows}, {"40000"})
        self.assertEqual(
            {int(row["preferred_physical_gpu"]) for row in rows}, {1, 3, 4, 5, 6, 7}
        )

    def test_terminal_stable(self) -> None:
        values = [(step, 1.0 - step / 10_000_000) for step in range(30250, 40001, 250)]
        train = [(step, 1.1 - step / 10_000_000) for step in range(30250, 40001, 250)]
        result = ASSESS.terminal_summary(values, train)
        self.assertTrue(result["all_absolute_changes_le_0p002"])
        self.assertEqual(result["shape_class"], "STABLE")
        self.assertFalse(result["instability"])

    def test_terminal_stable_but_accelerating(self) -> None:
        values = []
        for step in range(30250, 40001, 250):
            if step <= 34000:
                value = 1.0
            elif step <= 36000:
                value = 0.9996
            elif step <= 38000:
                value = 0.9990
            else:
                value = 0.9980
            values.append((step, value))
        train = [(step, 1.1 - step / 10_000_000) for step in range(30250, 40001, 250)]
        result = ASSESS.terminal_summary(values, train)
        self.assertTrue(result["all_absolute_changes_le_0p002"])
        self.assertEqual(result["shape_class"], "STABLE_BUT_ACCELERATING")

    def test_decision_branches_are_bounded(self) -> None:
        self.assertEqual(ASSESS.decision_branch(False, False, True, ["cell"])[3], 7)
        self.assertEqual(ASSESS.decision_branch(False, False, False, [])[3], 6)
        self.assertEqual(ASSESS.decision_branch(False, True, True, [])[3], 4)
        adequate = ASSESS.decision_branch(False, False, True, [])
        self.assertEqual(adequate[0], "FINAL_BUDGET_40K_ADEQUATE")
        self.assertEqual(adequate[1], 50000)
        self.assertEqual(adequate[3], 0)

    def test_protocol_forbids_selective_and_unbounded_extension(self) -> None:
        text = (ROOT / "PROTOCOL_ADDENDUM.v7.2.4.md").read_text(encoding="utf-8")
        self.assertIn("complete `3 models x 2 masks` factorial", text)
        self.assertIn("Continuing only the two failed cells is forbidden", text)
        self.assertIn("No run may continue beyond 40k", text)
        self.assertIn("cannot produce an architecture GO/NO-GO decision", text)

    def test_trainer_restores_optimizer_and_rng(self) -> None:
        text = (ROOT / "p0" / "train_final_budget_v7_2_4.py").read_text(
            encoding="utf-8"
        )
        self.assertIn("optimizer.load_state_dict", text)
        self.assertIn("torch.set_rng_state", text)
        self.assertIn("torch.cuda.set_rng_state_all", text)


if __name__ == "__main__":
    unittest.main()
