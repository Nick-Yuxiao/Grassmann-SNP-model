from __future__ import annotations

import csv
import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
FIREWALL_SCRIPT = ROOT / "p0" / "freeze_lr_validation_firewall_v7_2_0.py"
SPEC = importlib.util.spec_from_file_location("firewall_v720", FIREWALL_SCRIPT)
FIREWALL = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(FIREWALL)
SELECTOR_SCRIPT = ROOT / "p0" / "select_shared_lr_v7_2_0.py"
SELECTOR_SPEC = importlib.util.spec_from_file_location("selector_v720", SELECTOR_SCRIPT)
SELECTOR = importlib.util.module_from_spec(SELECTOR_SPEC)
assert SELECTOR_SPEC and SELECTOR_SPEC.loader
SELECTOR_SPEC.loader.exec_module(SELECTOR)


class TestV720(unittest.TestCase):
    def test_schedule_is_complete_factorial(self) -> None:
        with (ROOT / "LR_PILOT_SCHEDULE.v7.2.0.tsv").open(encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 24)
        cells = {(row["peak_lr"], row["model"], row["mask"]) for row in rows}
        self.assertEqual(len(cells), 24)
        self.assertEqual({int(row["preferred_physical_gpu"]) for row in rows}, {1, 3, 4, 5, 6, 7})
        self.assertNotIn(0, {int(row["preferred_physical_gpu"]) for row in rows})
        self.assertNotIn(2, {int(row["preferred_physical_gpu"]) for row in rows})
        for gpu in {1, 3, 4, 5, 6, 7}:
            block = [row for row in rows if int(row["preferred_physical_gpu"]) == gpu]
            self.assertEqual(len(block), 4)
            self.assertEqual(len({(row["model"], row["mask"]) for row in block}), 1)
            self.assertEqual({row["peak_lr"] for row in block}, {"0.0001", "0.0002", "0.0004", "0.0008"})

    def test_historical_ranking_is_deterministic(self) -> None:
        ids = [f"sample_{index:03d}" for index in range(249)]
        first = sorted(range(249), key=lambda index: FIREWALL.stable_int("validation_sample", 91001, ids[index]))[:32]
        second = sorted(range(249), key=lambda index: FIREWALL.stable_int("validation_sample", 91001, ids[index]))[:32]
        self.assertEqual(first, second)
        self.assertEqual(len(set(first)), 32)

    def test_protocol_keeps_architecture_decision_blocked(self) -> None:
        protocol = (ROOT / "PROTOCOL.v7.2.0.md").read_text(encoding="utf-8")
        self.assertIn("Formal A1-R remains blocked", protocol)
        self.assertIn("may not support GO/NO-GO", protocol)

    def test_selector_uses_lowest_near_best_lr(self) -> None:
        summaries = {
            "0.0001": {"eligible": True, "mean_gain_vs_1e_4": 0.0},
            "0.0002": {"eligible": True, "mean_gain_vs_1e_4": 0.0021},
            "0.0004": {"eligible": True, "mean_gain_vs_1e_4": 0.0025},
            "0.0008": {"eligible": False, "mean_gain_vs_1e_4": 0.0030},
        }
        selected, status = SELECTOR.choose_shared_lr(summaries)
        self.assertEqual(selected, 0.0002)
        self.assertEqual(status, "LR_PILOT_SHARED_PEAK_SELECTED")


if __name__ == "__main__":
    unittest.main()
