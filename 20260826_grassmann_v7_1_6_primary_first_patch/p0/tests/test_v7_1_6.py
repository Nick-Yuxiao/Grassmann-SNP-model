from __future__ import annotations

import importlib.util
import json
import unittest
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
MODULE_PATH = ROOT / "p0" / "build_t04_primary_first_v7_1_6.py"
SPEC = importlib.util.spec_from_file_location("primary_first", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class PrimaryFirstTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.grid = json.loads((ROOT / "PRIMARY_FIRST_GRID.v7.1.6.yaml").read_text(encoding="utf-8"))
        cls.rows = MODULE.make_primary_schedule(cls.grid)

    def test_primary_run_count_and_models(self) -> None:
        self.assertEqual(len(self.rows), 120)
        self.assertEqual(Counter(row["model"] for row in self.rows), Counter({model: 40 for model in MODULE.MODELS}))

    def test_blocks_keep_three_models_on_one_gpu(self) -> None:
        blocks = defaultdict(list)
        for row in self.rows:
            blocks[row["block_id"]].append(row)
        self.assertEqual(len(blocks), 40)
        for rows in blocks.values():
            self.assertEqual({row["model"] for row in rows}, set(MODULE.MODELS))
            self.assertEqual(len({row["preferred_physical_gpu"] for row in rows}), 1)

    def test_gpu_balance_and_gpu0_forbidden(self) -> None:
        block_first = [row for row in self.rows if row["order_in_block"] == 1]
        counts = Counter(row["preferred_physical_gpu"] for row in block_first)
        self.assertNotIn(0, counts)
        self.assertEqual(set(counts), {1, 2, 3, 4, 5, 6})
        self.assertLessEqual(max(counts.values()) - min(counts.values()), 1)

    def test_convergence_rule(self) -> None:
        convergence = self.grid["convergence_pilot"]
        self.assertEqual(convergence["candidate_common_steps"], [4000, 6000, 8000, 10000])
        self.assertEqual(convergence["hgdp_access"], "FORBIDDEN")
        self.assertEqual(convergence["budget_selection_blinding"], "WITHIN_CURVE_ONLY_NO_BETWEEN_MODEL_DELTAS")


if __name__ == "__main__":
    unittest.main()
