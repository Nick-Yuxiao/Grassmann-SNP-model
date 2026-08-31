from __future__ import annotations

import importlib.util
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
MODULE_PATH = ROOT / "p0" / "build_t04_contract_v7_1_5.py"
SPEC = importlib.util.spec_from_file_location("t04_contract", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class T04ContractTests(unittest.TestCase):
    def test_constants(self) -> None:
        self.assertEqual(len(MODULE.MODELS), 3)
        self.assertEqual(len(MODULE.EXPECTED_PROFILE_SHA256), 64)
        self.assertEqual(len(MODULE.EXPECTED_VARIANT_IDS_SHA256), 64)

    def test_pilot_grid(self) -> None:
        grid = json.loads((ROOT / "PILOT_GRID.v7.1.5.yaml").read_text(encoding="utf-8"))
        self.assertEqual(grid["sequence_length"], 154850)
        self.assertEqual(grid["train_steps_per_run"], 2000)
        self.assertEqual(len(grid["mask_cells"]), 4)
        self.assertEqual(len(grid["pilot_data_seeds"]), 2)
        self.assertEqual(grid["hgdp_access"], "FORBIDDEN")
        self.assertFalse(grid["preliminary_signal_rule"]["confirmatory_claim_allowed"])


if __name__ == "__main__":
    unittest.main()
