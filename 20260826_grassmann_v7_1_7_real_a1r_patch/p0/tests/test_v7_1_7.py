from __future__ import annotations

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader
    spec.loader.exec_module(module)
    return module


class V717Tests(unittest.TestCase):
    def test_grid_counts_and_policy(self):
        grid = json.loads((ROOT / "A1R_GRID.v7.1.7.yaml").read_text(encoding="utf-8"))
        self.assertEqual(grid["primary_A1R_100pct"]["runs"], 120)
        self.assertEqual(grid["sample_size_diagnostic"]["added_runs"], 24)
        self.assertEqual(grid["data_contract"]["hgdp_access"], "FORBIDDEN")
        self.assertNotIn(0, grid["allowed_physical_gpus"])

    def test_schedule(self):
        module = load("contract", ROOT / "p0" / "build_a1r_contract_v7_1_7.py")
        grid = json.loads((ROOT / "A1R_GRID.v7.1.7.yaml").read_text(encoding="utf-8"))
        rows = module.build_schedule(grid)
        self.assertEqual(len(rows), 144)
        self.assertEqual(sum(r["stage"] == "PRIMARY_100P" for r in rows), 120)
        self.assertEqual(sum(r["stage"] == "DIAGNOSTIC_SIZE" for r in rows), 24)
        self.assertTrue(all(r["preferred_physical_gpu"] in range(1, 7) for r in rows))
        for curve in {r["size_curve_id"] for r in rows if r["stage"] == "DIAGNOSTIC_SIZE"}:
            self.assertEqual(len({r["preferred_physical_gpu"] for r in rows if r["size_curve_id"] == curve}), 1)

    def test_nested_allocation(self):
        module = load("subsets", ROOT / "p0" / "freeze_a1r_subsets_v7_1_7.py")
        counts = {"A": 101, "B": 87, "C": 59}
        q50 = module.allocate(counts, 124)
        q25 = module.allocate(q50, 62)
        self.assertEqual(sum(q25.values()), 62)
        self.assertEqual(sum(q50.values()), 124)
        self.assertTrue(all(q25[p] <= q50[p] <= counts[p] for p in counts))


if __name__ == "__main__":
    unittest.main()
