from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "validate_manifests.py"
SPEC = importlib.util.spec_from_file_location("validate_manifests", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class ValidatorTests(unittest.TestCase):
    def test_family_component_crossing_split_fails(self) -> None:
        rows = [
            {"sample_id": "a", "component_id": "c1", "component_size": "2", "split": "train"},
            {"sample_id": "b", "component_id": "c1", "component_size": "2", "split": "test"},
        ]
        with self.assertRaises(ValueError):
            MODULE.validate_family(rows)

    def test_valid_family_manifest_passes(self) -> None:
        rows = [
            {"sample_id": "a", "component_id": "c1", "component_size": "2", "split": "train"},
            {"sample_id": "b", "component_id": "c1", "component_size": "2", "split": "train"},
            {"sample_id": "c", "component_id": "c2", "component_size": "1", "split": "test"},
        ]
        result = MODULE.validate_family(rows)
        self.assertEqual(result["component_overlap_count"], 0)
        self.assertEqual(result["max_component_size"], 2)

    def test_cross_split_guard_overlap_fails(self) -> None:
        def row(block_id: str, start0: int, end0: int, split: str, gs: int, ge: int):
            return {
                "block_id": block_id,
                "chrom": "22",
                "core_start0": str(start0),
                "core_end0": str(end0),
                "core_start1": str(start0 + 1),
                "core_end1": str(end0),
                "guard_start0": str(gs),
                "guard_end0": str(ge),
                "guard_start1": str(gs + 1),
                "guard_end1": str(ge),
                "split": split,
            }

        rows = [
            row("train", 0, 10, "train", 0, 10),
            row("test", 10, 20, "test", 5, 25),
        ]
        with self.assertRaises(ValueError):
            MODULE.validate_blocks(rows)


if __name__ == "__main__":
    unittest.main()

