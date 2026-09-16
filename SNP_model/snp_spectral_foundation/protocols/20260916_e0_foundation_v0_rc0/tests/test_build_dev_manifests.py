from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "build_dev_manifests.py"
SPEC = importlib.util.spec_from_file_location("build_dev_manifests", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class ManifestTests(unittest.TestCase):
    def test_related_samples_share_component_and_split(self) -> None:
        samples = {
            name: MODULE.Sample(name, "P1", "S1", "unknown")
            for name in ["a", "b", "c", "d", "e", "f"]
        }
        components = MODULE.build_components(samples, [("a", "b"), ("b", "c")])
        member_to_component = {
            member: component for component, members in components.items() for member in members
        }
        self.assertEqual(member_to_component["a"], member_to_component["b"])
        self.assertEqual(member_to_component["b"], member_to_component["c"])
        assignments = MODULE.assign_component_splits(components, samples, seed=7)
        self.assertEqual(
            assignments[member_to_component["a"]], assignments[member_to_component["c"]]
        )

    def test_component_assignment_is_deterministic(self) -> None:
        samples = {
            f"s{i}": MODULE.Sample(f"s{i}", "P1", "S1", "unknown") for i in range(30)
        }
        components = MODULE.build_components(samples, [])
        first = MODULE.assign_component_splits(components, samples, seed=11)
        second = MODULE.assign_component_splits(components, samples, seed=11)
        self.assertEqual(first, second)
        counts = {split: list(first.values()).count(split) for split in MODULE.SPLIT_RATIOS}
        self.assertEqual(sum(counts.values()), 30)
        for split, ratio in MODULE.SPLIT_RATIOS.items():
            self.assertLessEqual(abs(counts[split] - 30 * ratio), 1.0)

    def test_coordinate_contract_and_guard_isolation(self) -> None:
        blocks = MODULE.fixed_bp_blocks("22", contig_length=50_000_000, block_bp=5_000_000)
        rows = MODULE.assign_blocks_with_guards(blocks, guard_bp=1_000_000, seed=19)
        MODULE.validate_block_rows(rows)
        for row in rows:
            self.assertEqual(row["core_start1"], row["core_start0"] + 1)
            self.assertEqual(row["core_end1"], row["core_end0"])
        self.assertTrue(any(row["split"] == "test" for row in rows))
        self.assertTrue(any(row["split"] == "validation" for row in rows))
        self.assertTrue(any(str(row["split"]).startswith("buffer_") for row in rows))

    def test_invalid_or_overlapping_blocks_fail(self) -> None:
        with self.assertRaises(ValueError):
            MODULE.CoreBlock("bad", "22", 5, 5)
        blocks = [
            MODULE.CoreBlock("a", "22", 0, 10),
            MODULE.CoreBlock("b", "22", 9, 20),
        ]
        with self.assertRaises(ValueError):
            MODULE.assign_blocks_with_guards(blocks, guard_bp=1, seed=1)

    def test_cm_blocks_and_guards_obey_coordinate_contract(self) -> None:
        positions1 = [101, 201, 301, 401, 501, 601]
        map_cm = [0.0, 0.4, 0.8, 1.2, 1.6, 2.0]
        blocks = MODULE.fixed_cm_blocks("22", positions1, map_cm, block_cm=0.5)
        guards = MODULE.cm_guard_bounds(blocks, positions1, map_cm, guard_cm=0.5)
        self.assertGreaterEqual(len(blocks), 3)
        for block in blocks:
            guard_start0, guard_end0 = guards[block.block_id]
            self.assertLessEqual(guard_start0, block.start0)
            self.assertGreaterEqual(guard_end0, block.end0)
        self.assertEqual(blocks[0].start0, positions1[0] - 1)
        self.assertEqual(blocks[-1].end0, positions1[-1])


if __name__ == "__main__":
    unittest.main()
