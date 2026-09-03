from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


builder = load("long_range_pilot_v773", ROOT / "p0/build_architecture_pilot_readiness_v7_7_3.py")


class TestV773(unittest.TestCase):
    def test_four_factorial_cells_are_declared(self) -> None:
        spec = (ROOT / "LONG_RANGE_ARCHITECTURE_PILOT_SPEC.v7.7.3.yaml").read_text(encoding="utf-8")
        for cell in ("R0_LOCAL", "R1_ROUTER", "R2_GRASSMANN_ONLY", "R3_ROUTER_GRASSMANN"):
            self.assertIn(cell, spec)

    def test_primary_contrast_uses_same_router(self) -> None:
        protocol = (ROOT / "PROTOCOL_ADDENDUM.v7.7.3.md").read_text(encoding="utf-8")
        self.assertIn("NLL(R1_ROUTER) - NLL(R3_ROUTER_GRASSMANN)", protocol)
        self.assertIn("identical", protocol)

    def test_pilot_is_not_execution_authorized(self) -> None:
        spec = (ROOT / "LONG_RANGE_ARCHITECTURE_PILOT_SPEC.v7.7.3.yaml").read_text(encoding="utf-8")
        self.assertIn("gpu_authorized: false", spec)
        self.assertIn("pilot_execution_authorized: false", spec)

    def test_truth_seeds_are_disjoint_from_validity_execution(self) -> None:
        spec = (ROOT / "LONG_RANGE_ARCHITECTURE_PILOT_SPEC.v7.7.3.yaml").read_text(encoding="utf-8")
        self.assertNotIn("77201", spec)
        self.assertIn("77301", spec)

    def test_blinding_forbids_arm_means(self) -> None:
        spec = (ROOT / "LONG_RANGE_ARCHITECTURE_PILOT_SPEC.v7.7.3.yaml").read_text(encoding="utf-8")
        self.assertIn("release_pilot_arm_means: false", spec)

    def test_builder_requires_source_pass_and_cpu_only(self) -> None:
        text = (ROOT / "p0/build_architecture_pilot_readiness_v7_7_3.py").read_text(encoding="utf-8")
        self.assertIn('!= "LONG_RANGE_TASK_VALIDITY_PASS"', text)
        self.assertIn('get("gpu_used") is not False', text)

    def test_schedule_has_48_nonexecutable_cells(self) -> None:
        rows = builder.pilot_rows()
        self.assertEqual(len(rows), 48)
        self.assertEqual({row[2] for row in rows}, {"R0_LOCAL", "R1_ROUTER", "R2_GRASSMANN_ONLY", "R3_ROUTER_GRASSMANN"})
        self.assertEqual({row[0] for row in rows}, set(range(77301, 77307)))


if __name__ == "__main__":
    unittest.main()
