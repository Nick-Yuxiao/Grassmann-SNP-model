from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

import torch


P0 = Path(__file__).resolve().parents[1]
SERVER_OPS = P0.parent / "server_ops"


def load(name: str):
    spec = importlib.util.spec_from_file_location(name, P0 / f"{name}.py")
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


class P0Tests(unittest.TestCase):
    def test_model_shapes_and_parameter_parity(self) -> None:
        profiler = load("profile_models")
        cfg = profiler.Architecture()
        counts = []
        tokens = torch.randint(0, 4, (1, 130))
        pcs = torch.randn(1, cfg.pc_dim)
        for kind in ("local_attn", "local_attn_gpc", "grassmann_full"):
            model = profiler.MaskedGenotypeModel(kind, cfg)
            self.assertEqual(model(tokens, pcs).shape, (1, 130, 3))
            counts.append(profiler.parameter_count(model))
        self.assertLess(max(counts) / min(counts), 1.01)
        self.assertTrue(all(7_900_000 <= count <= 8_200_000 for count in counts))

    def test_audit_csv_parser(self) -> None:
        audit = load("audit_server")
        rows = audit.csv_rows("0, GPU-a, RTX 5090, 32768\n", ["index", "uuid", "name", "memory"])
        self.assertEqual(rows, [{"index": "0", "uuid": "GPU-a", "name": "RTX 5090", "memory": "32768"}])

    def test_frozen_constant_names_are_independent(self) -> None:
        root = P0.parent
        decisions = (root / "DECISIONS.v7.0.1.tsv").read_text(encoding="utf-8")
        metrics = (root / "METRIC_DEFINITIONS.md").read_text(encoding="utf-8")
        for name in ("delta_min", "delta_NI", "delta_LD", "overfit_thr", "pc_control_thr"):
            self.assertIn(name, decisions)
            self.assertIn(name, metrics)

    def test_gpu_version_parser(self) -> None:
        path = SERVER_OPS / "gpu_test_5090.py"
        spec = importlib.util.spec_from_file_location("gpu_test_5090", path)
        assert spec and spec.loader
        module = importlib.util.module_from_spec(spec)
        sys.modules["gpu_test_5090"] = module
        spec.loader.exec_module(module)
        self.assertEqual(module.pair("2.7.1+cu128"), (2, 7))
        self.assertEqual(module.pair("12.8"), (12, 8))
        self.assertIsNone(module.pair(None))


if __name__ == "__main__":
    unittest.main()
