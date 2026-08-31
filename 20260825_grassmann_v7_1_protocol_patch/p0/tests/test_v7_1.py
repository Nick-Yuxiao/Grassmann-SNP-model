from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path
from unittest.mock import mock_open, patch

import torch


P0 = Path(__file__).resolve().parents[1]


def load(name: str):
    spec = importlib.util.spec_from_file_location(name, P0 / f"{name}.py")
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


class V710Tests(unittest.TestCase):
    def test_local_window_shape_and_parameter_parity(self) -> None:
        profiler = load("profile_models_v7_1")
        cfg = profiler.Architecture()
        self.assertEqual(cfg.attention_window, 256)
        self.assertEqual(cfg.query_block, 128)
        tokens = torch.randint(0, 4, (1, 259))
        pcs = torch.randn(1, cfg.pc_dim)
        counts = []
        for kind in profiler.MODEL_KINDS:
            model = profiler.MaskedGenotypeModel(kind, cfg)
            logits = model(tokens, pcs)
            self.assertEqual(tuple(logits.shape), (1, 259, 3))
            logits.mean().backward()
            self.assertTrue(all(
                parameter.grad is None or torch.isfinite(parameter.grad).all()
                for parameter in model.parameters()
            ))
            counts.append(profiler.parameter_count(model))
        self.assertLess(max(counts) / min(counts), 1.01)
        self.assertTrue(all(7_800_000 <= count <= 8_300_000 for count in counts))

    def test_id_list_rejects_duplicates(self) -> None:
        panel = load("build_panel_manifest_v7_1")
        path = Path("ids.tsv")
        with patch("pathlib.Path.open", mock_open(read_data="sample_id\nA\nA\n")):
            with self.assertRaises(ValueError):
                panel.read_first_column(path, {"sample_id"})

    def test_cpu_profile_is_never_valid_t03(self) -> None:
        profiler = load("profile_models_v7_1")
        cfg = profiler.Architecture()
        cfg.validate()
        self.assertEqual(cfg.attention_window, 256)


if __name__ == "__main__":
    unittest.main()
