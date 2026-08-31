from __future__ import annotations

"""Unit tests for the V7 CPU learning-smoke.

These require torch and access to the frozen ``profile_models_v7_1.py``. Point
``V7_MODEL_DIR`` at the release ``p0/`` that holds it (the same dir the smoke is
dropped into), or run from that dir. Tests skip cleanly if either is missing so
``py_compile`` / collection never hard-fails on a torch-less box.
"""

import importlib
import os
import sys
import unittest
from pathlib import Path

P0_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(P0_DIR))

try:
    import torch  # noqa: F401

    HAVE_TORCH = True
except Exception:  # pragma: no cover
    HAVE_TORCH = False


def _load_smoke():
    return importlib.import_module("run_learning_smoke")


def _load_frozen(smoke):
    model_dir = os.environ.get("V7_MODEL_DIR")
    return smoke._import_frozen_model(model_dir)


@unittest.skipUnless(HAVE_TORCH, "torch not available")
class LearningSmokeTests(unittest.TestCase):
    def setUp(self) -> None:
        self.smoke = _load_smoke()
        try:
            self.frozen = _load_frozen(self.smoke)
        except SystemExit as exc:  # frozen model not locatable in this env
            self.skipTest(str(exc))
        self.cfg = self.frozen.Architecture()
        self.period = self.cfg.attention_window

    def test_dgp_is_deterministic_and_follows_formula(self) -> None:
        pat = self.smoke.build_pos_pattern(self.period, 70101, torch)
        kwargs = dict(
            torch=torch,
            length=512,
            batch_size=8,
            pc_dim=self.cfg.pc_dim,
            seed=123,
            mask_rate=0.15,
            pos_pattern=pat,
            period=self.period,
            device=torch.device("cpu"),
        )
        t1, p1, y1, m1 = self.smoke.structured_batch(**kwargs)
        t2, p2, y2, m2 = self.smoke.structured_batch(**kwargs)
        self.assertTrue(torch.equal(y1, y2))
        self.assertTrue(torch.equal(t1, t2))
        self.assertTrue(torch.equal(m1, m2))
        # target[b, j] == (pos_pattern[j mod period] + shift[b]) mod 3, with shift
        # recoverable from any single column.
        positions = torch.arange(512) % self.period
        base = pat[positions]
        shift = (y1[:, 0] - base[0]) % 3
        expected = (base.unsqueeze(0) + shift.unsqueeze(1)) % 3
        self.assertTrue(torch.equal(y1, expected))
        # masked positions carry the mask token (== genotype_states); unmasked equal target
        self.assertTrue(torch.all(t1[m1] == self.cfg.genotype_states))
        self.assertTrue(torch.equal(t1[~m1], y1[~m1]))

    def test_all_three_arms_forward_shape(self) -> None:
        pat = self.smoke.build_pos_pattern(self.period, 70101, torch)
        tokens, pcs, _y, _m = self.smoke.structured_batch(
            torch=torch, length=512, batch_size=4, pc_dim=self.cfg.pc_dim,
            seed=1, mask_rate=0.15, pos_pattern=pat, period=self.period,
            device=torch.device("cpu"),
        )
        for kind in self.frozen.MODEL_KINDS:
            model = self.frozen.MaskedGenotypeModel(kind, self.cfg)
            logits = model(tokens, pcs)
            self.assertEqual(tuple(logits.shape), (4, 512, self.cfg.genotype_states))

    def test_grassmann_wedge_degenerate_at_one_block(self) -> None:
        """Documents the L<=256 trap: a single 256-block => zero wedge channel."""
        mixer = self.frozen.GrassmannBlockMixer(self.cfg)
        rank = self.cfg.grassmann_rank
        # one block
        reduced1 = torch.randn(2, 1, rank)
        context1 = reduced1.mean(dim=1, keepdim=True)
        w1 = (
            reduced1[..., mixer.wedge_i] * context1[..., mixer.wedge_j]
            - reduced1[..., mixer.wedge_j] * context1[..., mixer.wedge_i]
        )
        self.assertEqual(float(w1.abs().max()), 0.0)
        # two distinct blocks
        reduced2 = torch.randn(2, 2, rank)
        context2 = reduced2.mean(dim=1, keepdim=True)
        w2 = (
            reduced2[..., mixer.wedge_i] * context2[..., mixer.wedge_j]
            - reduced2[..., mixer.wedge_j] * context2[..., mixer.wedge_i]
        )
        self.assertGreater(float(w2.abs().max()), 0.0)

    def test_cheapest_arm_learns(self) -> None:
        """End-to-end: loss drops meaningfully below its start on learnable data."""
        pat = self.smoke.build_pos_pattern(self.period, 70101, torch)
        row = self.smoke.train_one_arm(
            frozen=self.frozen, torch=torch, kind="local_attn_8m_w256",
            cfg=self.cfg, length=512, steps=80, batch_size=8, lr=1e-3,
            mask_rate=0.15, seed=70101, pos_pattern=pat, period=self.period,
            device=torch.device("cpu"), eval_every=20,
        )
        self.assertLess(row["final_eval_loss"], row["initial_eval_loss"] - 0.10)


if __name__ == "__main__":
    unittest.main()
