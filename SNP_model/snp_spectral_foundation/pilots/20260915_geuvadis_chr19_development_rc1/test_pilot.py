from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest

import numpy as np


MODULE_PATH = Path(__file__).with_name("run_pilot.py")
SPEC = importlib.util.spec_from_file_location("geuvadis_pilot", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
pilot = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(pilot)


class PilotUnitTests(unittest.TestCase):
    def test_group_split_has_no_group_overlap(self) -> None:
        samples = [f"S{i}" for i in range(30)]
        relationships = {sample: {"FID": f"F{i // 2}"} for i, sample in enumerate(samples)}
        train, val, test, groups = pilot.make_group_split(samples, relationships, 123)
        sets = [{groups[i] for i in idx} for idx in (train, val, test)]
        self.assertFalse(sets[0] & sets[1])
        self.assertFalse(sets[0] & sets[2])
        self.assertFalse(sets[1] & sets[2])
        self.assertEqual(len(set(train) | set(val) | set(test)), len(samples))

    def test_panel_selection_is_deterministic_unique_and_block_complete(self) -> None:
        first = pilot.choose_panel_indices(1321)
        second = pilot.choose_panel_indices(1321)
        np.testing.assert_array_equal(first, second)
        self.assertEqual(len(first), 1280)
        self.assertEqual(len(np.unique(first)), len(first))
        self.assertEqual(len(first) % 32, 0)

    def test_predictive_r2_uses_train_mean(self) -> None:
        y = np.asarray([1.0, 3.0])
        prediction = np.asarray([1.0, 3.0])
        self.assertEqual(pilot.predictive_r2(y, prediction, 0.0), 1.0)
        null = np.asarray([0.0, 0.0])
        self.assertEqual(pilot.predictive_r2(y, null, 0.0), 0.0)

    def test_bootstrap_is_reproducible_and_paired(self) -> None:
        y = np.asarray([0.0, 1.0, 2.0, 3.0])
        predictions = {arm: y + offset for arm, offset in zip("ABCDE", [1, 0.5, 0, 0.25, 0.75], strict=True)}
        first = pilot.bootstrap_statistics(y, predictions, -1.0, replicates=100, seed=7)
        second = pilot.bootstrap_statistics(y, predictions, -1.0, replicates=100, seed=7)
        self.assertEqual(first, second)
        self.assertGreater(first[1]["C-B"]["delta_r2"], 0)


if __name__ == "__main__":
    unittest.main()
