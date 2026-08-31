from __future__ import annotations

import importlib.util
import hashlib
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
MODULE_PATH = ROOT / "p0" / "build_materialization_audit_v7_1_4.py"
SPEC = importlib.util.spec_from_file_location("materialization_audit", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class MaterializationAuditTests(unittest.TestCase):
    def test_expected_contract(self) -> None:
        self.assertEqual(MODULE.EXPECTED["joint_release"], (154850, 3264))
        self.assertEqual(MODULE.EXPECTED["donor_train"], (154850, 2247))
        self.assertEqual(MODULE.EXPECTED["hgdp_primary"], (154850, 768))
        self.assertEqual(len(MODULE.SOURCE_SHA256), 64)

    def test_sha256(self) -> None:
        path = ROOT / "PANEL_SPEC.v7.1.4.yaml"
        expected = hashlib.sha256(path.read_bytes()).hexdigest()
        self.assertEqual(MODULE.sha256(path), expected)

    def test_panel_spec(self) -> None:
        panel = json.loads((ROOT / "PANEL_SPEC.v7.1.4.yaml").read_text(encoding="utf-8"))
        self.assertEqual(panel["expected_L"], 154850)
        self.assertEqual(panel["variant_match"], "CHROM_POS_REF_ALT_EXACT")
        self.assertEqual(panel["info_policy"]["retain"], [])
        self.assertFalse(panel["gpu_required"])


if __name__ == "__main__":
    unittest.main()
