from __future__ import annotations

import gzip
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
AUDITOR = ROOT / "p0" / "audit_site_stage1_v7_1_2.py"


class SiteAuditTest(unittest.TestCase):
    def test_nonempty_metrics_and_duplicate_audit(self) -> None:
        rows = [
            "chr22\t100\tA\tG\t10\t100\t0.1\t0.1\t0\n",
            "chr22\t200\tC\tT\t2\t100\t0.02\t0.02\t0.001\n",
            "chr22\t200\tC\tG\t3\t100\t0.03\t0.03\t0.01\n",
        ]
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parent) as tmp:
            tmp_path = Path(tmp)
            all_metrics = tmp_path / "all.tsv.gz"
            maf_sites = tmp_path / "maf.tsv.gz"
            for path in (all_metrics, maf_sites):
                with gzip.open(path, "wt", encoding="utf-8") as handle:
                    handle.writelines(rows)
            output = tmp_path / "audit.json"
            subprocess.run(
                [
                    sys.executable,
                    str(AUDITOR),
                    "--all-metrics",
                    str(all_metrics),
                    "--maf-sites",
                    str(maf_sites),
                    "--source-sha256",
                    "0" * 64,
                    "--output",
                    str(output),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            payload = json.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(payload["maf_gt_0.01_records"], 3)
            self.assertEqual(payload["duplicate_positions"], 1)
            self.assertEqual(payload["records_at_duplicate_positions"], 2)
            self.assertEqual(payload["status"], "PASS")


if __name__ == "__main__":
    unittest.main()
