from __future__ import annotations

import gzip
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
FINALIZER = ROOT / "p0" / "finalize_sites_v7_1_3.py"


class FinalizeSitesTest(unittest.TestCase):
    def test_duplicate_rejection_and_monomorphic_retention(self) -> None:
        donor_rows = [
            "chr22\t100\tA\tG\t2\t8\t0.25\t0.25\t0\n",
            "chr22\t200\tC\tT\t2\t8\t0.25\t0.25\t0\n",
            "chr22\t200\tC\tG\t3\t8\t0.375\t0.375\t0\n",
            "chr22\t300\tG\tA\t2\t8\t0.25\t0.25\t0\n",
        ]
        hgdp_rows = [
            "chr22\t100\tA\tG\t0\t4\t0\t0\t0\n",
            "chr22\t200\tC\tT\t1\t4\t0.25\t0.25\t0\n",
            "chr22\t200\tC\tG\t1\t4\t0.25\t0.25\t0\n",
            "chr22\t300\tG\tA\t1\t4\t0.25\t0.25\t0\n",
        ]
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parent) as tmp:
            tmp_path = Path(tmp)
            donor = tmp_path / "donor.tsv.gz"
            hgdp = tmp_path / "hgdp.tsv.gz"
            for path, rows in ((donor, donor_rows), (hgdp, hgdp_rows)):
                with gzip.open(path, "wt", encoding="utf-8") as handle:
                    handle.writelines(rows)
            output = tmp_path / "final"
            subprocess.run(
                [
                    sys.executable,
                    str(FINALIZER),
                    "--donor-sites", str(donor),
                    "--hgdp-metrics", str(hgdp),
                    "--output-dir", str(output),
                    "--expected-donor-candidates", "4",
                    "--expected-duplicate-positions", "1",
                    "--expected-duplicate-records", "2",
                    "--expected-final", "2",
                    "--expected-hgdp-monomorphic", "1",
                    "--expected-hgdp-an", "4",
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            audit = json.loads((output / "FINAL_SITE_FREEZE.v7.1.3.json").read_text(encoding="utf-8"))
            self.assertEqual(audit["final_L"], 2)
            self.assertEqual(audit["records_at_duplicate_positions_rejected"], 2)
            self.assertEqual(audit["hgdp_monomorphic_reference_sites_retained"], 1)
            self.assertFalse(audit["hgdp_AC_used_for_selection"])
            self.assertEqual(
                (output / "FINAL_VARIANT_IDS.v7.1.3.txt").read_text(encoding="utf-8").splitlines(),
                ["chr22:100:A:G", "chr22:300:G:A"],
            )


if __name__ == "__main__":
    unittest.main()
