from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "p0" / "freeze_samples_v7_1_1.py"


class FreezeSamplesTest(unittest.TestCase):
    def test_population_split_and_release_representative(self) -> None:
        fields = [
            "s",
            "hgdp_tgp_meta.Population",
            "hgdp_tgp_meta.Genetic.region",
            "subsets.tgp",
            "subsets.hgdp",
            "high_quality",
            "release",
            "sample_filters.release_related",
            "sample_filters.all_samples_related",
        ]
        rows = []
        for sample, population in (("T1", "P1"), ("T2", "P1"), ("T3", "P2"), ("T4", "P2")):
            rows.append([sample, population, "SP", "true", "false", "true", "true", "false", "false"])
        rows.extend([
            ["H1", "HP1", "HSP", "false", "true", "true", "true", "false", "true"],
            ["H2", "HP2", "HSP", "false", "true", "true", "true", "false", "false"],
        ])

        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parent) as tmp:
            tmp_path = Path(tmp)
            metadata = tmp_path / "metadata.tsv"
            with metadata.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(fields)
                writer.writerows(rows)
            samples = tmp_path / "samples.txt"
            samples.write_text("\n".join(row[0] for row in rows) + "\n", encoding="utf-8")
            output = tmp_path / "frozen"
            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--metadata",
                    str(metadata),
                    "--bcf-samples",
                    str(samples),
                    "--output-dir",
                    str(output),
                    "--expected-tgp",
                    "4",
                    "--expected-hgdp",
                    "2",
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            summary = json.loads((output / "SPLIT_SUMMARY.json").read_text(encoding="utf-8"))
            self.assertEqual(summary["partitions"]["donor_train"], 2)
            self.assertEqual(summary["partitions"]["donor_validation"], 2)
            self.assertEqual(summary["partitions"]["hgdp_primary"], 2)
            self.assertEqual(summary["release_all_samples_related"]["HGDP_ids"], ["H1"])
            calibration_lines = (output / "HGDP_SNPBAG_CALIBRATION.tsv").read_text(encoding="utf-8").splitlines()
            self.assertEqual(len(calibration_lines), 1)


if __name__ == "__main__":
    unittest.main()
