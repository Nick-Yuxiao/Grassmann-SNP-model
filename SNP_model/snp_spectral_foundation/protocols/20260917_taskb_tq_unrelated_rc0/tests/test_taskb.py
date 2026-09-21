#!/usr/bin/env python3
"""Synthetic end-to-end tests for the Task B qualification chain.

A small cohort is simulated with a known genetic component, written as real
PLINK 1 binary files, and pushed through B3 -> B4 -> B5. A trait with signal must
qualify; a trait that is pure noise must not.
"""

from __future__ import annotations

import csv
import io
import json
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

import numpy as np

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

import b2b_audit_flag_coding  # noqa: E402
import b3_build_unrelated_split  # noqa: E402
import b4_build_snp_panel  # noqa: E402
import b5_task_qualification  # noqa: E402
import b6_trait_preflight  # noqa: E402
import e2_bundle_task_gate  # noqa: E402

N_SAMPLES = 4000
N_VARIANTS = 600
N_CAUSAL = 40


def kinship_flag_value(index: int) -> str:
    """Every documented value of field 22021 plus a blank, at known counts."""
    if index % 10 == 1:
        return "1"            # at least one relative
    if index % 10 == 2:
        return "10"           # ten or more third-degree relatives
    if index % 20 == 3:
        return "-1"           # excluded from kinship inference: relatedness UNKNOWN
    if index % 20 == 7:
        return ""             # blank: relatedness UNKNOWN
    return "0"


def write_plink(directory: Path, stem: str, dosage: np.ndarray, chroms, positions) -> None:
    n_variants, n_samples = dosage.shape
    with (directory / f"{stem}.fam").open("w", encoding="utf-8") as handle:
        for index in range(n_samples):
            handle.write(f"F{index} IID{index} 0 0 {1 + index % 2} -9\n")
    with (directory / f"{stem}.bim").open("w", encoding="utf-8") as handle:
        for index in range(n_variants):
            handle.write(f"{chroms[index]} rs{index} 0 {positions[index]} A G\n")

    # count of A1 -> PLINK code: 2 -> 0b00, 1 -> 0b10, 0 -> 0b11
    code = np.select([dosage == 2, dosage == 1, dosage == 0], [0, 2, 3], default=1).astype(np.uint8)
    per_variant = (n_samples + 3) // 4
    payload = bytearray(b"\x6c\x1b\x01")
    for variant in range(n_variants):
        padded = np.zeros(per_variant * 4, dtype=np.uint8)
        padded[:n_samples] = code[variant]
        packed = (
            padded[0::4] | (padded[1::4] << 2) | (padded[2::4] << 4) | (padded[3::4] << 6)
        ).astype(np.uint8)
        payload.extend(packed.tobytes())
    (directory / f"{stem}.bed").write_bytes(bytes(payload))


def run_json(callable_, argv: list[str]) -> dict:
    buffer = io.StringIO()
    with redirect_stdout(buffer):
        code = callable_(argv)
    return {"exit_code": code, "stdout": json.loads(buffer.getvalue())}


class TaskBChain(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.temporary = tempfile.TemporaryDirectory()
        base = Path(cls.temporary.name)
        cls.base = base
        rng = np.random.default_rng(20260917)

        maf = rng.uniform(0.05, 0.5, N_VARIANTS)
        dosage = rng.binomial(2, maf[:, None], size=(N_VARIANTS, N_SAMPLES)).astype(np.int8)
        chroms = np.where(np.arange(N_VARIANTS) < N_VARIANTS // 2, "1", "2")
        positions = (np.arange(N_VARIANTS) % (N_VARIANTS // 2)) * 200_000 + 1_000_000

        cls.data = base / "data"
        cls.data.mkdir()
        write_plink(cls.data, "chrsim", dosage, chroms, positions)

        sex = np.array([1 + index % 2 for index in range(N_SAMPLES)], dtype=float)
        batch = rng.choice([-11, 2000], size=N_SAMPLES)
        pcs = rng.normal(size=(N_SAMPLES, 3))

        causal = rng.choice(N_VARIANTS, N_CAUSAL, replace=False)
        standardized = (dosage[causal].astype(float) - 2 * maf[causal, None]) / np.sqrt(
            2 * maf[causal, None] * (1 - maf[causal, None])
        )
        effects = rng.normal(0, 1.0, N_CAUSAL)
        genetic = effects @ standardized
        genetic = (genetic - genetic.mean()) / genetic.std()
        covariate_part = 0.6 * (sex - sex.mean()) + 0.4 * pcs[:, 0]
        signal_trait = covariate_part + 1.2 * genetic + rng.normal(0, 1.0, N_SAMPLES)
        noise_trait = covariate_part + rng.normal(0, 1.0, N_SAMPLES)

        flags = [kinship_flag_value(index) for index in range(N_SAMPLES)]
        cls.flag_counts = {value: flags.count(value) for value in {"0", "1", "10", "-1", ""}}
        cls.unrelated_count = cls.flag_counts["0"]
        cls.dropped_flag_count = N_SAMPLES - cls.unrelated_count

        cls.covariates = base / "covariates.tsv"
        with cls.covariates.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["eid", "iid", "n_22001_0_0", "n_22000_0_0",
                             "n_22009_0_1", "n_22009_0_2", "n_22009_0_3"])
            for index in range(N_SAMPLES):
                writer.writerow([f"E{index}", f"IID{index}", int(sex[index]), int(batch[index]),
                                 *[round(float(value), 6) for value in pcs[index]]])

        cls.kinship = base / "kinship_flag.tsv"
        with cls.kinship.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["eid", "n_22021_0_0", "n_22027_0_0"])
            for index in range(N_SAMPLES):
                writer.writerow([f"E{index}", flags[index], ""])

        cls.phenotype = base / "phenotype.tsv"
        with cls.phenotype.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["eid", "n_30020_0_0", "n_99999_0_0", "n_88888_0_0"])
            for index in range(N_SAMPLES):
                # n_88888 is deliberately degenerate: sentinel codes and heavy rounding.
                degenerate = -3 if index % 4 == 0 else round(float(signal_trait[index]))
                writer.writerow([f"E{index}", round(float(signal_trait[index]), 6),
                                 round(float(noise_trait[index]), 6), degenerate])

        cls.split_dir = base / "split"
        cls.split_result = run_json(b3_build_unrelated_split.main, [
            "--covariates", str(cls.covariates),
            "--covariate-id-column", "eid",
            "--genotype-id-column", "iid",
            "--fam", str(cls.data / "chrsim.fam"),
            "--kinship-flag", str(cls.kinship),
            "--qc-flag", "n_22027_0_0",
            "--strata-column", "n_22001_0_0",
            "--seed", "20260917",
            "--out-dir", str(cls.split_dir),
        ])

        cls.panel_dir = base / "panel"
        cls.panel_result = run_json(b4_build_snp_panel.main, [
            "--bed", str(cls.data / "chrsim.bed"),
            "--split-manifest", str(cls.split_dir / "tq_split_manifest.tsv"),
            "--out-dir", str(cls.panel_dir),
            "--thin-bp", "1000",
            "--maf-min", "0.01",
        ])

    @classmethod
    def tearDownClass(cls) -> None:
        cls.temporary.cleanup()

    def qualification(self, trait: str, out_name: str, allow_test: bool = True) -> dict:
        argv = [
            "--panel-dir", str(self.panel_dir),
            "--covariates", str(self.covariates),
            "--covariate-column", "n_22009_0_1",
            "--covariate-column", "n_22009_0_2",
            "--covariate-column", "n_22009_0_3",
            "--categorical-column", "n_22001_0_0",
            "--phenotype", str(self.phenotype),
            "--trait-column", trait,
            "--phenotype-transform", "zscore",
            "--bootstrap", "200",
            "--out-dir", str(self.base / out_name),
        ]
        if allow_test:
            argv.append("--allow-test")
        return run_json(b5_task_qualification.main, argv)["stdout"]

    def test_related_participants_are_dropped(self) -> None:
        summary = json.loads((self.split_dir / "SPLIT_SUMMARY.json").read_text(encoding="utf-8"))
        self.assertEqual(
            summary["filtering"]["dropped_related_or_unassessed"], self.dropped_flag_count
        )
        self.assertEqual(self.split_result["stdout"]["retained"], self.unrelated_count)
        self.assertFalse(summary["relatedness_axis"]["pairwise_edges_available"])

    def test_unknown_relatedness_is_dropped_not_kept(self) -> None:
        axis = json.loads(
            (self.split_dir / "SPLIT_SUMMARY.json").read_text(encoding="utf-8")
        )["relatedness_axis"]
        self.assertEqual(axis["kept_values"], ["0"])
        self.assertTrue(axis["unknown_relatedness_is_dropped_not_kept"])
        observed = axis["observed_value_counts"]
        self.assertEqual(observed["0"], self.flag_counts["0"])
        self.assertEqual(observed["<blank>"], self.flag_counts[""])
        # -1 means relatedness was never assessed; it must be dropped, not read as unrelated.
        self.assertEqual(axis["dropped_by_value"]["-1"], self.flag_counts["-1"])
        self.assertEqual(axis["dropped_by_value"]["<blank>"], self.flag_counts[""])
        self.assertEqual(axis["dropped_by_value"]["10"], self.flag_counts["10"])

    def test_flag_audit_separates_unknown_from_related(self) -> None:
        result = run_json(b2b_audit_flag_coding.main, [
            "--table", str(self.kinship),
            "--id-column", "eid",
            "--column", "n_22021_0_0",
            "--cross-reference", str(self.covariates),
            "--cross-id-column", "eid",
            "--out-dir", str(self.base / "flag_audit"),
        ])
        summary = result["stdout"]["n_22021_0_0"]
        self.assertEqual(summary["keep_known_unrelated"], self.flag_counts["0"])
        self.assertEqual(
            summary["drop_known_related"], self.flag_counts["1"] + self.flag_counts["10"]
        )
        self.assertEqual(
            summary["drop_relatedness_unknown"], self.flag_counts["-1"] + self.flag_counts[""]
        )
        self.assertEqual(summary["undocumented_values_present"], [])
        self.assertEqual(result["stdout"]["participant_overlap"]["in_both"], N_SAMPLES)

    def test_bridge_holdout_split_stays_reserved(self) -> None:
        result = run_json(b3_build_unrelated_split.main, [
            "--covariates", str(self.covariates),
            "--covariate-id-column", "eid",
            "--genotype-id-column", "iid",
            "--fam", str(self.data / "chrsim.fam"),
            "--kinship-flag", str(self.kinship),
            "--strata-column", "n_22001_0_0",
            "--reserve-bridge-holdout",
            "--seed", "20260917",
            "--out-dir", str(self.base / "split4"),
        ])
        counts = result["stdout"]["split_sample_counts"]
        self.assertEqual(
            sorted(counts), ["bridge_holdout", "tq_test", "train", "validation"]
        )
        total = sum(counts.values())
        self.assertAlmostEqual(counts["train"] / total, 0.55, delta=0.02)
        self.assertAlmostEqual(counts["bridge_holdout"] / total, 0.15, delta=0.02)

    def test_qualification_never_touches_the_bridge_holdout(self) -> None:
        panel_dir = self.base / "panel4"
        run_json(b4_build_snp_panel.main, [
            "--bed", str(self.data / "chrsim.bed"),
            "--split-manifest", str(self.base / "split4" / "tq_split_manifest.tsv"),
            "--out-dir", str(panel_dir), "--thin-bp", "1000",
        ])
        report = run_json(b5_task_qualification.main, [
            "--panel-dir", str(panel_dir),
            "--covariates", str(self.covariates),
            "--covariate-column", "n_22009_0_1",
            "--categorical-column", "n_22001_0_0",
            "--phenotype", str(self.phenotype),
            "--trait-column", "n_30020_0_0",
            "--test-split", "tq_test",
            "--bootstrap", "100",
            "--allow-test",
            "--out-dir", str(self.base / "tq4"),
        ])["stdout"]
        full = json.loads((self.base / "tq4" / "TQ_RESULTS.json").read_text(encoding="utf-8"))
        self.assertEqual(full["splits"]["never_touched"], ["bridge_holdout"])
        self.assertEqual(full["splits"]["test"], "tq_test")
        self.assertEqual(report["verdict"], "QUALIFIED")

    def test_preflight_screens_traits_and_keeps_test_closed(self) -> None:
        result = run_json(b6_trait_preflight.main, [
            "--phenotype", str(self.phenotype),
            "--phenotype-id-column", "eid",
            "--trait-column", "n_30020_0_0",
            "--trait-column", "n_88888_0_0",
            "--covariates", str(self.covariates),
            "--covariate-id-column", "eid",
            "--covariate-column", "n_22009_0_1",
            "--categorical-column", "n_22001_0_0",
            "--split-manifest", str(self.split_dir / "tq_split_manifest.tsv"),
            "--min-effective-n", "500",
            "--out-dir", str(self.base / "preflight"),
        ])["stdout"]
        self.assertEqual(result["splits_read"], ["train", "validation"])
        self.assertEqual(result["splits_left_closed"], ["test"])
        self.assertIn("n_30020_0_0", result["recommended_traits"])
        self.assertNotIn("n_88888_0_0", result["recommended_traits"])

        report = json.loads(
            (self.base / "preflight" / "PREFLIGHT.json").read_text(encoding="utf-8")
        )
        degenerate = next(item for item in report["traits"]
                          if item["trait_column"] == "n_88888_0_0")
        self.assertFalse(degenerate["usable"])
        # Heavy rounding plus a dominant value is what disqualifies it. The negative
        # integers are reported but not treated as sentinels, because this trait's median
        # is not positive, so going below zero is not by itself suspicious.
        self.assertTrue(any("distinct values" in reason for reason in degenerate["reasons"]))
        self.assertTrue(any("dominates" in reason for reason in degenerate["reasons"]))
        self.assertIn(-3.0, degenerate["negative_integer_modes"])
        self.assertEqual(degenerate["sentinel_suspects"], [])

    def test_split_is_singleton_and_proportional(self) -> None:
        with (self.split_dir / "tq_split_manifest.tsv").open(encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertTrue(all(row["component_size"] == "1" for row in rows))
        self.assertEqual(len({row["sample_id"] for row in rows}), len(rows))
        counts = self.split_result["stdout"]["split_sample_counts"]
        self.assertAlmostEqual(counts["train"] / len(rows), 0.70, delta=0.02)
        self.assertAlmostEqual(counts["test"] / len(rows), 0.15, delta=0.02)

    def test_panel_decodes_real_plink_bytes(self) -> None:
        panel = np.load(self.panel_dir / "panel.int8.npy", mmap_mode="r")
        summary = json.loads((self.panel_dir / "PANEL_SUMMARY.json").read_text(encoding="utf-8"))
        self.assertEqual(panel.shape[0], summary["shape"]["variants"])
        self.assertEqual(panel.shape[1], summary["shape"]["samples"])
        self.assertTrue(((panel >= 0) & (panel <= 2)).all())
        self.assertGreater(panel.shape[0], 400)

    def test_signal_trait_qualifies(self) -> None:
        report = self.qualification("n_30020_0_0", "tq_signal")
        self.assertEqual(report["verdict"], "QUALIFIED")
        self.assertGreater(report["test"]["delta_r2"], 0.005)
        self.assertGreater(report["test"]["delta_r2_ci95"][0], 0)
        # Scrambling the score against the outcome must not look like a gain; with real
        # signal the model becomes actively worse than covariates alone.
        negative = report["test"]["negative_control_permuted_score_delta_r2"]
        self.assertLess(negative, 0.0)
        self.assertLess(negative, report["test"]["delta_r2"])

    def test_noise_trait_does_not_qualify(self) -> None:
        report = self.qualification("n_99999_0_0", "tq_noise")
        self.assertNotEqual(report["verdict"], "QUALIFIED")

    def test_test_split_stays_closed_without_the_flag(self) -> None:
        report = self.qualification("n_30020_0_0", "tq_closed", allow_test=False)
        self.assertEqual(report["verdict"], "VALIDATION_ONLY_TEST_NOT_OPENED")
        self.assertNotIn("test", report)
        self.assertFalse((self.base / "tq_closed" / "TEST_OPENED.marker").exists())

    def test_export_bundle_is_auditable_and_keeps_test_shut(self) -> None:
        export = self.base / "export"
        for trait in ("n_30020_0_0", "n_99999_0_0"):
            run_json(b5_task_qualification.main, [
                "--panel-dir", str(self.panel_dir),
                "--covariates", str(self.covariates),
                "--covariate-column", "n_22009_0_1",
                "--categorical-column", "n_22001_0_0",
                "--phenotype", str(self.phenotype),
                "--trait-column", trait,
                "--split-manifest", str(self.split_dir / "tq_split_manifest.tsv"),
                "--export-dir", str(export),
                "--export-bootstrap", "50",
                "--out-dir", str(self.base / f"export_run_{trait}"),
            ])

        with (export / "n_30020_0_0_validation_curve.csv").open(encoding="utf-8") as handle:
            curve = list(csv.DictReader(handle))
        for column in ("delta_r2", "delta_r2_ci95_low", "delta_r2_ci95_high",
                       "pearson_A", "pearson_B", "n_snps", "selected"):
            self.assertIn(column, curve[0])
        # Exactly one threshold is flagged as the chosen one, and every row is
        # evaluated on validation rather than the held-out test split.
        self.assertEqual(sum(int(row["selected"]) for row in curve), 1)
        self.assertTrue(all(row["split_evaluated"] == "validation" for row in curve))
        for row in curve:
            self.assertLessEqual(float(row["delta_r2_ci95_low"]), float(row["delta_r2"]))
            self.assertGreaterEqual(float(row["delta_r2_ci95_high"]), float(row["delta_r2"]))

        # One participant keeps one surrogate id across traits, so the tables join,
        # and only validation rows are exported.
        def read_predictions(trait: str) -> list[dict[str, str]]:
            with (export / f"validation_predictions_{trait}.csv").open(encoding="utf-8") as fh:
                return list(csv.DictReader(fh))

        first, second = read_predictions("n_30020_0_0"), read_predictions("n_99999_0_0")
        self.assertEqual({row["surrogate_id"] for row in first},
                         {row["surrogate_id"] for row in second})
        self.assertTrue(all(row["split"] == "validation" for row in first))
        self.assertNotIn("IID0", (export / "validation_predictions_n_30020_0_0.csv")
                         .read_text(encoding="utf-8"))

        config = self.base / "CONFIG.json"
        config.write_text(json.dumps({
            "firewall": {"test_opened": False},
            "cohort": {"withdrawal_list_supplied": False},
            "panel": {"chromosomes_included": [1], "chromosomes_excluded": []},
        }), encoding="utf-8")
        panel_summary = self.panel_dir / "PANEL_SUMMARY.json"
        bundle = run_json(e2_bundle_task_gate.main, [
            "--config", str(config),
            "--split-summary", str(self.split_dir / "SPLIT_SUMMARY.json"),
            "--panel-summary", str(panel_summary),
            "--export-dir", str(export),
            "--trait", "n_30020_0_0",
            "--out-dir", str(self.base / "task_gate"),
            "--include-predictions",
        ])["stdout"]
        self.assertFalse(bundle["test_opened"])
        self.assertIn("hashes.txt", bundle["files"])
        self.assertIn("exclusions.json", bundle["files"])
        self.assertIn("panel_variants.tsv", bundle["files"])
        readme = (self.base / "task_gate" / "README.md").read_text(encoding="utf-8")
        # The bundle must not oversell what it proves.
        self.assertIn("not proof that the implementation is leak-free", readme)
        self.assertIn("pseudonymisation, not anonymisation", readme)

        # A trait whose b5 export has not finished must leave nothing behind, so a
        # half-written bundle cannot be tarred up and shipped as the real thing.
        partial = self.base / "task_gate_partial"
        with self.assertRaises(SystemExit) as caught:
            e2_bundle_task_gate.main([
                "--config", str(config),
                "--split-summary", str(self.split_dir / "SPLIT_SUMMARY.json"),
                "--panel-summary", str(panel_summary),
                "--export-dir", str(export),
                "--trait", "n_30020_0_0", "--trait", "n_30080_0_0",
                "--out-dir", str(partial),
            ])
        self.assertIn("n_30080_0_0", str(caught.exception))
        self.assertIn("b5 --export-dir", str(caught.exception))
        self.assertFalse(partial.exists())

    def test_reopening_the_test_split_is_refused(self) -> None:
        self.qualification("n_30020_0_0", "tq_once")
        with self.assertRaises(SystemExit):
            self.qualification("n_30020_0_0", "tq_once")


class RareVariantPanelTests(unittest.TestCase):
    """An imputed panel is mostly rare variants; windows must still fill."""

    WINDOWS = 20
    PER_WINDOW = 10
    SAMPLES = 2000
    WINDOW_BP = 100_000

    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        base = Path(self.temporary.name)
        self.base = base
        rng = np.random.default_rng(4242)

        n_variants = self.WINDOWS * self.PER_WINDOW
        dosage = np.zeros((n_variants, self.SAMPLES), dtype=np.int8)
        positions = []
        chroms = []
        self.common_rows = set()
        for window in range(self.WINDOWS):
            # Exactly one common variant per window, never the first candidate.
            common_slot = 3 + window % (self.PER_WINDOW - 3)
            for slot in range(self.PER_WINDOW):
                row = window * self.PER_WINDOW + slot
                positions.append(window * self.WINDOW_BP + slot * 5_000 + 1_000)
                chroms.append("1")
                if slot == common_slot:
                    dosage[row] = rng.binomial(2, 0.3, self.SAMPLES)
                    self.common_rows.add(row)
                else:
                    carriers = rng.choice(self.SAMPLES, 2, replace=False)
                    dosage[row, carriers] = 1

        self.data = base / "data"
        self.data.mkdir()
        write_plink(self.data, "rare", dosage, chroms, positions)

        self.manifest = base / "split.tsv"
        with self.manifest.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["sample_id", "eid", "split"])
            for index in range(self.SAMPLES):
                writer.writerow([f"IID{index}", f"E{index}",
                                 "train" if index % 4 else "tq_test"])

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_every_window_fills_by_retrying(self) -> None:
        result = run_json(b4_build_snp_panel.main, [
            "--bed", str(self.data / "rare.bed"),
            "--split-manifest", str(self.manifest),
            "--out-dir", str(self.base / "panel"),
            "--thin-bp", str(self.WINDOW_BP),
            "--maf-min", "0.01",
        ])["stdout"]
        # One kept variant per window, and it is the common one, not the leftmost.
        self.assertEqual(result["variants_kept"], self.WINDOWS)
        self.assertEqual(result["windows_unfilled"], 0)
        self.assertGreater(result["reads_per_kept_variant"], 1.0)
        self.assertGreater(result["dropped"]["maf"], 0)

        with (self.base / "panel" / "panel_variants.tsv").open(encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), self.WINDOWS)
        for row in rows:
            self.assertGreaterEqual(float(row["train_maf"]), 0.01)
        self.assertEqual(
            {int(row["bim_index"]) for row in rows}, self.common_rows
        )

    def test_truncated_bed_is_refused_before_any_genotype_read(self) -> None:
        # This cohort keeps duplicate copies of some chromosomes and at least one is
        # incomplete. Catching it up front costs a bim line count; catching it late
        # costs an hour of genotype reading.
        bed = self.data / "rare.bed"
        with bed.open("r+b") as handle:
            handle.truncate(bed.stat().st_size - 500)
        out_dir = self.base / "panel_truncated"
        with self.assertRaises(SystemExit) as caught:
            run_json(b4_build_snp_panel.main, [
                "--bed", str(bed),
                "--split-manifest", str(self.manifest),
                "--out-dir", str(out_dir),
                "--thin-bp", str(self.WINDOW_BP),
            ])
        message = str(caught.exception)
        self.assertIn("incomplete", message)
        self.assertIn("rare.bed", message)
        self.assertFalse((out_dir / "panel.tmp.npy").exists())

    def test_windows_go_unfilled_when_retries_are_capped(self) -> None:
        result = run_json(b4_build_snp_panel.main, [
            "--bed", str(self.data / "rare.bed"),
            "--split-manifest", str(self.manifest),
            "--out-dir", str(self.base / "panel_capped"),
            "--thin-bp", str(self.WINDOW_BP),
            "--maf-min", "0.01",
            "--max-tries-per-window", "1",
        ])["stdout"]
        self.assertLess(result["variants_kept"], self.WINDOWS)
        self.assertGreater(result["windows_unfilled"], 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
