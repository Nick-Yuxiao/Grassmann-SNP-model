#!/usr/bin/env python3
"""Tests for the complementarity Concept Gate package.

The fixture is a small synthetic cohort where the answer is known by construction:
the phenotype is built from covariates plus a handful of panel variants, and the
"functional" annotation marks exactly those variants. A correct gate must then show
F adding over covariates, and must not show a carrier adding when the annotation
points at noise.
"""

from __future__ import annotations

import csv
import json
import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import c0_prepare_annotations as c0  # noqa: E402
import c1_build_population_features as c1  # noqa: E402
import c1b_overlap_spectrum as c1b  # noqa: E402
import c2_build_functional_features as c2  # noqa: E402
import c3_complementarity_gate as c3  # noqa: E402
import lib_arms  # noqa: E402


def write_tsv(path: Path, fields: list[str], rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def build_cohort(root: Path, n_samples: int = 1400, n_variants: int = 120, seed: int = 7):
    """A panel, covariates and a phenotype with known genetic and covariate parts."""
    rng = np.random.default_rng(seed)
    root.mkdir(parents=True, exist_ok=True)
    panel_dir = root / "panel"
    panel_dir.mkdir(exist_ok=True)

    frequency = rng.uniform(0.15, 0.5, n_variants)
    dosage = rng.binomial(2, frequency[:, None], size=(n_variants, n_samples)).astype(np.int8)
    dosage[3, :5] = -1  # a few missing calls, the panel's sentinel
    np.save(panel_dir / "panel.int8.npy", dosage)

    # Two chromosomes so block building has a boundary to respect.
    variants = []
    for index in range(n_variants):
        chrom = "1" if index < n_variants // 2 else "2"
        position = 100000 + (index % (n_variants // 2)) * 5000
        observed = dosage[index][dosage[index] >= 0]
        variants.append({
            "panel_row": index, "bed": "synthetic.bed", "bim_index": index,
            "chrom": chrom, "pos": position, "variant_id": f"rs{index}",
            "a1": "A", "a2": "G",
            "train_a1_frequency": round(float(observed.mean()) / 2.0, 6),
            "train_maf": round(min(float(observed.mean()) / 2.0,
                                   1 - float(observed.mean()) / 2.0), 6),
            "train_missing_rate": 0.0,
        })
    write_tsv(panel_dir / "panel_variants.tsv", list(variants[0]), variants)

    splits = np.array(["train"] * n_samples, dtype=object)
    splits[int(n_samples * 0.6):int(n_samples * 0.8)] = "validation"
    splits[int(n_samples * 0.8):int(n_samples * 0.9)] = "tq_test"
    splits[int(n_samples * 0.9):] = "bridge_holdout"

    samples = [{
        "panel_column": column, "sample_id": f"S{column}", "eid": str(100000 + column),
        "split": splits[column], "stratum": "0|1",
    } for column in range(n_samples)]
    write_tsv(panel_dir / "panel_samples.tsv", list(samples[0]), samples)
    write_tsv(root / "split_manifest.tsv", ["sample_id", "split"],
              [{"sample_id": row["sample_id"], "split": row["split"]} for row in samples])

    age = rng.normal(58, 8, n_samples)
    sex = rng.integers(0, 2, n_samples).astype(float)
    centre = rng.integers(0, 3, n_samples)
    causal = np.arange(5, 25)  # the variants the annotation will mark
    scaled = lib_arms.standardise_block(dosage[causal].T, frequency[causal])
    # Same-sign effects, so an unweighted annotation aggregate can actually see them.
    # With random signs the sum would be near-orthogonal to the genetic value and the
    # F arm would correctly find nothing, which tests the wrong thing here.
    genetic = scaled @ (np.abs(rng.normal(0, 1.0, causal.size)) + 0.5)
    y = (0.05 * (age - 58) + 0.8 * sex + 0.6 * genetic / genetic.std()
         + rng.normal(0, 1.0, n_samples))

    write_tsv(root / "covariates.tsv", ["eid", "age", "sex", "centre"], [
        {"eid": samples[i]["eid"], "age": round(float(age[i]), 4),
         "sex": int(sex[i]), "centre": f"C{centre[i]}"} for i in range(n_samples)
    ])
    write_tsv(root / "phenotype.tsv", ["eid", "trait"], [
        {"eid": samples[i]["eid"], "trait": round(float(y[i]), 6)} for i in range(n_samples)
    ])

    annotation = root / "causal.bed"
    with annotation.open("w", encoding="utf-8") as handle:
        for index in causal:
            row = variants[index]
            handle.write(f"{row['chrom']}\t{row['pos'] - 1}\t{row['pos']}\tcausal\t1\n")
    write_tsv(root / "annotation_manifest.tsv",
              ["track", "path", "kind", "value_column", "phenotype_free", "source"],
              [{"track": "causal", "path": str(annotation), "kind": "bed",
                "value_column": "score", "phenotype_free": "true", "source": "synthetic"}])
    return panel_dir


class TestLibArms(unittest.TestCase):
    def test_missing_dosage_becomes_the_train_mean(self):
        dosage = np.array([[0, 1, 2, -1]], dtype=np.int8).T
        scaled = lib_arms.standardise_block(dosage, np.array([0.5]))
        self.assertAlmostEqual(float(scaled[3, 0]), 0.0)

    def test_ridge_leaves_the_intercept_unpenalised(self):
        design = np.hstack([np.ones((50, 1)), np.linspace(-1, 1, 50).reshape(-1, 1)])
        y = 7.0 + np.zeros(50)
        beta = lib_arms.fit_ridge(design.T @ design, design.T @ y, penalty=1e6)
        # A huge penalty crushes the slope but must leave the mean intact.
        self.assertAlmostEqual(float(beta[0]), 7.0, places=6)
        self.assertLess(abs(float(beta[1])), 1e-6)

    def test_inner_cv_prefers_a_light_penalty_for_a_strong_signal(self):
        rng = np.random.default_rng(3)
        design = np.hstack([np.ones((400, 1)), rng.normal(size=(400, 4))])
        y = design[:, 1] * 3.0 + rng.normal(0, 0.05, 400)
        penalty, trace = lib_arms.select_penalty_by_inner_cv(design, y, folds=4, seed=1)
        self.assertLess(penalty, 10.0)
        self.assertEqual(len(trace), len(lib_arms.LAMBDA_GRID))

    def test_paired_bootstrap_of_identical_predictions_is_zero(self):
        rng = np.random.default_rng(5)
        y = rng.normal(size=300)
        prediction = rng.normal(size=300)
        result = lib_arms.paired_bootstrap(y, prediction, prediction, 0.0, 60, 11)
        self.assertAlmostEqual(result["delta_r2"], 0.0, places=12)
        self.assertAlmostEqual(result["ci_low"], 0.0, places=12)

    def test_canonical_correlation_of_a_block_with_itself_is_one(self):
        rng = np.random.default_rng(9)
        block = rng.normal(size=(200, 3))
        values = lib_arms.canonical_correlations(block, block, top=3)
        self.assertTrue(all(abs(value - 1.0) < 1e-6 for value in values))

    def test_out_of_sample_r2_can_be_negative(self):
        y = np.array([1.0, 2.0, 3.0])
        self.assertLess(lib_arms.out_of_sample_r2(y, np.full(3, 10.0), float(y.mean())), 0.0)


class TestPopulationFeatures(unittest.TestCase):
    def test_blocks_never_span_a_chromosome(self):
        variants = [{"chrom": "1"}] * 7 + [{"chrom": "2"}] * 5
        ranges = c1.block_ranges(variants, block_size=4)
        for start, stop, chrom in ranges:
            self.assertTrue(all(variants[i]["chrom"] == chrom for i in range(start, stop)))
        self.assertEqual(sum(stop - start for start, stop, _ in ranges), len(variants))

    def test_end_to_end_writes_train_fitted_components(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            panel_dir = build_cohort(root)
            out = root / "P"
            self.assertEqual(c1.main([
                "--panel-dir", str(panel_dir),
                "--split-manifest", str(root / "split_manifest.tsv"),
                "--block-size", "30", "--pcs-per-block", "2",
                "--out-dir", str(out),
            ]), 0)
            features = np.load(out / "P_population.npy")
            self.assertEqual(features.shape[0], 1400)
            self.assertEqual(features.shape[1], 8)  # 4 blocks x 2 components
            summary = json.loads((out / "P_SUMMARY.json").read_text())
            self.assertFalse(summary["phenotype_read"])
            self.assertFalse(summary["identifiers_included"])


class TestOverlapSpectrum(unittest.TestCase):
    def test_a_pc_copied_from_P_shows_up_as_a_pair_at_one(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            panel_dir = build_cohort(root)
            out = root / "P"
            c1.main(["--panel-dir", str(panel_dir),
                     "--split-manifest", str(root / "split_manifest.tsv"),
                     "--block-size", "30", "--pcs-per-block", "2", "--out-dir", str(out)])
            features = np.load(out / "P_population.npy")
            samples = list(c1b.read_tsv(out / "P_samples.tsv"))

            # One "global PC" is a column of P exactly, the rest are noise. The first
            # canonical pair must then be 1, and the others clearly below it.
            rng = np.random.default_rng(31)
            noise = rng.normal(size=(features.shape[0], 2))
            write_tsv(root / "pcs.tsv",
                      ["eid", "n_22009_0_1", "n_22009_0_2", "n_22009_0_3"],
                      [{"eid": row["eid"],
                        "n_22009_0_1": float(features[i, 0]),
                        "n_22009_0_2": float(noise[i, 0]),
                        "n_22009_0_3": float(noise[i, 1])}
                       for i, row in enumerate(samples)])

            self.assertEqual(c1b.main([
                "--population-dir", str(out),
                "--covariates", str(root / "pcs.tsv"),
                "--global-pc-count", "3",
            ]), 0)
            report = json.loads((out / "P_OVERLAP_SPECTRUM.json").read_text())
            self.assertEqual(report["canonical_pairs"], 3)
            self.assertAlmostEqual(report["spectrum"][0], 1.0, places=4)
            self.assertGreaterEqual(report["pairs_above_threshold"]["above_0.99"], 1)
            self.assertFalse(report["phenotype_read"])
            self.assertTrue((out / "P_OVERLAP_SPECTRUM.csv").exists())


class TestFunctionalFeatures(unittest.TestCase):
    def test_bed_lookup_is_half_open(self):
        index = c2.load_bed(self._bed("1\t100\t110\tx\t2.5\n"))
        self.assertEqual(c2.bed_value(index, "1", 101), 2.5)  # 1-based 101 is 0-based 100
        self.assertEqual(c2.bed_value(index, "1", 110), 2.5)  # 0-based 109, inside
        self.assertEqual(c2.bed_value(index, "1", 111), 0.0)  # 0-based 110, past the end

    def _bed(self, content: str) -> Path:
        import tempfile

        handle = tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False)
        handle.write(content)
        handle.close()
        return Path(handle.name)

    def test_a_track_not_declared_phenotype_free_is_refused(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            panel_dir = build_cohort(root, n_samples=400, n_variants=40)
            manifest = root / "bad_manifest.tsv"
            write_tsv(manifest, ["track", "path", "kind", "value_column", "phenotype_free",
                                 "source"],
                      [{"track": "gwas_hits", "path": str(root / "causal.bed"), "kind": "bed",
                        "value_column": "score", "phenotype_free": "false",
                        "source": "gwas catalogue"}])
            with self.assertRaises(SystemExit) as caught:
                c2.main(["--panel-dir", str(panel_dir),
                         "--split-manifest", str(root / "split_manifest.tsv"),
                         "--annotation-manifest", str(manifest),
                         "--out-dir", str(root / "F")])
            self.assertIn("phenotype_free", str(caught.exception))

    def test_end_to_end_aggregates_dosage_under_the_annotation(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            panel_dir = build_cohort(root)
            out = root / "F"
            self.assertEqual(c2.main([
                "--panel-dir", str(panel_dir),
                "--split-manifest", str(root / "split_manifest.tsv"),
                "--annotation-manifest", str(root / "annotation_manifest.tsv"),
                "--group-by", "chromosome", "--out-dir", str(out),
            ]), 0)
            features = np.load(out / "F_functional.npy")
            self.assertEqual(features.shape[0], 1400)
            summary = json.loads((out / "F_SUMMARY.json").read_text())
            self.assertTrue(summary["all_tracks_declared_phenotype_free"])
            self.assertEqual(summary["join"],
                             "variant position only; allele orientation is not used")
            self.assertGreater(summary["tracks"][0]["variants_annotated"], 0)


class TestAnnotationPrep(unittest.TestCase):
    def test_an_association_derived_local_track_is_refused(self):
        with self.assertRaises(SystemExit) as caught:
            c0.main(["--panel-dir", ".", "--out-dir", ".", "--offline",
                     "--local", "bcx_gwas_hits=/tmp/x.bed"])
        self.assertIn("association-derived", str(caught.exception))

    def test_the_regex_catches_the_shapes_that_matter(self):
        for name in ["gwas_catalog", "ukb_PRS_weights", "whole_blood_eQTL", "twas_z",
                     "ldsc_baseline", "bcx_sumstats"]:
            self.assertIsNotNone(c0.ASSOCIATION_DERIVED.search(name), name)
        for name in ["k562_chromhmm_broad", "phastcons100way", "dnase_clusters",
                     "tss_proximity"]:
            self.assertIsNone(c0.ASSOCIATION_DERIVED.search(name), name)

    def test_ucsc_tables_drop_their_leading_bin_column(self):
        import tempfile

        handle = tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False)
        handle.write("585\tchr1\t100\t200\tTxn_Elongation\t1000\n")
        handle.close()
        intervals = c0.to_intervals(Path(handle.name), "ucsc_bed_bin")
        self.assertEqual(intervals, [("1", 100, 200, "Txn_Elongation", 1000.0)])

    def test_every_catalogue_track_carries_a_justification(self):
        for name, entry in c0.CATALOGUE.items():
            self.assertTrue(str(entry["why_phenotype_free"]).strip(), name)
            self.assertTrue(str(entry["caveat"]).strip(), name)
            self.assertTrue(str(entry["url"]).startswith("https://"), name)


class TestConceptGate(unittest.TestCase):
    def _run(self, root: Path, panel_dir: Path, extra: list[str] | None = None):
        c1.main(["--panel-dir", str(panel_dir),
                 "--split-manifest", str(root / "split_manifest.tsv"),
                 "--block-size", "30", "--pcs-per-block", "2", "--out-dir", str(root / "P")])
        c2.main(["--panel-dir", str(panel_dir),
                 "--split-manifest", str(root / "split_manifest.tsv"),
                 "--annotation-manifest", str(root / "annotation_manifest.tsv"),
                 "--out-dir", str(root / "F")])
        return c3.main([
            "--panel-dir", str(panel_dir),
            "--covariates", str(root / "covariates.tsv"),
            "--covariate-column", "age", "--covariate-column", "sex",
            "--square-column", "age", "--categorical-column", "centre",
            "--phenotype", str(root / "phenotype.tsv"), "--trait-column", "trait",
            "--population-features", str(root / "P" / "P_population.npy"),
            "--functional-features", str(root / "F" / "F_functional.npy"),
            "--frozen-threshold", "0.05",
            "--split-manifest", str(root / "split_manifest.tsv"),
            "--inner-folds", "3", "--bootstrap", "60",
            "--out-dir", str(root / "gate"), *(extra or []),
        ])

    def test_sealed_splits_cannot_be_evaluated(self):
        for split in ("tq_test", "bridge_holdout"):
            with self.assertRaises(SystemExit) as caught:
                c3.main(["--panel-dir", ".", "--covariates", ".", "--phenotype", ".",
                         "--trait-column", "t", "--population-features", ".",
                         "--functional-features", ".", "--out-dir", ".",
                         "--eval-split", split, "--frozen-threshold", "0.05"])
            self.assertIn("sealed", str(caught.exception))

    def test_exactly_one_threshold_source_is_required(self):
        base = ["--panel-dir", ".", "--covariates", ".", "--phenotype", ".",
                "--trait-column", "t", "--population-features", ".",
                "--functional-features", ".", "--out-dir", "."]
        with self.assertRaises(SystemExit) as caught:
            c3.main(base)
        self.assertIn("exactly one", str(caught.exception))
        with self.assertRaises(SystemExit) as caught:
            c3.main(base + ["--frozen-threshold", "0.05", "--tq-results", "x.json"])
        self.assertIn("exactly one", str(caught.exception))

    def test_end_to_end_reports_all_five_arms_and_no_individual_rows(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            panel_dir = build_cohort(root)
            self.assertEqual(self._run(root, panel_dir), 0)
            report = json.loads((root / "gate" / "GATE_RESULTS.json").read_text())

            self.assertEqual(report["status"], "DEVELOPMENTAL_NOT_CONFIRMATORY")
            self.assertFalse(report["individual_level_rows_exported"])
            self.assertEqual({row["arm"] for row in report["arms"]},
                             {"A", "B", "P", "F", "PF"})
            self.assertFalse(report["frozen_threshold"]["reselected_here"])

            # Sealed participants must never reach the design.
            self.assertEqual(report["participants"]["usable"],
                             report["participants"]["train"] + report["participants"]["eval"])
            self.assertGreater(report["participants"]["dropped"]["sealed_split"], 0)

            written = {path.name for path in (root / "gate").iterdir()}
            self.assertEqual(written, {"GATE_RESULTS.json", "PENALTY_TRACE.json",
                                       "arms.csv", "comparisons.csv"})
            with (root / "gate" / "comparisons.csv").open() as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual({row["comparison"] for row in rows},
                             {label for label, _, _, _ in c3.COMPARISONS})

    def test_a_carrier_built_on_the_causal_variants_beats_covariates(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            panel_dir = build_cohort(root)
            self.assertEqual(self._run(root, panel_dir), 0)
            report = json.loads((root / "gate" / "GATE_RESULTS.json").read_text())
            by_label = {row["comparison"]: row for row in report["comparisons"]}
            # The annotation marks exactly the causal variants, so F must add over A.
            self.assertGreater(by_label["F_minus_A"]["delta_r2"], 0.0)
            # And PF must not be worse than F by more than sampling noise allows.
            self.assertGreater(by_label["PF_minus_F"]["ci_high"], 0.0)

    def test_the_gate_refuses_a_qualification_file_for_another_trait(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            panel_dir = build_cohort(root, n_samples=400, n_variants=40)
            tq = root / "TQ_RESULTS.json"
            tq.write_text(json.dumps({"trait_column": "n_30020_0_0",
                                      "selected_p_threshold": 0.01}))
            with self.assertRaises(SystemExit) as caught:
                c3.main([
                    "--panel-dir", str(panel_dir),
                    "--covariates", str(root / "covariates.tsv"),
                    "--covariate-column", "age",
                    "--phenotype", str(root / "phenotype.tsv"), "--trait-column", "trait",
                    "--population-features", str(root / "missing.npy"),
                    "--functional-features", str(root / "missing.npy"),
                    "--tq-results", str(tq), "--out-dir", str(root / "gate"),
                ])
            self.assertIn("n_30020_0_0", str(caught.exception))


if __name__ == "__main__":
    unittest.main(verbosity=2)
