#!/usr/bin/env python3
"""Synthetic tests for the R1 inventory and R2 family-axis tools.

No real cohort data is touched. The fixtures are tiny hand-built filesets that
exercise detection, identifier hygiene, the read-only boundary, component
construction and deliberate leakage failure.
"""

from __future__ import annotations

import io
import json
import math
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

import r1_bind_contract  # noqa: E402
import r1_inventory  # noqa: E402
import r1_validate_inventory  # noqa: E402
import r2_build_family_manifest  # noqa: E402
import r2_validate_family_manifest  # noqa: E402

SAMPLE_COUNT = 60
VARIANT_COUNT = 400


def write_plink_fileset(directory: Path, stem: str) -> None:
    fam_lines = []
    for index in range(SAMPLE_COUNT):
        fam_lines.append(f"FID{index} IID{index} 0 0 {1 + index % 2} -9")
    (directory / f"{stem}.fam").write_text("\n".join(fam_lines) + "\n", encoding="utf-8")

    bim_lines = []
    for index in range(VARIANT_COUNT):
        position = 100_000 + index * 600_000
        bim_lines.append(f"1 rs{index} 0 {position} A G")
    # One marker beyond the GRCh38 chr1 length but inside GRCh37's.
    bim_lines.append("1 rs_tail 0 249000000 C T")
    (directory / f"{stem}.bim").write_text("\n".join(bim_lines) + "\n", encoding="utf-8")

    variants = VARIANT_COUNT + 1
    payload = b"\x6c\x1b\x01" + b"\x00" * (math.ceil(SAMPLE_COUNT / 4) * variants)
    (directory / f"{stem}.bed").write_bytes(payload)


def write_kinship(path: Path, trio: bool = True) -> None:
    rows = ["FID1 ID1 FID2 ID2 N_SNP HetHet IBS0 Kinship"]
    if trio:
        rows.append("0 IID0 0 IID1 5000 0.2 0.0 0.2500")
        rows.append("0 IID1 0 IID2 5000 0.2 0.0 0.2500")
        rows.append("0 IID3 0 IID4 5000 0.1 0.0 0.0600")
    rows.append("0 IID5 0 IID6 5000 0.05 0.01 0.0100")
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def write_genetic_map(path: Path) -> None:
    rows = ["position COMBINED_rate(cM/Mb) Genetic_Map(cM)"]
    for index in range(50):
        rows.append(f"{100000 + index * 500000} 1.0 {index * 0.5}")
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def run_json(callable_, argv: list[str]) -> dict:
    buffer = io.StringIO()
    with redirect_stdout(buffer):
        code = callable_(argv)
    return {"exit_code": code, "stdout": json.loads(buffer.getvalue())}


class FixtureCase(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        base = Path(self.temporary.name)
        self.root = base / "cohort"
        self.root.mkdir()
        self.out = base / "out"
        self.out.mkdir()
        write_plink_fileset(self.root, "ukb_cal_chr1_v2")
        write_kinship(self.root / "ukb_rel_a11111.dat")
        write_genetic_map(self.root / "genetic_map_chr1.txt")

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def build_inventory(self, extra: list[str] | None = None) -> dict:
        argv = ["--root", str(self.root), "--out-dir", str(self.out),
                "--cohort-id", "TEST_COHORT", "--authorization-id", "TEST_AUTH",
                "--array-variant-min", "100"]
        result = run_json(r1_inventory.main, argv + (extra or []))
        self.assertEqual(result["exit_code"], 0)
        return json.loads((self.out / "INVENTORY.json").read_text(encoding="utf-8"))


class InventoryTests(FixtureCase):
    def test_detects_fileset_counts_and_build_hint(self) -> None:
        payload = self.build_inventory()
        filesets = payload["genotype_filesets"]
        self.assertEqual(len(filesets), 1)
        group = filesets[0]
        self.assertEqual(group["format"], "plink1_bed")
        self.assertEqual(group["sample_count"], SAMPLE_COUNT)
        self.assertEqual(group["variant_count"], VARIANT_COUNT + 1)
        self.assertEqual(group["build_hint"]["hint"], "GRCh37")
        self.assertTrue(group["e0_eligibility"]["e0_primary_target_eligible"])

    def test_kinship_components_are_summarized_without_identifiers(self) -> None:
        payload = self.build_inventory()
        kinship = payload["kinship_files"]
        self.assertEqual(len(kinship), 1)
        parsed = kinship[0]["parsed"]
        self.assertTrue(parsed["parsed"])
        # IID0-IID1-IID2 form one component; IID3-IID4 another; the 0.01 pair is below threshold.
        components = parsed["components_at_0.0442"]
        self.assertEqual(components["multi_member_component_count"], 2)
        self.assertEqual(components["largest_component_size"], 3)
        self.assertEqual(parsed["pairs_at_or_above"]["third_degree_or_closer_0.0442"], 3)
        self.assertEqual(parsed["pairs_at_or_above"]["first_degree_or_closer_0.177"], 2)
        serialized = json.dumps(payload)
        self.assertNotIn("IID0", serialized)
        self.assertNotIn("rs17", serialized)

    def test_refuses_output_inside_scanned_root(self) -> None:
        with self.assertRaises(SystemExit):
            r1_inventory.main(
                ["--root", str(self.root), "--out-dir", str(self.root / "inside")]
            )

    def test_does_not_write_into_root(self) -> None:
        before = sorted(item.name for item in self.root.iterdir())
        self.build_inventory()
        after = sorted(item.name for item in self.root.iterdir())
        self.assertEqual(before, after)

    def test_validator_passes_and_flags_environment_warnings(self) -> None:
        self.build_inventory()
        result = run_json(
            r1_validate_inventory.main,
            ["--inventory", str(self.out / "INVENTORY.json")],
        )
        self.assertEqual(result["exit_code"], 0)
        report = result["stdout"]
        self.assertEqual(report["status"], "PASS")
        self.assertTrue(report["r2_allowed"])
        self.assertFalse(report["run_authorized"])
        self.assertEqual(report["filesets"]["e0_primary_eligible_count"], 1)

    def test_validator_rejects_identifier_leak(self) -> None:
        self.build_inventory()
        path = self.out / "INVENTORY.json"
        payload = json.loads(path.read_text(encoding="utf-8"))
        payload["kinship_files"][0]["parsed"]["members"] = ["IID0", "IID1", "IID2"]
        path.write_text(json.dumps(payload), encoding="utf-8")
        result = run_json(r1_validate_inventory.main, ["--inventory", str(path)])
        self.assertEqual(result["exit_code"], 1)
        self.assertEqual(result["stdout"]["status"], "FAIL")
        self.assertFalse(result["stdout"]["r2_allowed"])

    def test_imputed_scale_fileset_is_not_e0_eligible(self) -> None:
        group = {"fileset": "/data/ukb_imp_chr1_v3", "variant_count": 9_000_000}
        eligibility = r1_inventory.e0_eligibility(group)
        self.assertEqual(eligibility["role"], "imputed_candidate")
        self.assertFalse(eligibility["e0_primary_target_eligible"])

    def test_contract_draft_keeps_unbound_fields_null(self) -> None:
        self.build_inventory()
        package = Path(__file__).resolve().parents[1]
        template = package / "CONTRACT.template.json"
        if not template.exists():
            template = package.parent / "CONTRACT.template.json"
        output = self.out / "CONTRACT.R1_DRAFT.json"
        result = run_json(
            r1_bind_contract.main,
            [
                "--inventory", str(self.out / "INVENTORY.json"),
                "--template", str(template),
                "--output", str(output),
            ],
        )
        self.assertEqual(result["exit_code"], 0)
        draft = json.loads(output.read_text(encoding="utf-8"))
        self.assertIsNone(draft["data"]["genome_build"])
        self.assertFalse(draft["authorization"]["run_authorized"])
        self.assertIn("data.genome_build", draft["r1_unbound_fields"])
        self.assertEqual(draft["data"]["genotype_format"], "plink1_bed")


class FamilyManifestTests(FixtureCase):
    def build_family(self, out_dir: Path, seed: int = 20260916) -> dict:
        argv = [
            "--sample-source", str(self.root / "ukb_cal_chr1_v2.fam"),
            "--kinship", str(self.root / "ukb_rel_a11111.dat"),
            "--seed", str(seed),
            "--out-dir", str(out_dir),
        ]
        result = run_json(r2_build_family_manifest.main, argv)
        self.assertEqual(result["exit_code"], 0)
        return result["stdout"]

    def test_components_and_splits(self) -> None:
        out_dir = self.out / "family"
        summary = self.build_family(out_dir)
        self.assertEqual(summary["samples"], SAMPLE_COUNT)
        # 60 samples, one trio and one pair collapse into components: 60 - 2 - 1 = 57.
        self.assertEqual(summary["components"], 57)
        counts = summary["split_sample_counts"]
        self.assertEqual(sum(counts.values()), SAMPLE_COUNT)
        self.assertGreater(counts["train"], counts["validation"])

    def test_deterministic_for_same_seed(self) -> None:
        first = self.build_family(self.out / "a")
        second = self.build_family(self.out / "b")
        self.assertEqual(first["family_manifest_sha256"], second["family_manifest_sha256"])

    def test_exclusions_are_dropped(self) -> None:
        exclusion = self.out / "withdrawn.txt"
        exclusion.write_text("IID7\nIID8\n", encoding="utf-8")
        argv = [
            "--sample-source", str(self.root / "ukb_cal_chr1_v2.fam"),
            "--kinship", str(self.root / "ukb_rel_a11111.dat"),
            "--exclude", str(exclusion),
            "--seed", "1",
            "--out-dir", str(self.out / "excluded"),
        ]
        result = run_json(r2_build_family_manifest.main, argv)
        self.assertEqual(result["stdout"]["samples"], SAMPLE_COUNT - 2)
        manifest = (self.out / "excluded" / "family_manifest.tsv").read_text(encoding="utf-8")
        self.assertNotIn("\tIID7\t", manifest)
        self.assertNotIn("IID7\tFC_", manifest)

    def test_validator_passes_on_generated_manifest(self) -> None:
        out_dir = self.out / "family"
        self.build_family(out_dir)
        result = run_json(
            r2_validate_family_manifest.main,
            [
                "--family-manifest", str(out_dir / "family_manifest.tsv"),
                "--kinship", str(self.root / "ukb_rel_a11111.dat"),
            ],
        )
        self.assertEqual(result["exit_code"], 0)
        report = result["stdout"]
        self.assertEqual(report["status"], "PASS")
        self.assertEqual(report["component_cross_split_count"], 0)
        self.assertEqual(report["related_pairs_crossing_splits"], 0)
        self.assertTrue(report["partition_check"]["identical"])

    def test_validator_fails_on_deliberate_family_leak(self) -> None:
        out_dir = self.out / "family"
        self.build_family(out_dir)
        path = out_dir / "family_manifest.tsv"
        lines = path.read_text(encoding="utf-8").splitlines()
        header, rows = lines[0], lines[1:]
        patched = []
        moved = False
        target_component = None
        for row in rows:
            fields = row.split("\t")
            if not moved and int(fields[2]) == 3:
                if target_component is None:
                    target_component = fields[1]
                if fields[1] == target_component:
                    fields[3] = "test" if fields[3] != "test" else "train"
                    moved = True
            patched.append("\t".join(fields))
        self.assertTrue(moved, "fixture must contain a multi-member component")
        path.write_text("\n".join([header, *patched]) + "\n", encoding="utf-8")
        result = run_json(
            r2_validate_family_manifest.main,
            [
                "--family-manifest", str(path),
                "--kinship", str(self.root / "ukb_rel_a11111.dat"),
            ],
        )
        self.assertEqual(result["exit_code"], 1)
        self.assertEqual(result["stdout"]["status"], "FAIL")
        self.assertGreater(result["stdout"]["component_cross_split_count"], 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
