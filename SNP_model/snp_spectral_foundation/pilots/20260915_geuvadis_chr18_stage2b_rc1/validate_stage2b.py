from __future__ import annotations

import csv
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path


PILOT = Path(__file__).resolve().parent
ROOT = PILOT.parents[2]
RESULTS = PILOT / "results"
ASSET = ROOT / "snp_spectral_foundation" / "pilot_assets" / "geuvadis_tensorqtl_chr18"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    checks: dict[str, object] = {}
    final = json.loads((RESULTS / "FINAL_STATUS.json").read_text(encoding="utf-8"))
    task = json.loads((RESULTS / "TASK_GATE_RESULTS.json").read_text(encoding="utf-8"))
    sensitivity = json.loads((RESULTS / "DESIGN_SENSITIVITY.json").read_text(encoding="utf-8"))
    traits = json.loads((RESULTS / "TRAIT_PANEL.json").read_text(encoding="utf-8"))
    models = json.loads((RESULTS / "FROZEN_MODELS.json").read_text(encoding="utf-8"))
    binding = json.loads((RESULTS / "RUN_BINDING.json").read_text(encoding="utf-8"))
    with (PILOT / "SAMPLE_MANIFEST.tsv").open("r", encoding="utf-8", newline="") as handle:
        samples = list(csv.DictReader(handle, delimiter="\t"))

    assert final["status"] == "TASK-INELIGIBLE"
    assert final["bridge_gate_opened"] is False and final["test_open_count"] == 0
    assert task["pass"] is False and task["status"] == "TASK-INELIGIBLE"
    assert not (RESULTS / "BRIDGE_GATE_RESULTS.json").exists()
    checks["bridge_outcome_artifact_absent"] = True

    counts = Counter(row["primary_role"] for row in samples)
    assert counts == {"development": 150, "task_gate": 50, "bridge_test": 230, "EXCLUDED_RC2_TEST": 15}
    families = {
        role: {row["family"] for row in samples if row["primary_role"] == role}
        for role in ["development", "task_gate", "bridge_test"]
    }
    assert not families["development"] & families["task_gate"]
    assert not families["development"] & families["bridge_test"]
    assert not families["task_gate"] & families["bridge_test"]
    checks["sample_counts"] = dict(counts)
    checks["family_overlap_zero"] = True

    assert len(traits) == len(models) == 8
    for trait, model in zip(traits, models, strict=True):
        assert trait["gene_id"] == model["gene_id"]
        assert len(trait["variant_ids"]) == len(set(trait["variant_ids"])) == 256
        assert trait["positions_1based"] == sorted(trait["positions_1based"])
        assert all(abs(int(pos) - int(trait["tss_1based"])) <= 1_000_000 for pos in trait["positions_1based"])
        checkpoint = RESULTS / "checkpoints" / f"{trait['gene_id']}.pt"
        assert sha256_file(checkpoint) == model["checkpoint_sha256"]
        assert model["pretrain_validation"]["contextual_lift"] > 0
    checks["trait_and_panel_contract"] = True
    checks["all_encoder_contextual_lifts_positive"] = True

    macro_a = sum(float(row["r2_A"]) for row in task["trait_metrics"]) / 8
    macro_b = sum(float(row["r2_B"]) for row in task["trait_metrics"]) / 8
    assert abs(macro_a - float(task["macro_r2"]["A"])) < 1e-12
    assert abs(macro_b - float(task["macro_r2"]["B"])) < 1e-12
    assert abs((macro_b - macro_a) - float(task["macro_delta_r2"])) < 1e-12
    assert float(task["macro_delta_r2"]) > 0 and float(task["paired_bootstrap_ci95"][0]) <= 0
    checks["task_point_estimate_recomputed"] = True
    checks["task_gate_failure_reason"] = "paired_ci95_lower_not_gt_zero"
    assert sensitivity["eligible"] is True and float(sensitivity["estimated_power"]) >= 0.80
    checks["design_sensitivity_passed_before_task"] = True

    for filename, expected in binding["assets"].items():
        path = ASSET / filename
        assert path.stat().st_size == int(expected["bytes"])
        assert sha256_file(path) == expected["sha256"]
    checks["asset_binding_verified"] = True

    expression = ASSET / "GEUVADIS.445_samples.expression.bed.gz"
    uncompressed = 0
    with gzip.open(expression, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            uncompressed += len(block)
    assert uncompressed == 173_950_222
    checks["expression_gzip_crc_and_uncompressed_bytes"] = uncompressed

    old_manifest = RESULTS / "RUN_MANIFEST.sha256"
    for line in old_manifest.read_text(encoding="utf-8").splitlines():
        expected, relative = line.split("  ", 1)
        assert sha256_file(ROOT / Path(relative)) == expected
    checks["original_run_manifest_verified"] = True

    validation = {
        "validation_status": "PASS",
        "scientific_status": final["status"],
        "checks": checks,
        "important_boundary": "Task Gate failed; Bridge phenotype was not evaluated, so no C-B conclusion exists for Stage 2B rc1.",
    }
    validation_path = RESULTS / "VALIDATION.json"
    validation_path.write_text(json.dumps(validation, indent=2) + "\n", encoding="utf-8")
    fixed_inputs = [
        Path(__file__),
        PILOT / "run_stage2b.py",
        PILOT / "freeze_design.py",
        PILOT / "audit_pgen.py",
        PILOT / "FROZEN_PROTOCOL.zh-CN.md",
        PILOT / "CONFIG_FROZEN.json",
        PILOT / "DESIGN_FREEZE.json",
        PILOT / "SAMPLE_MANIFEST.tsv",
        PILOT / "IMPLEMENTATION_INCIDENT_01.md",
    ]
    result_files = sorted(path for path in RESULTS.rglob("*") if path.is_file() and path.name != "RUN_MANIFEST.sha256")
    lines = [f"{sha256_file(path)}  {path.relative_to(ROOT).as_posix()}" for path in fixed_inputs + result_files]
    old_manifest.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(validation, indent=2))


if __name__ == "__main__":
    main()
