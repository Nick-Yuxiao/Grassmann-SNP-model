from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


MODELS = (
    "local_attn_8m_w256",
    "local_attn_gpc_8m_w256",
    "grassmann_full_8m_w256",
)
EXPECTED_PROFILE_SHA256 = "99c897000ec458b7099df1bc2a0849a1f8a3ac8aca5e3ec56d9a7039eadb0a09"
EXPECTED_SOURCE_SHA256 = "09e50c698522d19bbc2c43a2ea2d29b290d63cf8c1d6983c5198dcb4289ff3fa"
EXPECTED_VARIANT_IDS_SHA256 = "69fb8133b63e995c64b79f33d339214732dd6bf9d4cdf0515c97efbe5640d0bd"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description="Build v7.1.5 two-tier T04 compute contract")
    parser.add_argument("--grid", type=Path, required=True)
    parser.add_argument("--profile-report", type=Path, required=True)
    parser.add_argument("--materialization-audit", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    if args.output_dir.exists():
        raise SystemExit(f"refusing to overwrite output directory: {args.output_dir}")
    grid = json.loads(args.grid.read_text(encoding="utf-8"))
    profile = json.loads(args.profile_report.read_text(encoding="utf-8"))
    material = json.loads(args.materialization_audit.read_text(encoding="utf-8"))

    if grid.get("protocol_version") != "v7.1.5" or grid.get("sequence_length") != 154850:
        raise SystemExit("pilot grid protocol or sequence length mismatch")
    if sha256(args.profile_report) != EXPECTED_PROFILE_SHA256:
        raise SystemExit("T03 profile report SHA-256 mismatch")
    if not profile.get("valid_t03_measurement") or str(profile.get("cuda_visible_devices")) != "1":
        raise SystemExit("T03 is not a valid physical-GPU1 measurement")
    if material.get("status") != "PASS" or material.get("source_bcf_sha256") != EXPECTED_SOURCE_SHA256:
        raise SystemExit("materialization audit is not PASS or source hash differs")
    if material.get("final_variant_ids_sha256") != EXPECTED_VARIANT_IDS_SHA256:
        raise SystemExit("materialized variant IDs differ from frozen v7.1.3 IDs")
    if material.get("info_ids") != [] or material.get("format_ids") != ["GT"]:
        raise SystemExit("materialized annotation leakage contract failed")

    final_rows = {
        row.get("model"): row
        for row in profile.get("profiles", [])
        if row.get("sequence_length") == 154850 and row.get("status") == "PASS"
    }
    if set(final_rows) != set(MODELS):
        raise SystemExit(f"missing final-length PASS profiles: {sorted(final_rows)}")

    masks = len(grid["mask_cells"])
    data_seeds = len(grid["pilot_data_seeds"])
    init_seeds = len(grid["pilot_init_seeds"])
    pilot_steps = int(grid["train_steps_per_run"])
    pilot_runs_per_model = masks * data_seeds * init_seeds
    margin = float(grid["engineering_margin_multiplier"])
    raw_pilot_gpu_hours = sum(
        float(final_rows[model]["seconds_per_step"]) * pilot_steps * pilot_runs_per_model
        for model in MODELS
    ) / 3600.0
    raw_pilot_storage_gib = sum(
        int(final_rows[model]["checkpoint_bytes"]) * pilot_runs_per_model
        for model in MODELS
    ) / (2**30)
    pilot_gpu_hours = raw_pilot_gpu_hours * margin
    pilot_storage_gib = raw_pilot_storage_gib * margin

    envelope = grid["pilot_capacity_envelope"]
    usage = float(envelope["capacity_usage_limit"])
    concurrency = int(envelope["concurrency_limit"])
    pilot_wall_hours = pilot_gpu_hours / concurrency
    pilot_checks = {
        "gpu_hours_within_envelope": pilot_gpu_hours <= usage * float(envelope["gpu_hours"]),
        "storage_within_envelope": pilot_storage_gib <= usage * float(envelope["storage_gib"]),
        "wall_time_within_envelope": pilot_wall_hours <= usage * float(envelope["available_window_hours"]),
        "physical_gpu_is_not_zero": int(envelope["physical_gpu"]) != 0,
        "hgdp_is_forbidden": grid["hgdp_access"] == "FORBIDDEN",
    }
    pilot_status = "PILOT_CAPACITY_GO_IDLE_ONLY" if all(pilot_checks.values()) else "PILOT_CAPACITY_NO_GO"

    full = grid["full_a1_scenario"]
    full_runs_per_model = int(full["cells_per_model"]) * int(full["runs_per_cell"])
    full_steps = int(full["train_steps_per_run"])
    raw_full_gpu_hours = sum(
        float(final_rows[model]["seconds_per_step"]) * full_steps * full_runs_per_model
        for model in MODELS
    ) / 3600.0
    raw_full_storage_gib = sum(
        int(final_rows[model]["checkpoint_bytes"]) * full_runs_per_model
        for model in MODELS
    ) / (2**30)
    full_gpu_hours = raw_full_gpu_hours * margin
    full_storage_gib = raw_full_storage_gib * margin

    contract = {
        "schema_version": "1.0",
        "protocol_version": "v7.1.5",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "T04_PILOT_GO_FULL_A1_UNFUNDED" if pilot_status.startswith("PILOT_CAPACITY_GO") else "T04_NO_GO",
        "T02": {
            "status": "GO",
            "materialization_audit_sha256": sha256(args.materialization_audit),
            "sequence_length": 154850,
        },
        "T03": {
            "status": "GO",
            "profile_report_sha256": sha256(args.profile_report),
            "physical_gpu": 1,
            "all_final_length_cells_pass": True,
            "profiles": {
                model: {
                    "seconds_per_step": final_rows[model]["seconds_per_step"],
                    "peak_allocated_gib": final_rows[model]["peak_allocated_bytes"] / (2**30),
                    "peak_reserved_gib": final_rows[model]["peak_reserved_bytes"] / (2**30),
                    "checkpoint_mib": final_rows[model]["checkpoint_bytes"] / (2**20),
                }
                for model in MODELS
            },
        },
        "pilot": {
            "status": pilot_status,
            "runs": pilot_runs_per_model * len(MODELS),
            "runs_per_model": pilot_runs_per_model,
            "train_steps_per_run": pilot_steps,
            "raw_gpu_hours": raw_pilot_gpu_hours,
            "gpu_hours_with_margin": pilot_gpu_hours,
            "storage_gib_with_margin": pilot_storage_gib,
            "wall_hours_with_margin": pilot_wall_hours,
            "checks": pilot_checks,
            "inference": "EXPLORATORY_NO_CONFIRMATORY_CLAIM",
            "hgdp_access": "FORBIDDEN",
        },
        "full_A1_scenario": {
            "status": "UNFUNDED_UNRESERVED",
            "runs": full_runs_per_model * len(MODELS),
            "train_steps_per_run": full_steps,
            "raw_gpu_hours": raw_full_gpu_hours,
            "gpu_hours_with_margin": full_gpu_hours,
            "storage_gib_with_margin_final_checkpoint_only": full_storage_gib,
            "wall_hours_with_margin_1_gpu": full_gpu_hours,
            "wall_hours_with_margin_6_gpus": full_gpu_hours / 6.0,
            "excluded_from_storage_estimate": ["synthetic_corpus", "intermediate_checkpoints"],
        },
        "pilot_grid_sha256": sha256(args.grid),
    }
    args.output_dir.mkdir(parents=True, exist_ok=False)
    output = args.output_dir / "COMPUTE_CONTRACT.v7.1.5.json"
    output.write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    manifest = args.output_dir / "T04_CONTRACT.v7.1.5.sha256"
    manifest.write_text(f"{sha256(output)}  {output.name}\n", encoding="utf-8")
    print(json.dumps(contract, indent=2, sort_keys=True))
    if contract["status"] == "T04_NO_GO":
        raise SystemExit(4)


if __name__ == "__main__":
    main()
