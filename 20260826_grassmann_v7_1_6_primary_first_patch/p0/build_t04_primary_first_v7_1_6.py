from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
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


def block_rank(seed: int, values: tuple[object, ...]) -> str:
    return hashlib.sha256((str(seed) + "\t" + "\t".join(map(str, values))).encode()).hexdigest()


def make_primary_schedule(grid: dict) -> list[dict[str, object]]:
    primary = grid["primary_A1"]
    blocks = [
        (mask, fairness, data_seed, init_seed)
        for mask in primary["masks"]
        for fairness in primary["fairness"]
        for data_seed in primary["data_seeds"]
        for init_seed in primary["init_seeds"]
    ]
    seed = int(grid["schedule_seed"])
    blocks.sort(key=lambda values: block_rank(seed, values))
    gpus = list(grid["allowed_physical_gpus"])
    rows: list[dict[str, object]] = []
    for rank, (mask, fairness, data_seed, init_seed) in enumerate(blocks):
        gpu = gpus[rank % len(gpus)]
        rotation = rank % len(MODELS)
        ordered_models = MODELS[rotation:] + MODELS[:rotation]
        block_id = f"P1_B{rank + 1:03d}"
        for order, model in enumerate(ordered_models, start=1):
            rows.append({
                "run_id": f"{block_id}_O{order}_{model}_D{data_seed}_I{init_seed}",
                "block_id": block_id,
                "schedule_rank": rank + 1,
                "preferred_physical_gpu": gpu,
                "order_in_block": order,
                "model": model,
                "mask": mask,
                "fairness": fairness,
                "data_seed": data_seed,
                "init_seed": init_seed,
                "sequence_length": grid["sequence_length"],
                "steps_source": primary["steps"],
            })
    return rows


def envelope(seconds_sum: float, checkpoint_bytes_sum: int, runs_per_model: int, steps: int, margin: float) -> dict:
    raw_gpu_hours = seconds_sum * steps * runs_per_model / 3600.0
    raw_checkpoint_gib = checkpoint_bytes_sum * runs_per_model / (2**30)
    return {
        "steps_per_run": steps,
        "runs": runs_per_model * len(MODELS),
        "raw_gpu_hours": raw_gpu_hours,
        "gpu_hours_with_margin": raw_gpu_hours * margin,
        "final_checkpoint_gib_with_margin": raw_checkpoint_gib * margin,
        "wall_hours_with_margin_6_gpus": raw_gpu_hours * margin / 6.0,
        "required_window_hours_at_80pct": raw_gpu_hours * margin / 6.0 / 0.80,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Build v7.1.6 primary-first T04 contract and schedule")
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
    if grid.get("protocol_version") != "v7.1.6" or grid.get("sequence_length") != 154850:
        raise SystemExit("grid protocol or sequence length mismatch")
    if grid.get("allowed_physical_gpus") != [1, 2, 3, 4, 5, 6] or 0 not in grid.get("forbidden_physical_gpus", []):
        raise SystemExit("physical GPU policy mismatch")
    if sha256(args.profile_report) != EXPECTED_PROFILE_SHA256:
        raise SystemExit("T03 profile report SHA-256 mismatch")
    if not profile.get("valid_t03_measurement") or str(profile.get("cuda_visible_devices")) != "1":
        raise SystemExit("T03 profile is not a valid physical-GPU1 measurement")
    if material.get("status") != "PASS" or material.get("source_bcf_sha256") != EXPECTED_SOURCE_SHA256:
        raise SystemExit("T02 materialization audit mismatch")
    if material.get("final_variant_ids_sha256") != EXPECTED_VARIANT_IDS_SHA256:
        raise SystemExit("T02 variant ID hash mismatch")

    final_rows = {
        row.get("model"): row
        for row in profile.get("profiles", [])
        if row.get("sequence_length") == 154850 and row.get("status") == "PASS"
    }
    if set(final_rows) != set(MODELS):
        raise SystemExit("three final-length PASS profiles are required")
    seconds_sum = sum(float(final_rows[model]["seconds_per_step"]) for model in MODELS)
    checkpoint_sum = sum(int(final_rows[model]["checkpoint_bytes"]) for model in MODELS)
    margin = float(grid["engineering_margin_multiplier"])

    convergence = grid["convergence_pilot"]
    c0_runs_per_model = len(convergence["masks"]) * len(convergence["data_seeds"]) * len(convergence["init_seeds"])
    c0 = envelope(seconds_sum, checkpoint_sum, c0_runs_per_model, int(convergence["maximum_steps"]), margin)

    primary = grid["primary_A1"]
    p1_runs_per_model = len(primary["masks"]) * len(primary["fairness"]) * len(primary["data_seeds"]) * len(primary["init_seeds"])
    if p1_runs_per_model * len(MODELS) != 120:
        raise SystemExit("primary schedule is not 120 runs")
    primary_scenarios = {
        str(steps): envelope(seconds_sum, checkpoint_sum, p1_runs_per_model, steps, margin)
        for steps in convergence["candidate_common_steps"]
    }
    full_scenarios = {
        str(steps): envelope(seconds_sum, checkpoint_sum, 120, steps, margin)
        for steps in convergence["candidate_common_steps"]
    }

    worst_primary = primary_scenarios[str(max(convergence["candidate_common_steps"]))]
    cap = grid["primary_capacity_envelope"]
    usage = float(grid["capacity_usage_limit"])
    required_storage = (
        worst_primary["final_checkpoint_gib_with_margin"]
        + float(cap["synthetic_corpus_storage_cap_gib"])
        + float(cap["logs_and_metadata_storage_cap_gib"])
    )
    checks = {
        "primary_gpu_hours_within_80pct": worst_primary["gpu_hours_with_margin"] <= usage * float(cap["gpu_hours"]),
        "primary_wall_time_within_80pct": worst_primary["wall_hours_with_margin_6_gpus"] <= usage * float(cap["available_window_hours"]),
        "primary_storage_within_80pct": required_storage <= usage * float(cap["storage_gib"]),
        "six_nonzero_gpus": cap["gpu_count"] == 6 and cap["concurrency_limit"] == 6 and 0 not in grid["allowed_physical_gpus"],
        "hgdp_forbidden_in_c0": convergence["hgdp_access"] == "FORBIDDEN",
    }

    schedule = make_primary_schedule(grid)
    gpu_block_counts = Counter(row["preferred_physical_gpu"] for row in schedule if row["order_in_block"] == 1)
    if set(gpu_block_counts) != set(grid["allowed_physical_gpus"]):
        raise SystemExit("not all allowed GPUs receive blocks")
    if max(gpu_block_counts.values()) - min(gpu_block_counts.values()) > 1:
        raise SystemExit("GPU block allocation is imbalanced")

    args.output_dir.mkdir(parents=True, exist_ok=False)
    schedule_path = args.output_dir / "PRIMARY_RUN_SCHEDULE.v7.1.6.tsv"
    fields = list(schedule[0])
    with schedule_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(schedule)

    contract = {
        "schema_version": "1.0",
        "protocol_version": "v7.1.6",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "T04_PRIMARY_FIRST_CONDITIONAL_GO" if all(checks.values()) else "T04_NO_GO",
        "condition": "all six scheduled GPUs must be idle and locked; otherwise wait",
        "T02_status": "GO",
        "T03_status": "GO",
        "profile_report_sha256": sha256(args.profile_report),
        "materialization_audit_sha256": sha256(args.materialization_audit),
        "sequence_length": 154850,
        "convergence_pilot": c0,
        "convergence_rule": {
            "candidate_common_steps": convergence["candidate_common_steps"],
            "maximum_nll_improvement": convergence["maximum_nll_improvement_for_convergence"],
            "selection_blinding": convergence["budget_selection_blinding"],
            "per_run_early_stopping": "FORBIDDEN",
        },
        "primary_A1": {
            "runs": len(schedule),
            "scenarios": primary_scenarios,
            "capacity_required_storage_gib_at_10000": required_storage,
            "checks": checks,
            "schedule_sha256": sha256(schedule_path),
            "gpu_block_counts": {str(k): v for k, v in sorted(gpu_block_counts.items())},
        },
        "full_360_run_scenarios": full_scenarios,
        "full_360_status": "DEFERRED_NOT_REQUIRED_FOR_PRIMARY_DECISION",
        "allowed_physical_gpus": grid["allowed_physical_gpus"],
        "physical_gpu_0": "FORBIDDEN",
    }
    contract_path = args.output_dir / "COMPUTE_CONTRACT.v7.1.6.json"
    contract_path.write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    manifest_path = args.output_dir / "T04_PRIMARY_FIRST.v7.1.6.sha256"
    with manifest_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(f"{sha256(contract_path)}  {contract_path.name}\n")
        handle.write(f"{sha256(schedule_path)}  {schedule_path.name}\n")
    print(json.dumps(contract, indent=2, sort_keys=True))
    if contract["status"] == "T04_NO_GO":
        raise SystemExit(4)


if __name__ == "__main__":
    main()
