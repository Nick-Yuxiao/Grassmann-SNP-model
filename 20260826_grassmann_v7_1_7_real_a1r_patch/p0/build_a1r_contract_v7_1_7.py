from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


MODELS = ("local_attn_8m_w256", "local_attn_gpc_8m_w256", "grassmann_full_8m_w256")
EXPECTED_PROFILE_SHA256 = "99c897000ec458b7099df1bc2a0849a1f8a3ac8aca5e3ec56d9a7039eadb0a09"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def rank(seed: int, *parts: object) -> str:
    return hashlib.sha256((str(seed) + "\t" + "\t".join(map(str, parts))).encode()).hexdigest()


def model_rows(block: dict[str, object], rotation: int, grid: dict) -> list[dict[str, object]]:
    ordered = MODELS[rotation:] + MODELS[:rotation]
    rows = []
    for order, model in enumerate(ordered, 1):
        row = dict(block)
        row.update({
            "run_id": f"{block['block_id']}_O{order}_{model}",
            "order_in_block": order,
            "model": model,
            "sequence_length": grid["sequence_length"],
            "steps_source": "C0_SELECTED_COMMON_K",
        })
        rows.append(row)
    return rows


def build_schedule(grid: dict) -> list[dict[str, object]]:
    p = grid["primary_A1R_100pct"]
    seed = int(grid["schedule_seed"])
    keys = [(m, f, s, i) for m in p["masks"] for f in p["fairness"] for s in p["mask_seeds"] for i in p["init_seeds"]]
    keys.sort(key=lambda x: rank(seed, "P1", *x))
    gpus = grid["allowed_physical_gpus"]
    rows: list[dict[str, object]] = []
    anchors: dict[tuple[str, int], int] = {}
    for n, (mask, fairness, mask_seed, init_seed) in enumerate(keys, 1):
        gpu = gpus[(n - 1) % len(gpus)]
        block = {
            "stage": "PRIMARY_100P",
            "block_id": f"P1_B{n:03d}",
            "preferred_physical_gpu": gpu,
            "sample_fraction": "1.00",
            "sample_count": 2247,
            "mask": mask,
            "fairness": fairness,
            "mask_seed": mask_seed,
            "init_seed": init_seed,
            "size_curve_id": "",
        }
        if fairness == "matched_compute" and mask_seed in (1, 2) and init_seed == 1:
            curve = f"S_{mask}_M{mask_seed}_I1"
            block["size_curve_id"] = curve
            anchors[(mask, mask_seed)] = gpu
        rows.extend(model_rows(block, (n - 1) % 3, grid))

    d = grid["sample_size_diagnostic"]
    diagnostic = []
    dkeys = [(fraction, count, mask, mask_seed) for fraction, count in zip(d["sample_fractions"], d["sample_counts"]) for mask in d["masks"] for mask_seed in d["mask_seeds"]]
    dkeys.sort(key=lambda x: rank(seed, "D1", *x))
    for n, (fraction, count, mask, mask_seed) in enumerate(dkeys, 1):
        gpu = anchors[(mask, mask_seed)]
        block = {
            "stage": "DIAGNOSTIC_SIZE",
            "block_id": f"D1_B{n:03d}",
            "preferred_physical_gpu": gpu,
            "sample_fraction": f"{fraction:.2f}",
            "sample_count": count,
            "mask": mask,
            "fairness": "matched_compute",
            "mask_seed": mask_seed,
            "init_seed": 1,
            "size_curve_id": f"S_{mask}_M{mask_seed}_I1",
        }
        diagnostic.extend(model_rows(block, (n - 1) % 3, grid))
    return rows + diagnostic


def build_pilot(grid: dict) -> list[dict[str, object]]:
    c = grid["convergence_pilot"]
    rows = []
    n = 0
    for mask in c["masks"]:
        for mask_seed in c["mask_seeds"]:
            n += 1
            gpu = grid["allowed_physical_gpus"][(n - 1) % 6]
            block = {
                "stage": "C0_CONVERGENCE",
                "block_id": f"C0_B{n:03d}",
                "preferred_physical_gpu": gpu,
                "sample_fraction": "1.00",
                "sample_count": 2247,
                "mask": mask,
                "fairness": "pilot_profile",
                "mask_seed": mask_seed,
                "init_seed": c["init_seeds"][0],
                "size_curve_id": "",
            }
            rows.extend(model_rows(block, (n - 1) % 3, grid))
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Build v7.1.7 real-1KGP A1-R contract")
    parser.add_argument("--grid", type=Path, required=True)
    parser.add_argument("--profile-report", type=Path, required=True)
    parser.add_argument("--subset-freeze", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise SystemExit(f"refusing to overwrite output directory: {args.output_dir}")
    grid = json.loads(args.grid.read_text(encoding="utf-8"))
    profile = json.loads(args.profile_report.read_text(encoding="utf-8"))
    subset = json.loads(args.subset_freeze.read_text(encoding="utf-8"))
    if grid.get("protocol_version") != "v7.1.7" or grid.get("data_contract", {}).get("hgdp_access") != "FORBIDDEN":
        raise SystemExit("invalid v7.1.7 grid")
    if sha256(args.profile_report) != EXPECTED_PROFILE_SHA256 or not profile.get("valid_t03_measurement"):
        raise SystemExit("T03 profile mismatch")
    if subset.get("status") != "PASS" or subset.get("counts") != {"100P_2247": 2247, "25P_562": 562, "50P_1124": 1124}:
        raise SystemExit("subset freeze mismatch")
    schedule = build_schedule(grid)
    pilot = build_pilot(grid)
    if len(schedule) != 144 or len(pilot) != 12:
        raise SystemExit("run-count mismatch")
    if any(int(r["preferred_physical_gpu"]) == 0 for r in schedule + pilot):
        raise SystemExit("GPU0 scheduled")
    primary = [r for r in schedule if r["stage"] == "PRIMARY_100P"]
    diagnostic = [r for r in schedule if r["stage"] == "DIAGNOSTIC_SIZE"]
    if len(primary) != 120 or len(diagnostic) != 24:
        raise SystemExit("stage-count mismatch")
    for curve in {r["size_curve_id"] for r in diagnostic}:
        gpus = {int(r["preferred_physical_gpu"]) for r in schedule if r["size_curve_id"] == curve}
        if len(gpus) != 1:
            raise SystemExit(f"size curve crosses GPUs: {curve}")

    args.output_dir.mkdir(parents=True, exist_ok=False)
    fields = list(schedule[0])
    paths = []
    for name, values in (("A1R_RUN_SCHEDULE.v7.1.7.tsv", schedule), ("A1R_C0_SCHEDULE.v7.1.7.tsv", pilot)):
        path = args.output_dir / name
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(values)
        paths.append(path)
    final_profiles = [r for r in profile["profiles"] if r.get("sequence_length") == 154850 and r.get("status") == "PASS"]
    seconds = {r["model"]: float(r["seconds_per_step"]) for r in final_profiles}
    margin = float(grid["engineering_margin_multiplier"])
    raw_pilot_hours_at_10000 = sum(seconds[r["model"]] for r in pilot) * 10000 / 3600
    scenarios = {}
    for steps in grid["convergence_pilot"]["candidate_common_steps"]:
        raw_hours = sum(seconds[r["model"]] for r in schedule) * steps / 3600
        scenarios[str(steps)] = {
            "raw_gpu_hours": raw_hours,
            "gpu_hours_with_engineering_margin": raw_hours * margin,
            "wall_hours_at_6_gpus_with_engineering_margin": raw_hours * margin / 6,
        }
    contract = {
        "schema_version": "1.0",
        "protocol_version": "v7.1.7",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "A1R_CONTRACT_READY",
        "pilot_runs": len(pilot),
        "post_pilot_runs": len(schedule),
        "primary_runs": len(primary),
        "diagnostic_runs": len(diagnostic),
        "engineering_margin_multiplier": margin,
        "post_pilot_scenarios": scenarios,
        "pilot_10000_step_gpu_hours_with_engineering_margin": raw_pilot_hours_at_10000 * margin,
        "pilot_10000_step_wall_hours_at_6_gpus_with_engineering_margin": raw_pilot_hours_at_10000 * margin / 6,
        "allowed_physical_gpus": grid["allowed_physical_gpus"],
        "physical_gpu_0": "FORBIDDEN",
        "hgdp_access": "FORBIDDEN",
        "subset_freeze_sha256": sha256(args.subset_freeze),
        "profile_report_sha256": sha256(args.profile_report),
        "schedule_sha256": sha256(paths[0]),
        "pilot_schedule_sha256": sha256(paths[1]),
        "gpu_primary_block_counts": dict(sorted(Counter(int(r["preferred_physical_gpu"]) for r in primary if r["order_in_block"] == 1).items())),
        "negative_default": "INCONCLUSIVE_SAMPLE_LIMITED_HAPNEST",
    }
    contract_path = args.output_dir / "A1R_COMPUTE_CONTRACT.v7.1.7.json"
    contract_path.write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    paths.append(contract_path)
    manifest = args.output_dir / "A1R_T04.v7.1.7.sha256"
    manifest.write_text("".join(f"{sha256(p)}  {p.name}\n" for p in paths), encoding="utf-8")
    print(json.dumps(contract, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
