from __future__ import annotations

import argparse
import json
import statistics
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


LRS = [0.0001, 0.0002, 0.0004, 0.0008]
MODELS = ["local_attn_8m_w256", "local_attn_gpc_8m_w256", "grassmann_full_8m_w256"]
MASKS = ["ld_block_0p90", "within_chrom_longrange_0p90"]


def tail_mean(curve: Path) -> float:
    rows = [json.loads(line) for line in curve.read_text(encoding="utf-8").splitlines() if line]
    if len(rows) != 16 or int(rows[-1]["step"]) != 4000:
        raise ValueError(f"incomplete pilot curve: {curve}")
    return statistics.fmean(float(row["validation_masked_nll"]) for row in rows[-5:])


def choose_shared_lr(summaries: dict[str, dict[str, float | bool]]) -> tuple[float | None, str]:
    eligible_lrs = [lr for lr in LRS if bool(summaries[f"{lr:.4g}"]["eligible"])]
    if not eligible_lrs:
        return None, "LR_PILOT_NO_ELIGIBLE_LR_REPLAN"
    best = max(float(summaries[f"{lr:.4g}"]["mean_gain_vs_1e_4"]) for lr in eligible_lrs)
    near = [
        lr for lr in eligible_lrs
        if float(summaries[f"{lr:.4g}"]["mean_gain_vs_1e_4"]) >= best - 0.0005
    ]
    selected = min(near)
    status = "LR_PILOT_SHARED_PEAK_SELECTED" if selected > 0.0001 else "LR_PILOT_NO_ACCELERATION_REPLAN"
    return selected, status


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=Path, required=True)
    args = parser.parse_args()
    output = args.run_root / "LR_PILOT_DECISION.v7.2.0.json"
    if output.exists():
        raise SystemExit(f"refusing to overwrite: {output}")

    values: dict[tuple[float, str, str], float] = {}
    failures = sorted(args.run_root.glob("*/FAILURE.json"))
    for result_path in sorted(args.run_root.glob("*/RESULT.json")):
        result = json.loads(result_path.read_text(encoding="utf-8"))
        if result.get("status") != "PASS" or result.get("decision_holdout_read") is not False:
            raise ValueError(f"invalid result contract: {result_path}")
        key = (float(result["peak_lr"]), str(result["model"]), str(result["mask"]))
        if key in values:
            raise ValueError(f"duplicate LR/model/mask cell: {key}")
        values[key] = tail_mean(result_path.parent / "PILOT_CURVE.jsonl")

    expected = {(lr, model, mask) for lr in LRS for model in MODELS for mask in MASKS}
    missing = sorted(expected - set(values))
    if failures or missing:
        decision = {
            "schema_version": "1.0", "protocol_version": "v7.2.0",
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "status": "LR_PILOT_INCOMPLETE_REPLAN", "selected_shared_peak_lr": None,
            "result_count": len(values), "failure_count": len(failures),
            "missing_cells": [list(item) for item in missing], "hgdp_used": False,
            "architecture_decision_permitted": False,
        }
        output.write_text(json.dumps(decision, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(json.dumps(decision, indent=2, sort_keys=True))
        raise SystemExit(4)

    baseline = {(model, mask): values[(0.0001, model, mask)] for model in MODELS for mask in MASKS}
    summaries = {}
    for lr in LRS:
        gains = [baseline[(model, mask)] - values[(lr, model, mask)] for model in MODELS for mask in MASKS]
        eligible = min(gains) >= -0.002
        summaries[f"{lr:.4g}"] = {
            "eligible": eligible, "mean_gain_vs_1e_4": statistics.fmean(gains),
            "median_gain_vs_1e_4": statistics.median(gains),
            "worst_cell_gain_vs_1e_4": min(gains), "best_cell_gain_vs_1e_4": max(gains),
            "cell_gains": {
                f"{model}|{mask}": baseline[(model, mask)] - values[(lr, model, mask)]
                for model in MODELS for mask in MASKS
            },
        }
    selected, status = choose_shared_lr(summaries)
    decision = {
        "schema_version": "1.0", "protocol_version": "v7.2.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(), "status": status,
        "selected_shared_peak_lr": selected, "result_count": len(values), "failure_count": 0,
        "selection_rule": "LOWEST_ELIGIBLE_WITHIN_0.0005_OF_BEST_SIX_CELL_MEAN_GAIN",
        "cell_degradation_limit_nats": 0.002, "summaries": summaries,
        "validation_role": "HISTORICALLY_VIEWED_TUNING_ONLY",
        "decision_holdout_read": False, "hgdp_used": False,
        "architecture_decision_permitted": False,
        "next_authorized_stage": "ALL_12_RUN_C0_SCHEDULE_CONFIRMATION_FROM_STEP_ZERO" if status == "LR_PILOT_SHARED_PEAK_SELECTED" else "REPLAN",
    }
    output.write_text(json.dumps(decision, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(decision, indent=2, sort_keys=True))
    if status != "LR_PILOT_SHARED_PEAK_SELECTED":
        raise SystemExit(4)


if __name__ == "__main__":
    main()
