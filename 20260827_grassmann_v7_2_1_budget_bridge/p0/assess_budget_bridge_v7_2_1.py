from __future__ import annotations

import argparse
import json
import statistics
from datetime import datetime, timezone
from pathlib import Path


MODELS = ["local_attn_8m_w256", "local_attn_gpc_8m_w256", "grassmann_full_8m_w256"]
MASKS = ["ld_block_0p90", "within_chrom_longrange_0p90"]
LRS = [0.0001, 0.0004]


def rows(path: Path, field: str) -> list[tuple[int, float]]:
    result = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line:
            row = json.loads(line)
            result.append((int(row["step"]), float(row[field])))
    return result


def tail_mean(values: list[tuple[int, float]], step: int, count: int = 5) -> float:
    eligible = [value for current, value in values if current <= step]
    if len(eligible) < count:
        raise ValueError(f"fewer than {count} values through step {step}")
    return statistics.fmean(eligible[-count:])


def isotonic_nonincreasing(values: list[float]) -> list[float]:
    blocks: list[list[float | int]] = []
    for index, value in enumerate(values):
        blocks.append([index, index + 1, -float(value), 1])
        while len(blocks) >= 2 and float(blocks[-2][2]) > float(blocks[-1][2]):
            right = blocks.pop(); left = blocks.pop()
            weight = int(left[3]) + int(right[3])
            mean = (float(left[2]) * int(left[3]) + float(right[2]) * int(right[3])) / weight
            blocks.append([int(left[0]), int(right[1]), mean, weight])
    fitted = [0.0] * len(values)
    for start, stop, mean, _ in blocks:
        for index in range(int(start), int(stop)):
            fitted[index] = -float(mean)
    return fitted


def acceleration(control: list[tuple[int, float]], selected: list[tuple[int, float]]) -> dict[str, object]:
    c_steps = [step for step, _ in control]
    s_steps = [step for step, _ in selected]
    if c_steps != s_steps:
        raise ValueError("paired curves have different evaluation steps")
    c_fit = isotonic_nonincreasing([value for _, value in control])
    s_fit = isotonic_nonincreasing([value for _, value in selected])
    low = max(min(c_fit), min(s_fit))
    high = min(max(c_fit), max(s_fit))
    if not low < high:
        raise ValueError("paired curves lack a common NLL range")
    targets = [high - (high - low) * fraction / 6 for fraction in range(1, 6)]
    target_rows = []
    for target in targets:
        c_time = next(step for step, value in zip(c_steps, c_fit) if value <= target)
        s_time = next(step for step, value in zip(s_steps, s_fit) if value <= target)
        target_rows.append({
            "target_nll": target, "control_step": c_time, "selected_step": s_time,
            "acceleration_ratio": c_time / s_time,
        })
    ratios = [float(row["acceleration_ratio"]) for row in target_rows]
    return {
        "common_nll_low": low, "common_nll_high": high, "targets": target_rows,
        "median_acceleration_ratio": statistics.median(ratios),
        "min_acceleration_ratio": min(ratios), "max_acceleration_ratio": max(ratios),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--lr-pilot-run-root", type=Path, required=True)
    args = parser.parse_args()
    output = args.run_root / "BUDGET_BRIDGE_DECISION.v7.2.1.json"
    if output.exists():
        raise SystemExit(f"refusing to overwrite: {output}")

    bridge: dict[tuple[float, str, str], Path] = {}
    failures = sorted(args.run_root.glob("*/FAILURE.json"))
    for result_path in sorted(args.run_root.glob("*/RESULT.json")):
        result = json.loads(result_path.read_text(encoding="utf-8"))
        if result.get("status") != "PASS" or result.get("decision_holdout_read") is not False:
            raise ValueError(f"invalid bridge result: {result_path}")
        key = (float(result["learning_rate"]), str(result["model"]), str(result["mask"]))
        if key in bridge:
            raise ValueError(f"duplicate bridge cell: {key}")
        bridge[key] = result_path.parent / "BUDGET_CURVE.jsonl"
    expected = {(lr, model, mask) for lr in LRS for model in MODELS for mask in MASKS}

    pilot: dict[tuple[float, str, str], Path] = {}
    for result_path in sorted(args.lr_pilot_run_root.glob("*/RESULT.json")):
        result = json.loads(result_path.read_text(encoding="utf-8"))
        key = (float(result["peak_lr"]), str(result["model"]), str(result["mask"]))
        if key[0] in LRS:
            pilot[key] = result_path.parent / "PILOT_CURVE.jsonl"
    missing_bridge = sorted(expected - set(bridge))
    missing_pilot = sorted(expected - set(pilot))
    if failures or missing_bridge or missing_pilot:
        decision = {
            "schema_version": "1.0", "protocol_version": "v7.2.1",
            "status": "BUDGET_BRIDGE_REPLAN_INSTABILITY", "result_count": len(bridge),
            "failure_count": len(failures), "missing_bridge": [list(x) for x in missing_bridge],
            "missing_pilot": [list(x) for x in missing_pilot], "hgdp_used": False,
            "decision_holdout_read": False, "architecture_decision_permitted": False,
        }
        output.write_text(json.dumps(decision, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(json.dumps(decision, indent=2, sort_keys=True))
        raise SystemExit(4)

    reproduction = {}
    acceleration_by_cell = {}
    terminal = {}
    reproduction_failed = False
    instability = False
    terminal_adequate = True
    for model in MODELS:
        for mask in MASKS:
            cell = f"{model}|{mask}"
            reproduction[cell] = {}
            paired_curves = {}
            for lr in LRS:
                key = (lr, model, mask)
                bridge_val = rows(bridge[key], "validation_masked_nll")
                pilot_val = rows(pilot[key], "validation_masked_nll")
                difference = tail_mean(bridge_val, 4000) - tail_mean(pilot_val, 4000)
                passed = abs(difference) <= 0.001
                reproduction_failed |= not passed
                reproduction[cell][f"{lr:.4g}"] = {
                    "bridge_tail5": tail_mean(bridge_val, 4000),
                    "pilot_tail5": tail_mean(pilot_val, 4000),
                    "difference": difference, "pass": passed,
                }
                paired_curves[lr] = bridge_val
            acceleration_by_cell[cell] = acceleration(paired_curves[0.0001], paired_curves[0.0004])

            selected_curve = paired_curves[0.0004]
            changes = []
            for start, stop in ((14000, 16000), (16000, 18000), (18000, 20000)):
                change = tail_mean(selected_curve, start) - tail_mean(selected_curve, stop)
                changes.append({"start": start, "stop": stop, "nll_drop": change, "absolute_change": abs(change)})
            all_stable = all(float(row["absolute_change"]) <= 0.002 for row in changes)
            terminal_adequate &= all_stable
            rolling = [
                statistics.fmean([value for _, value in selected_curve][index - 4:index + 1])
                for index in range(4, len(selected_curve))
            ]
            terminal_tail = tail_mean(selected_curve, 20000)
            degradation = terminal_tail - min(rolling)
            cell_instability = degradation > 0.002
            instability |= cell_instability
            train_curve = rows(bridge[(0.0004, model, mask)], "train_masked_nll_interval")
            terminal[cell] = {
                "changes": changes, "all_absolute_changes_le_0p002": all_stable,
                "terminal_tail5": terminal_tail, "best_rolling_tail5": min(rolling),
                "best_to_terminal_degradation": degradation,
                "instability": cell_instability,
                "train_nll_previous_8_mean": statistics.fmean(value for _, value in train_curve[-16:-8]),
                "train_nll_final_8_mean": statistics.fmean(value for _, value in train_curve[-8:]),
            }

    if reproduction_failed or instability:
        status = "BUDGET_BRIDGE_REPLAN_INSTABILITY"
        proposed_total = None
        next_stage = "REPLAN"
        exit_code = 4
    elif terminal_adequate:
        status = "BUDGET_BRIDGE_20K_ADEQUATE"
        proposed_total = 30000
        next_stage = "FREEZE_WARMUP500_STABLE_TO20K_COSINE_TO30K"
        exit_code = 0
    else:
        status = "BUDGET_BRIDGE_EXTEND_ALL_TO_30K"
        proposed_total = None
        next_stage = "RESUME_ALL_12_CONSTANT_LR_RUNS_TO30K_AFTER_FRESH_GPU_AUDIT"
        exit_code = 5
    ratios = [float(row["median_acceleration_ratio"]) for row in acceleration_by_cell.values()]
    decision = {
        "schema_version": "1.0", "protocol_version": "v7.2.1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(), "status": status,
        "result_count": len(bridge), "failure_count": 0,
        "reproduction": reproduction, "acceleration_by_cell": acceleration_by_cell,
        "median_of_cell_median_acceleration_ratios": statistics.median(ratios),
        "terminal_4e_4": terminal, "proposed_formal_total_steps": proposed_total,
        "next_authorized_stage": next_stage, "decision_holdout_read": False,
        "hgdp_used": False, "architecture_decision_permitted": False,
    }
    output.write_text(json.dumps(decision, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(decision, indent=2, sort_keys=True))
    raise SystemExit(exit_code)


if __name__ == "__main__":
    main()

