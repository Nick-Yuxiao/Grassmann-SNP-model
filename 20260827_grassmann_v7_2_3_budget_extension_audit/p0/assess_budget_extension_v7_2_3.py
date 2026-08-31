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
            right = blocks.pop()
            left = blocks.pop()
            weight = int(left[3]) + int(right[3])
            mean = (
                float(left[2]) * int(left[3]) + float(right[2]) * int(right[3])
            ) / weight
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
    median_ratio = statistics.median(ratios)
    return {
        "common_nll_low": low, "common_nll_high": high, "targets": target_rows,
        "median_acceleration_ratio": median_ratio,
        "implied_selected_step_reduction_fraction": 1.0 - 1.0 / median_ratio,
        "min_acceleration_ratio": min(ratios), "max_acceleration_ratio": max(ratios),
        "target_support": "OBSERVED_COMMON_REACHABLE_NLL_RANGE_INTERIOR",
        "censoring": "NONE_BY_CONSTRUCTION_COMMON_SUPPORT_ONLY",
    }


def terminal_summary(
    values: list[tuple[int, float]], train_values: list[tuple[int, float]]
) -> dict[str, object]:
    changes = []
    for start, stop in ((24000, 26000), (26000, 28000), (28000, 30000)):
        change = tail_mean(values, start) - tail_mean(values, stop)
        changes.append({
            "start": start, "stop": stop, "nll_drop": change,
            "absolute_change": abs(change),
        })
    primary_stable = all(float(row["absolute_change"]) <= 0.002 for row in changes)
    all_positive = all(float(row["nll_drop"]) > 0.0 for row in changes)
    first = float(changes[0]["absolute_change"])
    last = float(changes[-1]["absolute_change"])
    last_to_first_ratio = last / first if first > 0.0 else None
    acceleration_flag = bool(
        primary_stable and all_positive and last_to_first_ratio is not None
        and last_to_first_ratio > 1.5
    )
    if not primary_stable:
        shape_class = "NOT_STABLE"
    elif acceleration_flag:
        shape_class = "STABLE_BUT_ACCELERATING"
    else:
        shape_class = "STABLE"
    rolling = [
        statistics.fmean([value for _, value in values][index - 4:index + 1])
        for index in range(4, len(values))
    ]
    terminal_tail = tail_mean(values, 30000)
    degradation = terminal_tail - min(rolling)
    return {
        "changes": changes, "all_absolute_changes_le_0p002": primary_stable,
        "all_nll_drops_positive": all_positive,
        "last_to_first_absolute_change_ratio": last_to_first_ratio,
        "shape_class": shape_class, "shape_flag_is_non_primary": True,
        "terminal_tail5": terminal_tail, "best_rolling_tail5": min(rolling),
        "best_to_terminal_degradation": degradation,
        "instability": degradation > 0.002,
        "train_nll_previous_8_mean": statistics.fmean(
            value for _, value in train_values[-16:-8]
        ),
        "train_nll_final_8_mean": statistics.fmean(
            value for _, value in train_values[-8:]
        ),
    }


def decision_branch(
    sequence_failure: bool,
    instability: bool,
    terminal_adequate: bool,
    selected_shape_flags: list[str],
) -> tuple[str, int | None, str, int]:
    if sequence_failure or instability:
        return "BUDGET_EXTENSION_REPLAN_INSTABILITY", None, "REPLAN", 4
    if not terminal_adequate:
        return (
            "BUDGET_EXTENSION_30K_NOT_ADEQUATE_REPLAN", None,
            "REPLAN_NO_AUTOMATIC_EXTENSION", 6,
        )
    if selected_shape_flags:
        return (
            "BUDGET_EXTENSION_30K_PRIMARY_ADEQUATE_SHAPE_REVIEW", None,
            "BUDGET_REPLAN_SHAPE_REVIEW", 7,
        )
    return (
        "BUDGET_EXTENSION_30K_ADEQUATE", 40000,
        "FREEZE_WARMUP500_STABLE_TO30K_COSINE_TO40K", 0,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--source-root", type=Path, required=True)
    args = parser.parse_args()
    output = args.run_root / "BUDGET_EXTENSION_DECISION.v7.2.3.json"
    if output.exists():
        raise SystemExit(f"refusing to overwrite: {output}")

    source_decision = json.loads(
        (args.source_root / "BUDGET_BRIDGE_DECISION.v7.2.1.json").read_text(encoding="utf-8")
    )
    source_authorized = (
        source_decision.get("status") == "BUDGET_BRIDGE_EXTEND_ALL_TO_30K"
        and source_decision.get("decision_holdout_read") is False
        and source_decision.get("architecture_decision_permitted") is False
    )

    source: dict[tuple[float, str, str], tuple[Path, dict[str, object]]] = {}
    for result_path in sorted(args.source_root.glob("*/RESULT.json")):
        result = json.loads(result_path.read_text(encoding="utf-8"))
        key = (float(result["learning_rate"]), str(result["model"]), str(result["mask"]))
        source[key] = (result_path.parent / "BUDGET_CURVE.jsonl", result)

    extension: dict[tuple[float, str, str], tuple[Path, dict[str, object]]] = {}
    failures = sorted(args.run_root.glob("*/FAILURE.json"))
    lineage_failures = []
    for result_path in sorted(args.run_root.glob("*/RESULT.json")):
        result = json.loads(result_path.read_text(encoding="utf-8"))
        key = (float(result["learning_rate"]), str(result["model"]), str(result["mask"]))
        source_item = source.get(key)
        valid = (
            result.get("status") == "PASS"
            and result.get("protocol_version") == "v7.2.3"
            and result.get("source_step") == 20000
            and result.get("target_step") == 30000
            and result.get("final_step") == 30000
            and result.get("optimizer_resumed") is True
            and result.get("decision_holdout_read") is False
            and result.get("hgdp_used") is False
            and source_item is not None
            and result.get("source_curve_sha256") == source_item[1].get("curve_sha256")
            and result.get("source_checkpoint_sha256") == source_item[1].get("checkpoint_sha256")
        )
        if not valid:
            lineage_failures.append(str(result_path))
        if key in extension:
            lineage_failures.append(f"duplicate extension cell: {key}")
        extension[key] = (result_path.parent / "BUDGET_EXTENSION_CURVE.jsonl", result)

    expected = {(lr, model, mask) for lr in LRS for model in MODELS for mask in MASKS}
    missing_source = sorted(expected - set(source))
    missing_extension = sorted(expected - set(extension))
    if failures or missing_source or missing_extension or lineage_failures or not source_authorized:
        decision = {
            "schema_version": "1.0", "protocol_version": "v7.2.3",
            "status": "BUDGET_EXTENSION_REPLAN_INSTABILITY",
            "result_count": len(extension), "failure_count": len(failures),
            "missing_source": [list(item) for item in missing_source],
            "missing_extension": [list(item) for item in missing_extension],
            "lineage_failures": lineage_failures, "source_authorized": source_authorized,
            "decision_holdout_read": False, "hgdp_used": False,
            "architecture_decision_permitted": False,
        }
        output.write_text(json.dumps(decision, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(json.dumps(decision, indent=2, sort_keys=True))
        raise SystemExit(4)

    combined: dict[tuple[float, str, str], list[tuple[int, float]]] = {}
    train_extension: dict[tuple[float, str, str], list[tuple[int, float]]] = {}
    sequence_failure = False
    for key in sorted(expected):
        source_values = rows(source[key][0], "validation_masked_nll")
        extension_values = rows(extension[key][0], "validation_masked_nll")
        expected_source_steps = list(range(250, 20001, 250))
        expected_extension_steps = list(range(20250, 30001, 250))
        sequence_failure |= [step for step, _ in source_values] != expected_source_steps
        sequence_failure |= [step for step, _ in extension_values] != expected_extension_steps
        combined[key] = source_values + extension_values
        train_extension[key] = rows(extension[key][0], "train_masked_nll_interval")

    acceleration_by_cell = {}
    terminal_by_lr = {"0.0001": {}, "0.0004": {}}
    terminal_adequate = True
    instability = False
    selected_shape_flags = []
    for model in MODELS:
        for mask in MASKS:
            cell = f"{model}|{mask}"
            acceleration_by_cell[cell] = acceleration(
                combined[(0.0001, model, mask)], combined[(0.0004, model, mask)]
            )
            for lr in LRS:
                summary = terminal_summary(
                    combined[(lr, model, mask)], train_extension[(lr, model, mask)]
                )
                terminal_by_lr[f"{lr:.4g}"][cell] = summary
                if lr == 0.0004:
                    terminal_adequate &= bool(summary["all_absolute_changes_le_0p002"])
                    instability |= bool(summary["instability"])
                    if summary["shape_class"] == "STABLE_BUT_ACCELERATING":
                        selected_shape_flags.append(cell)

    status, proposed_total, next_stage, exit_code = decision_branch(
        sequence_failure, instability, terminal_adequate, selected_shape_flags
    )
    ratios = [
        float(row["median_acceleration_ratio"]) for row in acceleration_by_cell.values()
    ]
    decision = {
        "schema_version": "1.0", "protocol_version": "v7.2.3",
        "created_at_utc": datetime.now(timezone.utc).isoformat(), "status": status,
        "source_status": source_decision.get("status"), "source_authorized": source_authorized,
        "result_count": len(extension), "failure_count": 0,
        "sequence_failure": sequence_failure,
        "acceleration_by_cell": acceleration_by_cell,
        "median_of_cell_median_acceleration_ratios": statistics.median(ratios),
        "implied_step_reduction_from_median_ratio": 1.0 - 1.0 / statistics.median(ratios),
        "acceleration_estimand": {
            "targets_per_cell": 5,
            "target_support": "OBSERVED_COMMON_REACHABLE_NLL_RANGE_INTERIOR",
            "censoring": "NONE_BY_CONSTRUCTION_COMMON_SUPPORT_ONLY",
            "global_training_time_multiplier_claim_permitted": False,
        },
        "terminal_1e_4": terminal_by_lr["0.0001"],
        "terminal_4e_4": terminal_by_lr["0.0004"],
        "primary_terminal_cells_pass": sum(
            bool(value["all_absolute_changes_le_0p002"])
            for value in terminal_by_lr["0.0004"].values()
        ),
        "selected_shape_flags": selected_shape_flags,
        "control_terminal_summary_is_descriptive": True,
        "bridge_endpoint_equals_formal_wsd_endpoint": False,
        "proposed_formal_total_steps": proposed_total,
        "next_authorized_stage": next_stage, "decision_holdout_read": False,
        "hgdp_used": False, "architecture_decision_permitted": False,
    }
    output.write_text(json.dumps(decision, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(decision, indent=2, sort_keys=True))
    raise SystemExit(exit_code)


if __name__ == "__main__":
    main()
