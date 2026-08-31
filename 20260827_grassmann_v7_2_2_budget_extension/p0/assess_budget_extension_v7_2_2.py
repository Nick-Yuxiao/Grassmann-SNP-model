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
    return {
        "common_nll_low": low, "common_nll_high": high, "targets": target_rows,
        "median_acceleration_ratio": statistics.median(ratios),
        "min_acceleration_ratio": min(ratios), "max_acceleration_ratio": max(ratios),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--source-root", type=Path, required=True)
    args = parser.parse_args()
    output = args.run_root / "BUDGET_EXTENSION_DECISION.v7.2.2.json"
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
            and result.get("protocol_version") == "v7.2.2"
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
            "schema_version": "1.0", "protocol_version": "v7.2.2",
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
    terminal = {}
    terminal_adequate = True
    instability = False
    for model in MODELS:
        for mask in MASKS:
            cell = f"{model}|{mask}"
            acceleration_by_cell[cell] = acceleration(
                combined[(0.0001, model, mask)], combined[(0.0004, model, mask)]
            )
            selected = combined[(0.0004, model, mask)]
            changes = []
            for start, stop in ((24000, 26000), (26000, 28000), (28000, 30000)):
                change = tail_mean(selected, start) - tail_mean(selected, stop)
                changes.append({
                    "start": start, "stop": stop, "nll_drop": change,
                    "absolute_change": abs(change),
                })
            all_stable = all(float(row["absolute_change"]) <= 0.002 for row in changes)
            terminal_adequate &= all_stable
            rolling = [
                statistics.fmean([value for _, value in selected][index - 4:index + 1])
                for index in range(4, len(selected))
            ]
            terminal_tail = tail_mean(selected, 30000)
            degradation = terminal_tail - min(rolling)
            cell_instability = degradation > 0.002
            instability |= cell_instability
            train_values = train_extension[(0.0004, model, mask)]
            terminal[cell] = {
                "changes": changes, "all_absolute_changes_le_0p002": all_stable,
                "terminal_tail5": terminal_tail, "best_rolling_tail5": min(rolling),
                "best_to_terminal_degradation": degradation,
                "instability": cell_instability,
                "train_nll_previous_8_mean": statistics.fmean(
                    value for _, value in train_values[-16:-8]
                ),
                "train_nll_final_8_mean": statistics.fmean(
                    value for _, value in train_values[-8:]
                ),
            }

    if sequence_failure or instability:
        status = "BUDGET_EXTENSION_REPLAN_INSTABILITY"
        proposed_total = None
        next_stage = "REPLAN"
        exit_code = 4
    elif terminal_adequate:
        status = "BUDGET_EXTENSION_30K_ADEQUATE"
        proposed_total = 40000
        next_stage = "FREEZE_WARMUP500_STABLE_TO30K_COSINE_TO40K"
        exit_code = 0
    else:
        status = "BUDGET_EXTENSION_30K_NOT_ADEQUATE_REPLAN"
        proposed_total = None
        next_stage = "REPLAN_NO_AUTOMATIC_EXTENSION"
        exit_code = 6
    ratios = [
        float(row["median_acceleration_ratio"]) for row in acceleration_by_cell.values()
    ]
    decision = {
        "schema_version": "1.0", "protocol_version": "v7.2.2",
        "created_at_utc": datetime.now(timezone.utc).isoformat(), "status": status,
        "source_status": source_decision.get("status"), "source_authorized": source_authorized,
        "result_count": len(extension), "failure_count": 0,
        "sequence_failure": sequence_failure,
        "acceleration_by_cell": acceleration_by_cell,
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
