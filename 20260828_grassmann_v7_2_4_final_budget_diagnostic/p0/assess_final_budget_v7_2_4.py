from __future__ import annotations

import argparse
import json
import statistics
from datetime import datetime, timezone
from pathlib import Path


MODELS = ["local_attn_8m_w256", "local_attn_gpc_8m_w256", "grassmann_full_8m_w256"]
MASKS = ["ld_block_0p90", "within_chrom_longrange_0p90"]
LEARNING_RATE = 0.0004


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


def terminal_summary(
    values: list[tuple[int, float]], train_values: list[tuple[int, float]]
) -> dict[str, object]:
    changes = []
    for start, stop in ((34000, 36000), (36000, 38000), (38000, 40000)):
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
    terminal_tail = tail_mean(values, 40000)
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
    shape_flags: list[str],
) -> tuple[str, int | None, str, int]:
    if sequence_failure or instability:
        return "FINAL_BUDGET_REPLAN_INSTABILITY", None, "REPLAN", 4
    if not terminal_adequate:
        return (
            "FINAL_BUDGET_40K_NOT_ADEQUATE_STOP", None,
            "NO_AUTOMATIC_EXTENSION_REPLAN_OR_DESCRIPTIVE_A1R", 6,
        )
    if shape_flags:
        return (
            "FINAL_BUDGET_40K_PRIMARY_ADEQUATE_SHAPE_REVIEW", None,
            "FINAL_BUDGET_REPLAN_SHAPE_REVIEW", 7,
        )
    return (
        "FINAL_BUDGET_40K_ADEQUATE", 50000,
        "AUTHORIZE_FORMAL_SCHEDULE_CONTRACT_WARMUP500_STABLE_TO40K_COSINE_TO50K", 0,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--source-root", type=Path, required=True)
    args = parser.parse_args()
    output = args.run_root / "FINAL_BUDGET_DECISION.v7.2.4.json"
    if output.exists():
        raise SystemExit(f"refusing to overwrite: {output}")

    source_decision = json.loads(
        (args.source_root / "BUDGET_EXTENSION_DECISION.v7.2.3.json").read_text(
            encoding="utf-8"
        )
    )
    source_authorized = (
        source_decision.get("status") == "BUDGET_EXTENSION_30K_NOT_ADEQUATE_REPLAN"
        and source_decision.get("result_count") == 12
        and source_decision.get("failure_count") == 0
        and source_decision.get("primary_terminal_cells_pass") == 4
        and source_decision.get("next_authorized_stage") == "REPLAN_NO_AUTOMATIC_EXTENSION"
        and source_decision.get("decision_holdout_read") is False
        and source_decision.get("architecture_decision_permitted") is False
    )

    source: dict[tuple[str, str], tuple[Path, dict[str, object]]] = {}
    source_duplicates = []
    for result_path in sorted(args.source_root.glob("*/RESULT.json")):
        result = json.loads(result_path.read_text(encoding="utf-8"))
        if float(result.get("learning_rate", -1.0)) != LEARNING_RATE:
            continue
        key = (str(result["model"]), str(result["mask"]))
        if key in source:
            source_duplicates.append(key)
        source[key] = (result_path.parent / "BUDGET_EXTENSION_CURVE.jsonl", result)

    final: dict[tuple[str, str], tuple[Path, dict[str, object]]] = {}
    failures = sorted(args.run_root.glob("*/FAILURE.json"))
    lineage_failures = []
    for result_path in sorted(args.run_root.glob("*/RESULT.json")):
        result = json.loads(result_path.read_text(encoding="utf-8"))
        key = (str(result.get("model")), str(result.get("mask")))
        source_item = source.get(key)
        valid = (
            result.get("status") == "PASS"
            and result.get("protocol_version") == "v7.2.4"
            and float(result.get("learning_rate", -1.0)) == LEARNING_RATE
            and result.get("source_step") == 30000
            and result.get("target_step") == 40000
            and result.get("final_step") == 40000
            and result.get("optimizer_resumed") is True
            and result.get("rng_resumed") is True
            and result.get("decision_holdout_read") is False
            and result.get("hgdp_used") is False
            and source_item is not None
            and result.get("source_curve_sha256") == source_item[1].get("curve_sha256")
            and result.get("source_checkpoint_sha256")
            == source_item[1].get("checkpoint_sha256")
        )
        if not valid:
            lineage_failures.append(str(result_path))
        if key in final:
            lineage_failures.append(f"duplicate final cell: {key}")
        final[key] = (result_path.parent / "FINAL_BUDGET_CURVE.jsonl", result)

    expected = {(model, mask) for model in MODELS for mask in MASKS}
    missing_source = sorted(expected - set(source))
    missing_final = sorted(expected - set(final))
    if (
        failures or missing_source or missing_final or lineage_failures
        or source_duplicates or not source_authorized
    ):
        decision = {
            "schema_version": "1.0", "protocol_version": "v7.2.4",
            "status": "FINAL_BUDGET_REPLAN_INSTABILITY",
            "result_count": len(final), "failure_count": len(failures),
            "missing_source": [list(item) for item in missing_source],
            "missing_final": [list(item) for item in missing_final],
            "source_duplicates": [list(item) for item in source_duplicates],
            "lineage_failures": lineage_failures, "source_authorized": source_authorized,
            "decision_holdout_read": False, "hgdp_used": False,
            "architecture_decision_permitted": False,
            "automatic_extension_beyond_40k_permitted": False,
        }
        output.write_text(json.dumps(decision, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(json.dumps(decision, indent=2, sort_keys=True))
        raise SystemExit(4)

    summaries = {}
    sequence_failure = False
    terminal_adequate = True
    instability = False
    shape_flags = []
    for key in sorted(expected):
        source_values = rows(source[key][0], "validation_masked_nll")
        final_values = rows(final[key][0], "validation_masked_nll")
        train_values = rows(final[key][0], "train_masked_nll_interval")
        sequence_failure |= [step for step, _ in source_values] != list(
            range(20250, 30001, 250)
        )
        sequence_failure |= [step for step, _ in final_values] != list(
            range(30250, 40001, 250)
        )
        summary = terminal_summary(final_values, train_values)
        cell = f"{key[0]}|{key[1]}"
        summaries[cell] = summary
        terminal_adequate &= bool(summary["all_absolute_changes_le_0p002"])
        instability |= bool(summary["instability"])
        if summary["shape_class"] == "STABLE_BUT_ACCELERATING":
            shape_flags.append(cell)

    status, proposed_total, next_stage, exit_code = decision_branch(
        sequence_failure, instability, terminal_adequate, shape_flags
    )
    decision = {
        "schema_version": "1.0", "protocol_version": "v7.2.4",
        "created_at_utc": datetime.now(timezone.utc).isoformat(), "status": status,
        "source_status": source_decision.get("status"),
        "source_authorized": source_authorized,
        "result_count": len(final), "failure_count": 0,
        "sequence_failure": sequence_failure,
        "learning_rate": LEARNING_RATE,
        "terminal_4e_4": summaries,
        "primary_terminal_cells_pass": sum(
            bool(value["all_absolute_changes_le_0p002"])
            for value in summaries.values()
        ),
        "selected_shape_flags": shape_flags,
        "proposed_formal_total_steps": proposed_total,
        "proposed_formal_schedule": (
            {
                "warmup_steps": 500, "stable_learning_rate": 0.0004,
                "stable_through_step": 40000, "cosine_end_step": 50000,
                "cosine_end_learning_rate": 0.00004,
            }
            if status == "FINAL_BUDGET_40K_ADEQUATE" else None
        ),
        "next_authorized_stage": next_stage,
        "automatic_extension_beyond_40k_permitted": False,
        "decision_holdout_read": False, "hgdp_used": False,
        "architecture_decision_permitted": False,
    }
    output.write_text(json.dumps(decision, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(decision, indent=2, sort_keys=True))
    raise SystemExit(exit_code)


if __name__ == "__main__":
    main()

