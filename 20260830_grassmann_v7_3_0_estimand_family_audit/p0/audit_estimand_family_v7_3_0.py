from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import statistics
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


MODELS = (
    "local_attn_8m_w256",
    "local_attn_gpc_8m_w256",
    "grassmann_full_8m_w256",
)
MASKS = ("ld_block_0p90", "within_chrom_longrange_0p90")
SELECTED_LR = 0.0004
EXPECTED_CELLS = {(model, mask) for model in MODELS for mask in MASKS}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def curve_rows(path: Path, field: str) -> list[tuple[int, float]]:
    rows: list[tuple[int, float]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line:
            continue
        row = json.loads(line)
        rows.append((int(row["step"]), float(row[field])))
    if not rows:
        raise ValueError(f"empty curve: {path}")
    return rows


def tail_mean(rows: list[tuple[int, float]], through_step: int, count: int = 5) -> float:
    eligible = [value for step, value in rows if step <= through_step]
    if len(eligible) < count:
        raise ValueError(f"fewer than {count} values through step {through_step}")
    return statistics.fmean(eligible[-count:])


def collect_selected_runs(
    root: Path,
    protocol_version: str,
    curve_name: str,
    final_step: int,
    expected_steps: list[int],
) -> dict[tuple[str, str], dict[str, Any]]:
    collected: dict[tuple[str, str], dict[str, Any]] = {}
    for result_path in sorted(root.glob("*/RESULT.json")):
        result = load_json(result_path)
        if float(result.get("learning_rate", -1.0)) != SELECTED_LR:
            continue
        key = (str(result.get("model")), str(result.get("mask")))
        if key in collected:
            raise ValueError(f"duplicate selected-LR cell in {root}: {key}")
        curve_path = result_path.parent / curve_name
        if not curve_path.is_file():
            raise ValueError(f"missing curve: {curve_path}")
        rows = curve_rows(curve_path, "validation_masked_nll")
        actual_steps = [step for step, _ in rows]
        expected = {
            "status": "PASS",
            "protocol_version": protocol_version,
            "model": key[0],
            "mask": key[1],
            "mask_seed": 92001,
            "init_seed": 82001,
            "final_step": final_step,
            "decision_holdout_read": False,
            "hgdp_used": False,
        }
        mismatches = {
            field: (result.get(field), value)
            for field, value in expected.items()
            if result.get(field) != value
        }
        if mismatches:
            raise ValueError(f"result contract mismatch for {key}: {mismatches}")
        if actual_steps != expected_steps:
            raise ValueError(
                f"curve sequence mismatch for {key}: "
                f"observed {actual_steps[:2]}...{actual_steps[-2:]}"
            )
        if sha256(curve_path) != result.get("curve_sha256"):
            raise ValueError(f"curve hash mismatch for {key}: {curve_path}")
        collected[key] = {
            "result": result,
            "result_path": result_path,
            "curve_path": curve_path,
            "rows": rows,
        }
    if set(collected) != EXPECTED_CELLS:
        raise ValueError(
            f"selected-LR factorial mismatch in {root}: "
            f"missing={sorted(EXPECTED_CELLS - set(collected))}, "
            f"extra={sorted(set(collected) - EXPECTED_CELLS)}"
        )
    return collected


def validate_decisions(
    bridge_root: Path, extension_root: Path, final_root: Path
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    bridge = load_json(bridge_root / "BUDGET_BRIDGE_DECISION.v7.2.1.json")
    extension = load_json(extension_root / "BUDGET_EXTENSION_DECISION.v7.2.3.json")
    final = load_json(final_root / "FINAL_BUDGET_DECISION.v7.2.4.json")
    expected = (
        (
            bridge,
            "BUDGET_BRIDGE_EXTEND_ALL_TO_30K",
            12,
            "v7.2.1 bridge",
        ),
        (
            extension,
            "BUDGET_EXTENSION_30K_NOT_ADEQUATE_REPLAN",
            12,
            "v7.2.3 extension",
        ),
        (
            final,
            "FINAL_BUDGET_40K_NOT_ADEQUATE_STOP",
            6,
            "v7.2.4 final",
        ),
    )
    for decision, status, count, label in expected:
        if (
            decision.get("status") != status
            or decision.get("result_count") != count
            or decision.get("failure_count") != 0
            or decision.get("architecture_decision_permitted") is not False
        ):
            raise ValueError(f"{label} decision contract mismatch")
    return bridge, extension, final


def validate_lineage(
    bridge: dict[tuple[str, str], dict[str, Any]],
    extension: dict[tuple[str, str], dict[str, Any]],
    final: dict[tuple[str, str], dict[str, Any]],
) -> None:
    for key in sorted(EXPECTED_CELLS):
        bridge_result = bridge[key]["result"]
        extension_result = extension[key]["result"]
        final_result = final[key]["result"]
        if extension_result.get("source_curve_sha256") != bridge_result.get(
            "curve_sha256"
        ):
            raise ValueError(f"20k-to-30k curve lineage mismatch for {key}")
        if extension_result.get("source_checkpoint_sha256") != bridge_result.get(
            "checkpoint_sha256"
        ):
            raise ValueError(f"20k-to-30k checkpoint lineage mismatch for {key}")
        if final_result.get("source_curve_sha256") != extension_result.get(
            "curve_sha256"
        ):
            raise ValueError(f"30k-to-40k curve lineage mismatch for {key}")
        if final_result.get("source_checkpoint_sha256") != extension_result.get(
            "checkpoint_sha256"
        ):
            raise ValueError(f"30k-to-40k checkpoint lineage mismatch for {key}")


def is_primary_pass(summary: dict[str, Any]) -> bool:
    return bool(summary.get("all_absolute_changes_le_0p002")) and not bool(
        summary.get("instability")
    )


def family_eligibility(
    extension_decision: dict[str, Any], final_decision: dict[str, Any]
) -> tuple[dict[str, dict[str, Any]], list[str]]:
    at_30k = extension_decision["terminal_4e_4"]
    at_40k = final_decision["terminal_4e_4"]
    false_plateaus = sorted(
        cell
        for cell in at_30k
        if is_primary_pass(at_30k[cell]) and not is_primary_pass(at_40k[cell])
    )
    ld_cells = sorted(cell for cell in at_40k if cell.endswith("|ld_block_0p90"))
    long_cells = sorted(
        cell
        for cell in at_40k
        if cell.endswith("|within_chrom_longrange_0p90")
    )
    ld_persistence = all(
        is_primary_pass(at_30k[cell])
        and is_primary_pass(at_40k[cell])
        and at_40k[cell].get("shape_class") == "STABLE"
        for cell in ld_cells
    )
    longrange_still_learning = all(
        not is_primary_pass(at_40k[cell])
        and all(float(change["nll_drop"]) > 0.0 for change in at_40k[cell]["changes"])
        for cell in long_cells
    )
    if len(ld_cells) != 3 or len(long_cells) != 3:
        raise ValueError("unexpected final decision cell count by mask family")
    if not ld_persistence or not longrange_still_learning:
        raise ValueError("observed family pattern does not satisfy v7.3.0 frozen fork")
    families = {
        "ld_block_0p90": {
            "capacity_status": "ELIGIBLE_FOR_CONFIRMATORY_30K",
            "confirmatory_horizon": 30000,
            "basis": "OUT_OF_SAMPLE_PERSISTENCE_30K_TO40K",
            "current_single_seed_is_capacity_evidence": False,
        },
        "within_chrom_longrange_0p90": {
            "capacity_status": "INCONCLUSIVE_BUDGET_NOT_ESTIMABLE",
            "efficiency_status": "DESCRIPTIVE_ONLY",
            "registered_efficiency_budget_points": [20000, 30000, 40000],
            "positive_capacity_upgrade_permitted": False,
            "negative_capacity_upgrade_permitted": False,
        },
    }
    return families, false_plateaus


def boundary_audit(
    extension: dict[tuple[str, str], dict[str, Any]],
    final: dict[tuple[str, str], dict[str, Any]],
) -> dict[str, Any]:
    values: dict[str, float] = {}
    for key in sorted(EXPECTED_CELLS):
        source_rows = extension[key]["rows"]
        final_rows = final[key]["rows"]
        if source_rows[-1][0] != 30000 or final_rows[0][0] != 30250:
            raise ValueError(f"resume boundary sequence mismatch for {key}")
        values[f"{key[0]}|{key[1]}"] = source_rows[-1][1] - final_rows[0][1]
    absolute = [abs(value) for value in values.values()]
    return {
        "role": "BOUNDARY_ASSOCIATED_CHANGE_UPPER_BOUND",
        "pure_resume_discontinuity_estimate": False,
        "confounded_with_optimizer_steps": 250,
        "per_cell_nll_drop": values,
        "maximum_absolute_change": max(absolute),
        "minimum_signed_change": min(values.values()),
        "maximum_signed_change": max(values.values()),
        "future_pre_update_same_checkpoint_revalidation_required": True,
    }


def post_hoc_tail_audit(final_decision: dict[str, Any]) -> dict[str, Any]:
    means = {
        cell: statistics.fmean(
            float(change["nll_drop"]) for change in summary["changes"]
        )
        for cell, summary in sorted(final_decision["terminal_4e_4"].items())
    }
    return {
        "registered_in_v7_2_4": False,
        "role": "POST_HOC_STOP_AND_PLANNING_ONLY",
        "go_direction_permitted": False,
        "mean_last_three_2k_drops_by_cell": means,
        "h5_0p3_delta_min_gate_role": "NOT_REGISTERED_NOT_EXECUTABLE_FOR_GO",
    }


def historical_variance_proxy(c0_root: Path) -> dict[str, Any]:
    grouped: dict[tuple[str, str], list[dict[str, Any]]] = {
        key: [] for key in EXPECTED_CELLS
    }
    for result_path in sorted(c0_root.glob("*/RESULT.json")):
        result = load_json(result_path)
        key = (str(result.get("model")), str(result.get("mask")))
        if key not in grouped or result.get("status") != "PASS":
            continue
        curve_path = result_path.parent / "VALIDATION_CURVE_STEP000250_TO020000.jsonl"
        rows = curve_rows(curve_path, "masked_nll")
        if rows[-1][0] != 20000 or result.get("steps") != 20000:
            raise ValueError(f"historical C0 horizon mismatch: {result_path}")
        grouped[key].append(
            {
                "terminal_tail5": tail_mean(rows, 20000),
                "mask_seed": result.get("mask_seed"),
                "init_seed": result.get("init_seed"),
                "data_seed": result.get("data_seed"),
            }
        )
    return variance_proxy_from_grouped(grouped)


def variance_proxy_from_grouped(
    grouped: dict[tuple[str, str], list[dict[str, Any]]]
) -> dict[str, Any]:
    cells: dict[str, dict[str, Any]] = {}
    sd_by_key: dict[tuple[str, str], float | None] = {}
    for key in sorted(EXPECTED_CELLS):
        observations = grouped[key]
        terminal = [float(row["terminal_tail5"]) for row in observations]
        sd = statistics.stdev(terminal) if len(terminal) >= 2 else None
        sd_by_key[key] = sd
        data_seeds = {row["data_seed"] for row in observations if row["data_seed"] is not None}
        cells[f"{key[0]}|{key[1]}"] = {
            "observation_count": len(observations),
            "terminal_tail5_sd": sd,
            "mask_seed_count": len({row["mask_seed"] for row in observations}),
            "init_seed_count": len({row["init_seed"] for row in observations}),
            "independent_data_seed_count": len(data_seeds),
            "data_seed_recorded": bool(data_seeds),
            "role": "PLANNING_PROXY_ONLY",
        }
    pair_bounds: dict[str, float | None] = {}
    for mask in MASKS:
        for model_a, model_b in itertools.combinations(MODELS, 2):
            sd_a = sd_by_key[(model_a, mask)]
            sd_b = sd_by_key[(model_b, mask)]
            label = f"{model_a}_vs_{model_b}|{mask}"
            pair_bounds[label] = None if sd_a is None or sd_b is None else sd_a + sd_b
    return {
        "source_protocol": "v7.1.13_constant_lr_1e_4",
        "role": "PLANNING_PROXY_ONLY",
        "selected_lr_or_schedule_match": False,
        "replaces_n5_data_seed_pilot": False,
        "sqrt2_is_unconditional_bound": False,
        "unknown_covariance_conservative_bound": "SD_A_PLUS_SD_B",
        "cells": cells,
        "paired_delta_sd_conservative_bounds": pair_bounds,
    }


def write_family_tsv(path: Path, families: dict[str, dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "mask_family",
                "capacity_status",
                "efficiency_status",
                "confirmatory_horizon",
                "basis",
                "architecture_decision_permitted",
            ]
        )
        for mask in MASKS:
            row = families[mask]
            writer.writerow(
                [
                    mask,
                    row["capacity_status"],
                    row.get("efficiency_status", "CONFIRMATORY_DESIGN_REQUIRED"),
                    row.get("confirmatory_horizon", ""),
                    row.get("basis", "TRAINING_BUDGET_NOT_ESTIMABLE"),
                    "false",
                ]
            )


def write_variance_tsv(path: Path, proxy: dict[str, Any]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "cell",
                "observation_count",
                "terminal_tail5_sd",
                "mask_seed_count",
                "init_seed_count",
                "independent_data_seed_count",
                "role",
            ]
        )
        for cell, row in sorted(proxy["cells"].items()):
            writer.writerow(
                [
                    cell,
                    row["observation_count"],
                    "" if row["terminal_tail5_sd"] is None else f'{row["terminal_tail5_sd"]:.12g}',
                    row["mask_seed_count"],
                    row["init_seed_count"],
                    row["independent_data_seed_count"],
                    row["role"],
                ]
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bridge-20k-root", type=Path, required=True)
    parser.add_argument("--extension-30k-root", type=Path, required=True)
    parser.add_argument("--final-40k-root", type=Path, required=True)
    parser.add_argument("--c0-20k-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    audit_path = args.output_dir / "ESTIMAND_FAMILY_AUDIT.v7.3.0.json"
    family_tsv = args.output_dir / "FAMILY_ELIGIBILITY.v7.3.0.tsv"
    variance_tsv = args.output_dir / "VARIANCE_PLANNING_PROXY.v7.3.0.tsv"
    for path in (audit_path, family_tsv, variance_tsv):
        if path.exists():
            raise SystemExit(f"refusing to overwrite: {path}")

    bridge_decision, extension_decision, final_decision = validate_decisions(
        args.bridge_20k_root, args.extension_30k_root, args.final_40k_root
    )
    bridge = collect_selected_runs(
        args.bridge_20k_root,
        "v7.2.1",
        "BUDGET_CURVE.jsonl",
        20000,
        list(range(250, 20001, 250)),
    )
    extension = collect_selected_runs(
        args.extension_30k_root,
        "v7.2.3",
        "BUDGET_EXTENSION_CURVE.jsonl",
        30000,
        list(range(20250, 30001, 250)),
    )
    final = collect_selected_runs(
        args.final_40k_root,
        "v7.2.4",
        "FINAL_BUDGET_CURVE.jsonl",
        40000,
        list(range(30250, 40001, 250)),
    )
    validate_lineage(bridge, extension, final)
    families, false_plateaus = family_eligibility(extension_decision, final_decision)
    boundary = boundary_audit(extension, final)
    variance = historical_variance_proxy(args.c0_20k_root)
    post_hoc_tail = post_hoc_tail_audit(final_decision)

    audit = {
        "schema_version": "1.0",
        "protocol_version": "v7.3.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "ESTIMAND_FAMILY_AUDIT_PASS_NO_GPU_AUTHORIZED",
        "input_decision_statuses": {
            "v7.2.1": bridge_decision["status"],
            "v7.2.3": extension_decision["status"],
            "v7.2.4": final_decision["status"],
        },
        "selected_lr_factorial_cells": 6,
        "selected_learning_rate": SELECTED_LR,
        "lineage_20k_30k_40k": "PASS",
        "families": families,
        "false_plateau_cells_30k_pass_40k_fail": false_plateaus,
        "post_hoc_tail_audit": post_hoc_tail,
        "resume_boundary_audit": boundary,
        "historical_variance_proxy": variance,
        "independent_replication_unit": "DATA_SEED",
        "registered_efficiency_budget_points": [20000, 30000, 40000],
        "efficiency_positive_capacity_upgrade_permitted": False,
        "efficiency_negative_capacity_upgrade_permitted": False,
        "constant_lr_gate_reusable_for_decay": False,
        "gpu_used": False,
        "gpu_authorized": False,
        "formal_a1r_authorized": False,
        "n5_data_seed_pilot_authorized": False,
        "hapnest_authorized": False,
        "decision_holdout_read": False,
        "hgdp_used": False,
        "architecture_ranking_emitted": False,
        "architecture_decision_permitted": False,
        "next_authorized_stage": "DRAFT_NEW_EXPERIMENTAL_CONTRACT_ONLY",
        "input_decision_sha256": {
            "v7.2.1": sha256(
                args.bridge_20k_root / "BUDGET_BRIDGE_DECISION.v7.2.1.json"
            ),
            "v7.2.3": sha256(
                args.extension_30k_root / "BUDGET_EXTENSION_DECISION.v7.2.3.json"
            ),
            "v7.2.4": sha256(
                args.final_40k_root / "FINAL_BUDGET_DECISION.v7.2.4.json"
            ),
        },
    }
    write_family_tsv(family_tsv, families)
    write_variance_tsv(variance_tsv, variance)
    audit_path.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
