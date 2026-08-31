from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def add_parent_sources() -> None:
    metadata = json.loads((ROOT / "PARENT_EVIDENCE.json").read_text(encoding="utf-8"))
    rank_parent = ROOT.parent / metadata["implementation_reference"]
    if sha256_file(rank_parent / "MANIFEST.sha256") != metadata["implementation_reference_manifest_sha256"]:
        raise RuntimeError("implementation-reference manifest mismatch")
    r1_parent = rank_parent.parent / json.loads(
        (rank_parent / "PARENT_EVIDENCE.json").read_text(encoding="utf-8")
    )["code_parent"]
    gc_parent = r1_parent.parent / json.loads(
        (r1_parent / "PARENT_EVIDENCE.json").read_text(encoding="utf-8")
    )["code_parent"]
    sys.path[:0] = [str(ROOT / "src"), str(gc_parent / "src")]


def atomic_json(path: Path, value: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n"
    )
    temporary.replace(path)


def quantiles(values: list[float]) -> dict[str, float]:
    array = np.asarray(values, dtype=float)
    return {
        "q10": float(np.quantile(array, 0.10)),
        "median": float(np.quantile(array, 0.50)),
        "q90": float(np.quantile(array, 0.90)),
    }


def planning_label(fraction: float, config: dict[str, object]) -> str:
    thresholds = config["global_classification"]
    if fraction < thresholds["stop_review_below"]:
        return "STOP_REVIEW"
    if fraction < thresholds["caution_below"]:
        return "CAUTION"
    return "PROVISIONALLY_FEASIBLE"


def markdown_report(summary: dict[str, object], rows: list[dict[str, object]]) -> str:
    selected = [
        row
        for row in rows
        if row["population_relative_gap"] == 0.10 and row["true_angle_deg"] == 20.0
    ]
    selected.sort(key=lambda row: (row["n"], row["maf"]))
    lines = [
        "# Grassmann v6.1 非正式设计可行性结果",
        "",
        "_Classification: `EXPLORATORY_NON_EVIDENCE`; this report is permanently excluded from the formal evidence chain._",
        "",
        "---",
        "",
        "## 📋 结论",
        "",
        f"机器生成的非正式 planning label 是 `{summary['planning_label']}`。"
        f"在预先定义的 planning subset 中，workable grid-cell fraction 为 "
        f"`{summary['workable_planning_cell_fraction']:.3f}`。这只是等权设计格比例，不是真实候选比例。",
        "",
        "> ⚠️ **证据边界：** `null_exceedance_rate` 是 12-replicate Monte Carlo detectability proxy，"
        "不是 p 值、正式 power、FWER 或方法排名。",
        "",
        "## 📊 中心边界切片",
        "",
        "下表固定 population gap `0.10`、true angle `20°`，展示样本量和 MAF 对可估性的影响。完整网格见 `summary.tsv`。",
        "",
        "| n | MAF | Expected n2 | Group pass | Gap pass | Joint estimable | Null exceedance | Median fitted angle |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in selected:
        lines.append(
            f"| {row['n']} | {row['maf']:.2f} | {row['expected_n2']:.1f} | "
            f"{row['group_count_pass_rate']:.3f} | {row['gap_pass_rate']:.3f} | "
            f"{row['joint_estimable_rate']:.3f} | {row['null_exceedance_rate']:.3f} | "
            f"{row['estimated_angle_median_deg']:.2f}° |"
        )
    lines.extend(
        [
            "",
            "## 🔍 解释路径",
            "",
            "```mermaid",
            "flowchart LR",
            "    accTitle: Feasibility result interpretation",
            "    accDescr: The result can trigger a stop or broad redesign, while an evidence firewall prevents it from selecting formal settings or entering calibration and power claims.",
            "",
            "    result([📊 NON_EVIDENCE result]) --> review{🔍 Planning review}",
            "    review -->|Weak| stop([⚠️ Stop or redesign])",
            "    review -->|Plausible| design([📝 Continue T09 design])",
            "    result -.-> firewall[🔒 Evidence firewall]",
            "    firewall -.-> blocked([🚫 No formal ranking])",
            "",
            "    classDef action fill:#dbeafe,stroke:#2563eb,stroke-width:2px,color:#1e3a5f",
            "    classDef warning fill:#fef9c3,stroke:#ca8a04,stroke-width:2px,color:#713f12",
            "    classDef danger fill:#fee2e2,stroke:#dc2626,stroke-width:2px,color:#7f1d1d",
            "    classDef neutral fill:#f3f4f6,stroke:#6b7280,stroke-width:2px,color:#1f2937",
            "",
            "    class result,design action",
            "    class review,stop warning",
            "    class firewall neutral",
            "    class blocked danger",
            "```",
            "",
            "## 🚫 不可作出的结论",
            "",
            "- 不能把 grid-cell fraction 解释为真实候选中可检测者的比例",
            "- 不能据此选择正式样本量、MAF cut、effect size、seed、rank 或 gap threshold",
            "- 不能把结果并入 GC-screen、GC-final、T14、T16 或真实 phenotype 证据",
            "- 不能把继续设计解释为方法已通过 calibration",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    add_parent_sources()

    from design_feasibility import generate_scenario, fit_scenario, wilson_interval

    config = json.loads((ROOT / "config" / "FEASIBILITY_CONFIG.json").read_text(encoding="utf-8"))
    firewall = json.loads((ROOT / "EVIDENCE_FIREWALL.json").read_text(encoding="utf-8"))
    if config["classification"] != "EXPLORATORY_NON_EVIDENCE" or firewall["formal_evidence_eligible"]:
        raise RuntimeError("evidence firewall is not active")
    if os.environ.get("CUDA_VISIBLE_DEVICES") not in (None, "", "-1"):
        raise RuntimeError("feasibility run is CPU-only")

    output = Path(args.output_dir).resolve()
    results_root = (ROOT / "results").resolve()
    if output == results_root or results_root not in output.parents:
        raise ValueError("output must be a new subdirectory under this package's results directory")
    if output.exists():
        raise FileExistsError("duplicate output directory refused")
    output.mkdir(parents=True)

    started = time.perf_counter()
    raw_rows: list[dict[str, object]] = []
    dimensions = (
        len(config["sample_sizes"]),
        len(config["mafs"]),
        len(config["population_relative_gaps"]),
    )
    for n_index, n in enumerate(config["sample_sizes"]):
        for maf_index, maf in enumerate(config["mafs"]):
            for gap_index, population_gap in enumerate(config["population_relative_gaps"]):
                stratum_index = (n_index * dimensions[1] + maf_index) * dimensions[2] + gap_index
                for replicate in range(config["replicates_per_cell"]):
                    seed = config["seed_base"] + stratum_index * config["replicates_per_cell"] + replicate
                    for angle in config["true_max_principal_angles_deg"]:
                        scenario = generate_scenario(
                            seed=seed,
                            n=n,
                            maf=maf,
                            population_relative_gap=population_gap,
                            true_max_principal_angle_deg=angle,
                            effect_scale=config["effect_scale"],
                            residual_sd=config["residual_sd"],
                            conditional_ld_rhos_by_dosage=tuple(config["conditional_ld_rhos_by_dosage"]),
                            residual_scales_by_dosage=tuple(config["residual_scales_by_dosage"]),
                            region_features=config["region_features"],
                            traits=config["traits"],
                            covariate_count=config["covariates"],
                        )
                        fitted = fit_scenario(
                            scenario,
                            rank=config["rank"],
                            ridge_lambda=config["ridge_lambda"],
                            minimum_group_count=config["minimum_group_count"],
                            minimum_fitted_rank_gap=config["minimum_fitted_rank_gap"],
                        )
                        raw_rows.append(
                            {
                                "classification": "EXPLORATORY_NON_EVIDENCE",
                                "run_id": f"feas:{n}:{maf:.2f}:{population_gap:.2f}:{angle}:{replicate}",
                                "seed": seed,
                                "replicate": replicate,
                                "n": n,
                                "maf": maf,
                                "population_relative_gap": population_gap,
                                "true_angle_deg": float(angle),
                                "n0": fitted.genotype_counts[0],
                                "n1": fitted.genotype_counts[1],
                                "n2": fitted.genotype_counts[2],
                                "minimum_group_count": fitted.minimum_group_count,
                                "group_count_eligible": fitted.group_count_eligible,
                                "fitted_minimum_gap": fitted.fitted_minimum_gap,
                                "gap_eligible": fitted.gap_eligible,
                                "jointly_estimable": fitted.jointly_estimable,
                                "true_minimum_gap": fitted.true_minimum_gap,
                                "true_direction_score": fitted.true_direction_score,
                                "true_max_principal_angle_deg": fitted.true_max_principal_angle_deg,
                                "fitted_direction_score": fitted.fitted_direction_score,
                                "gated_direction_score": fitted.gated_direction_score,
                                "fitted_max_principal_angle_deg": fitted.fitted_max_principal_angle_deg,
                            }
                        )

    grouped: dict[tuple[int, float, float], list[dict[str, object]]] = defaultdict(list)
    for row in raw_rows:
        grouped[(row["n"], row["maf"], row["population_relative_gap"])].append(row)
    null_thresholds = {}
    for stratum, stratum_rows in grouped.items():
        null_scores = [row["gated_direction_score"] for row in stratum_rows if row["true_angle_deg"] == 0.0]
        null_thresholds[stratum] = float(
            np.quantile(null_scores, config["null_reference_quantile"], method="higher")
        )

    cell_rows: list[dict[str, object]] = []
    cells: dict[tuple[int, float, float, float], list[dict[str, object]]] = defaultdict(list)
    for row in raw_rows:
        cells[(row["n"], row["maf"], row["population_relative_gap"], row["true_angle_deg"])].append(row)
    for (n, maf, population_gap, angle), values in sorted(cells.items()):
        total = len(values)
        group_successes = sum(row["group_count_eligible"] for row in values)
        gap_successes = sum(row["gap_eligible"] for row in values)
        joint_successes = sum(row["jointly_estimable"] for row in values)
        threshold = null_thresholds[(n, maf, population_gap)]
        exceedances = sum(
            row["jointly_estimable"] and row["gated_direction_score"] > threshold for row in values
        )
        joint_ci = wilson_interval(joint_successes, total)
        exceedance_ci = wilson_interval(exceedances, total)
        angle_q = quantiles([row["fitted_max_principal_angle_deg"] for row in values])
        gap_q = quantiles([row["fitted_minimum_gap"] for row in values])
        cell_rows.append(
            {
                "classification": "EXPLORATORY_NON_EVIDENCE",
                "n": n,
                "maf": maf,
                "population_relative_gap": population_gap,
                "true_angle_deg": angle,
                "replicates": total,
                "expected_n0": n * (1.0 - maf) ** 2,
                "expected_n1": n * 2.0 * maf * (1.0 - maf),
                "expected_n2": n * maf**2,
                "group_count_pass_rate": group_successes / total,
                "gap_pass_rate": gap_successes / total,
                "joint_estimable_rate": joint_successes / total,
                "joint_estimable_ci_low": joint_ci[0],
                "joint_estimable_ci_high": joint_ci[1],
                "null_reference_95": threshold,
                "null_exceedance_rate": exceedances / total,
                "null_exceedance_ci_low": exceedance_ci[0],
                "null_exceedance_ci_high": exceedance_ci[1],
                "estimated_angle_q10_deg": angle_q["q10"],
                "estimated_angle_median_deg": angle_q["median"],
                "estimated_angle_q90_deg": angle_q["q90"],
                "fitted_gap_q10": gap_q["q10"],
                "fitted_gap_median": gap_q["median"],
                "fitted_gap_q90": gap_q["q90"],
            }
        )

    planning = config["planning_subset"]
    workable_rule = config["workable_cell_rule"]
    planning_rows = [
        row
        for row in cell_rows
        if row["population_relative_gap"] >= planning["minimum_population_gap"]
        and row["true_angle_deg"] >= planning["minimum_true_angle_deg"]
    ]
    workable = [
        row
        for row in planning_rows
        if row["joint_estimable_rate"] >= workable_rule["minimum_joint_estimable_rate"]
        and row["null_exceedance_rate"] >= workable_rule["minimum_null_exceedance_rate"]
    ]
    workable_fraction = len(workable) / len(planning_rows)
    summary = {
        "classification": "EXPLORATORY_NON_EVIDENCE",
        "planning_label": planning_label(workable_fraction, config),
        "workable_planning_cells": len(workable),
        "planning_cells": len(planning_rows),
        "workable_planning_cell_fraction": workable_fraction,
        "grid_cells": len(cell_rows),
        "angle_specific_fits": len(raw_rows),
        "unique_paired_base_families": len({row["seed"] for row in raw_rows}),
        "replicates_per_cell": config["replicates_per_cell"],
        "paired_angle_seed_reuse": True,
        "note_on_replication": "Each cell has 12 independent family replicates; angle cells are paired through 576 shared base seeds and are not 2,880 globally independent families.",
        "grid_weight_interpretation": config["grid_weight_interpretation"],
        "formal_evidence_eligible": False,
        "does_not_authorize": config["does_not_authorize"],
        "elapsed_seconds": time.perf_counter() - started,
    }

    with (output / "per_run.jsonl").open("x", encoding="utf-8", newline="\n") as stream:
        for row in raw_rows:
            stream.write(json.dumps(row, sort_keys=True) + "\n")
    with (output / "summary.tsv").open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(cell_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(cell_rows)
    atomic_json(output / "SUMMARY.json", summary)
    (output / "FEASIBILITY_REPORT.md").write_text(
        markdown_report(summary, cell_rows), encoding="utf-8", newline="\n"
    )
    atomic_json(
        output / "ENVIRONMENT.json",
        {
            "classification": "EXPLORATORY_NON_EVIDENCE",
            "python": sys.version,
            "numpy": np.__version__,
            "platform": platform.platform(),
            "accelerator_used": False,
            "config_sha256": sha256_file(ROOT / "config" / "FEASIBILITY_CONFIG.json"),
            "package_manifest_sha256": sha256_file(ROOT / "MANIFEST.sha256"),
        },
    )
    result_files = sorted(path for path in output.iterdir() if path.name != "RESULT_MANIFEST.sha256")
    manifest = "".join(f"{sha256_file(path)}  {path.name}\n" for path in result_files)
    (output / "RESULT_MANIFEST.sha256").write_text(manifest, encoding="utf-8", newline="\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
