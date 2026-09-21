#!/usr/bin/env python3
"""B6: trait preflight before any trait is frozen for qualification.

For each candidate trait this reports effective N, missingness, distribution and
extremes, sentinel-value suspicion, genotyped counts per split, and how much of
the trait the non-local covariates already explain. It then recommends which
traits are usable.

Preflight reads the training-side splits only. The qualification test split and
any bridge holdout stay closed, so choosing traits here cannot leak the answer.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_split import open_text, sha256_file  # noqa: E402


def numeric(value: str) -> float:
    text = (value or "").strip()
    if not text:
        return float("nan")
    try:
        return float(text)
    except ValueError:
        return float("nan")


def read_table(path: Path, id_column: str) -> dict[str, dict[str, str]]:
    table: dict[str, dict[str, str]] = {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit(f"No header in {path}")
        if id_column not in reader.fieldnames:
            raise SystemExit(f"{path.name} has no column {id_column}")
        for row in reader:
            key = (row.get(id_column) or "").strip()
            if key:
                table[key] = row
    return table


def design_matrix(
    keys: list[str],
    covariates: dict[str, dict[str, str]],
    numeric_columns: list[str],
    categorical_columns: list[str],
    squared_columns: list[str],
) -> tuple[np.ndarray, list[str]]:
    blocks = [np.ones((len(keys), 1))]
    names = ["intercept"]
    for column in numeric_columns:
        values = np.array([numeric(covariates[key].get(column, "")) for key in keys])
        finite = values[np.isfinite(values)]
        values = np.where(np.isfinite(values), values, float(finite.mean()) if finite.size else 0.0)
        centre, scale = float(values.mean()), float(values.std())
        scaled = (values - centre) / (scale if scale > 0 else 1.0)
        blocks.append(scaled.reshape(-1, 1))
        names.append(column)
        if column in squared_columns:
            blocks.append((scaled ** 2).reshape(-1, 1))
            names.append(f"{column}^2")
    for column in categorical_columns:
        raw = [str(covariates[key].get(column, "") or "NA").strip() for key in keys]
        levels = sorted(set(raw))
        for level in levels[1:]:
            blocks.append(np.array([1.0 if value == level else 0.0 for value in raw]).reshape(-1, 1))
            names.append(f"{column}={level}")
    return np.hstack(blocks), names


def covariate_fit(design: np.ndarray, y: np.ndarray) -> dict[str, object]:
    gram = design.T @ design
    gram.flat[:: gram.shape[0] + 1] += 1e-8 * float(np.trace(gram)) / gram.shape[0]
    inverse = np.linalg.inv(gram)
    beta = inverse @ (design.T @ y)
    fitted = design @ beta
    residual = y - fitted
    sse = float(residual @ residual)
    sst = float(((y - y.mean()) ** 2).sum())
    degrees = max(len(y) - design.shape[1], 1)
    variance = sse / degrees
    standard_error = np.sqrt(np.maximum(np.diag(inverse) * variance, 0.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        t_statistic = np.where(standard_error > 0, beta / standard_error, 0.0)
    return {
        "r2": 1.0 - sse / sst if sst > 0 else float("nan"),
        "t": np.nan_to_num(t_statistic),
    }


def describe(values: np.ndarray) -> dict[str, float]:
    quantiles = np.percentile(values, [0.1, 1, 25, 50, 75, 99, 99.9])
    centred = values - values.mean()
    scale = values.std()
    skew = float((centred ** 3).mean() / scale ** 3) if scale > 0 else float("nan")
    kurtosis = float((centred ** 4).mean() / scale ** 4 - 3.0) if scale > 0 else float("nan")
    return {
        "mean": float(values.mean()),
        "sd": float(scale),
        "min": float(values.min()),
        "p0.1": float(quantiles[0]),
        "p1": float(quantiles[1]),
        "q1": float(quantiles[2]),
        "median": float(quantiles[3]),
        "q3": float(quantiles[4]),
        "p99": float(quantiles[5]),
        "p99.9": float(quantiles[6]),
        "max": float(values.max()),
        "skew": skew,
        "excess_kurtosis": kurtosis,
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phenotype", type=Path, required=True)
    parser.add_argument("--phenotype-id-column", default="n_eid")
    parser.add_argument("--trait-column", action="append", required=True)
    parser.add_argument("--covariates", type=Path, required=True)
    parser.add_argument("--covariate-id-column", default="eid")
    parser.add_argument("--covariate-column", action="append", default=[])
    parser.add_argument("--categorical-column", action="append", default=[])
    parser.add_argument("--square-column", action="append", default=[])
    parser.add_argument("--split-manifest", type=Path, required=True)
    parser.add_argument("--preflight-split", action="append", default=[],
                        help="Splits preflight may read; defaults to train and validation.")
    parser.add_argument("--min-effective-n", type=int, default=50_000)
    parser.add_argument("--max-missing-rate", type=float, default=0.30)
    parser.add_argument("--min-distinct-values", type=int, default=50)
    parser.add_argument("--max-single-value-share", type=float, default=0.05)
    parser.add_argument("--max-covariate-r2", type=float, default=0.50)
    parser.add_argument("--reject-negative-values", action="store_true",
                        help="Also reject any trait with negative values. Off by default so an "
                             "already standardised trait is not rejected for being negative.")
    parser.add_argument("--trait-cap", type=int, default=3)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    report_path = out_dir / "PREFLIGHT.json"
    if report_path.exists() and not args.overwrite:
        raise SystemExit(f"Refusing to overwrite: {report_path}")

    preflight_splits = set(args.preflight_split or ["train", "validation"])
    with open_text(args.split_manifest) as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    split_counts_all = Counter(row["split"] for row in manifest)
    closed = sorted(set(split_counts_all) - preflight_splits)

    eid_of_sample = {row["sample_id"]: (row.get("eid") or "").strip() for row in manifest}
    preflight_keys = [
        eid_of_sample[row["sample_id"]]
        for row in manifest
        if row["split"] in preflight_splits and eid_of_sample[row["sample_id"]]
    ]

    covariates = read_table(args.covariates, args.covariate_id_column)
    phenotypes = read_table(args.phenotype, args.phenotype_id_column)

    keys = [key for key in preflight_keys if key in covariates]
    genotyped_with_covariates = len(keys)
    per_split_genotyped = Counter(
        row["split"] for row in manifest
        if eid_of_sample[row["sample_id"]] in covariates
    )

    design, design_names = design_matrix(
        keys, covariates, args.covariate_column, args.categorical_column, args.square_column
    )

    traits: list[dict[str, object]] = []
    for trait in args.trait_column:
        raw = np.array([numeric(phenotypes.get(key, {}).get(trait, "")) for key in keys])
        observed = np.isfinite(raw)
        values = raw[observed]
        entry: dict[str, object] = {
            "trait_column": trait,
            "in_phenotype_file": any(trait in row for row in list(phenotypes.values())[:1]),
            "preflight_candidates": len(keys),
            "non_missing": int(observed.sum()),
            "missing_rate": round(1.0 - observed.sum() / max(len(keys), 1), 6),
        }
        if values.size < 100:
            entry["usable"] = False
            entry["reasons"] = ["fewer than 100 observed values in the preflight splits"]
            traits.append(entry)
            continue

        counts = Counter(np.round(values, 6).tolist())
        top_value, top_count = counts.most_common(1)[0]
        mean, sd = float(values.mean()), float(values.std())
        entry["distribution"] = describe(values)
        entry["distinct_values"] = len(counts)
        entry["most_frequent_value"] = {"value": float(top_value), "count": int(top_count),
                                        "share": round(top_count / values.size, 6)}
        entry["negative_values"] = int((values < 0).sum())
        entry["beyond_5sd"] = int((np.abs(values - mean) > 5 * sd).sum()) if sd > 0 else 0
        # Repeated negative integers are worth printing, but they are only evidence of a
        # sentinel code when the trait is otherwise positive: a standardised trait goes
        # below zero legitimately, and a rounded one repeats every integer it covers.
        negative_modes = sorted(
            float(value) for value in counts
            if value < 0
            and float(value) == round(float(value))
            and counts[value] >= max(10, 0.001 * values.size)
        )
        entry["negative_integer_modes"] = negative_modes[:10]
        median = float(np.median(values))
        entry["sentinel_suspects"] = (
            [value for value in negative_modes if counts[value] / values.size >= 0.005]
            if median > 0
            else []
        )

        fit = covariate_fit(design[observed], values)
        t_statistic = fit["t"]
        entry["covariate_model"] = {
            "r2": round(float(fit["r2"]), 6),
            "terms": len(design_names),
            "largest_abs_t": {
                name: round(float(abs(t_statistic[index])), 3)
                for index, name in sorted(
                    enumerate(design_names), key=lambda item: -abs(t_statistic[item[0]])
                )[:8]
            },
        }

        reasons: list[str] = []
        if entry["non_missing"] < args.min_effective_n:
            reasons.append(f"effective N {entry['non_missing']} below {args.min_effective_n}")
        if entry["missing_rate"] > args.max_missing_rate:
            reasons.append(f"missing rate {entry['missing_rate']} above {args.max_missing_rate}")
        if entry["distinct_values"] < args.min_distinct_values:
            reasons.append(f"only {entry['distinct_values']} distinct values; heavily rounded")
        if entry["most_frequent_value"]["share"] > args.max_single_value_share:  # type: ignore[index]
            reasons.append("a single value dominates; check for a sentinel or a floor effect")
        if entry["sentinel_suspects"]:
            reasons.append(
                f"repeated negative integer values {entry['sentinel_suspects']} look like "
                "sentinel codes rather than measurements"
            )
        if float(fit["r2"]) > args.max_covariate_r2:
            reasons.append(f"covariates already explain R2={float(fit['r2']):.3f}")
        if args.reject_negative_values and entry["negative_values"] > 0:
            reasons.append(f"{entry['negative_values']} negative values with the strict rule on")
        entry["usable"] = not reasons
        entry["reasons"] = reasons
        traits.append(entry)

    usable = [item for item in traits if item["usable"]]
    usable.sort(key=lambda item: -int(item["non_missing"]))
    recommended = [str(item["trait_column"]) for item in usable[: args.trait_cap]]

    report = {
        "classification": "TASK_B_TRAIT_PREFLIGHT",
        "identifiers_included": False,
        "splits_read": sorted(preflight_splits),
        "splits_left_closed": closed,
        "split_sizes": dict(sorted(split_counts_all.items())),
        "genotyped_with_covariates_per_split": dict(sorted(per_split_genotyped.items())),
        "preflight_participants": genotyped_with_covariates,
        "covariate_terms": design_names,
        "thresholds": {
            "min_effective_n": args.min_effective_n,
            "max_missing_rate": args.max_missing_rate,
            "min_distinct_values": args.min_distinct_values,
            "max_single_value_share": args.max_single_value_share,
            "max_covariate_r2": args.max_covariate_r2,
            "reject_negative_values": args.reject_negative_values,
            "trait_cap": args.trait_cap,
        },
        "traits": traits,
        "recommended_traits": recommended,
        "inputs": {
            "phenotype_sha256": sha256_file(args.phenotype),
            "covariates_sha256": sha256_file(args.covariates),
            "split_manifest_sha256": sha256_file(args.split_manifest),
        },
        "note": (
            "Recommendation is a screen on data quality, not on effect size. Freezing which "
            "traits enter qualification must happen before the test split is opened."
        ),
    }
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    lines = [
        "# Trait preflight",
        "",
        f"- 只读取 {sorted(preflight_splits)}；保持关闭：{closed or '无'}",
        f"- preflight 可用参与者 {genotyped_with_covariates}",
        "",
        "| trait | 非缺失 N | 缺失率 | 唯一值 | 众数占比 | >5SD | 协变量 R² | 可用 |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for item in traits:
        if "distribution" not in item:
            lines.append(f"| `{item['trait_column']}` | {item['non_missing']} | "
                         f"{item['missing_rate']} | - | - | - | - | ❌ |")
            continue
        lines.append(
            f"| `{item['trait_column']}` | {item['non_missing']} | {item['missing_rate']:.4f} | "
            f"{item['distinct_values']} | {item['most_frequent_value']['share']:.4f} | "  # type: ignore[index]
            f"{item['beyond_5sd']} | {item['covariate_model']['r2']:.4f} | "  # type: ignore[index]
            f"{'✅' if item['usable'] else '❌'} |"
        )
    lines += ["", "## 不可用原因", ""]
    flagged = [item for item in traits if not item["usable"]]
    if not flagged:
        lines.append("- 无")
    for item in flagged:
        lines.append(f"- `{item['trait_column']}`：" + "；".join(item["reasons"]))  # type: ignore[arg-type]
    lines += ["", f"## 建议冻结（上限 {args.trait_cap}）", "",
              ", ".join(f"`{name}`" for name in recommended) or "无", "",
              "> 这是数据质量筛选，不是效应量筛选。trait 冻结必须在开 test 之前完成。", ""]
    (out_dir / "PREFLIGHT.md").write_text("\n".join(lines), encoding="utf-8")

    print(json.dumps({
        "status": "OK",
        "report": str(report_path),
        "splits_read": sorted(preflight_splits),
        "splits_left_closed": closed,
        "recommended_traits": recommended,
        "usable_count": len(usable),
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
