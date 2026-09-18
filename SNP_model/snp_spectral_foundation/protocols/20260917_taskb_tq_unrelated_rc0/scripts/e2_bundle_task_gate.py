#!/usr/bin/env python3
"""Assemble the Task Gate audit bundle from artifacts already on disk.

The bundle lets someone re-derive the qualification numbers without the raw
cohort: the frozen config, what the splits contain, the per-threshold validation
curves with paired CIs, and optionally per-participant predictions under a
pseudonym.

Nothing here opens the test split. The bundle records that fact so an auditor can
check it rather than take it on trust.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_split import open_text, sha256_file  # noqa: E402

TRAIT_LABELS = {
    "n_30020_0_0": "HGB",
    "n_30010_0_0": "RBC",
    "n_30040_0_0": "MCV",
}


def write_split_summary(summary: dict, destination: Path) -> None:
    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["stratum", "split", "participants"])
        per_stratum = summary.get("split", {}).get("per_stratum_sample_counts", {})
        for stratum, counts in sorted(per_stratum.items()):
            for split, value in sorted(counts.items()):
                writer.writerow([stratum, split, value])
        for split, value in sorted(summary.get("split", {}).get("sample_counts", {}).items()):
            writer.writerow(["ALL", split, value])


def write_exclusions(summary: dict, config: dict, destination: Path) -> None:
    destination.write_text(json.dumps({
        "order_of_operations": [
            "1. covariates table restricted to participants aligned to the chr1 fam",
            "2. drop negative-sign identifiers (UKB's marker for withdrawn participants)",
            "3. drop participants absent from the kinship flag table",
            "4. keep only field 22021 == 0, dropping known-related and unknown-relatedness",
            "5. drop QC flags 22027 and 22019",
            "6. split the survivors, stratified, with a frozen seed",
            "7. per trait, drop participants whose trait value is missing",
        ],
        "counts": summary.get("filtering", {}),
        "relatedness_axis": summary.get("relatedness_axis", {}),
        "withdrawal": {
            "explicit_list_supplied": config.get("cohort", {}).get("withdrawal_list_supplied"),
            "partial_coverage": config.get("cohort", {}).get("withdrawal_partial_coverage"),
            "coverage_gap": config.get("cohort", {}).get("withdrawal_coverage_gap"),
            "resolution_required_before": config.get("cohort", {}).get(
                "withdrawal_resolution_required_before"
            ),
            "applied_at_step": 2,
            "note": (
                "No explicit withdrawal list exists on this server. Step 2 removes the "
                "negative-sign identifiers only, which predate the 2020-01 genotype "
                "conversion and therefore miss later withdrawals."
            ),
        },
    }, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def merge_predictions(sources: list[Path], destination: Path) -> int:
    fieldnames: list[str] = []
    rows: list[dict[str, str]] = []
    for source in sources:
        with open_text(source) as handle:
            for row in csv.DictReader(handle):
                for name in row:
                    if name not in fieldnames:
                        fieldnames.append(name)
                rows.append(row)
    if not rows:
        return 0
    rows.sort(key=lambda item: (item.get("trait", ""), item.get("surrogate_id", "")))
    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, restval="", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--split-summary", type=Path, required=True)
    parser.add_argument("--panel-summary", type=Path, required=True)
    parser.add_argument("--preflight", type=Path, default=None)
    parser.add_argument("--export-dir", type=Path, required=True,
                        help="Directory b5 --export-dir wrote into.")
    parser.add_argument("--trait", action="append", required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--include-predictions", action="store_true",
                        help="Include per-participant rows. They are individual-level "
                             "phenotype values under a pseudonym; confirm your data "
                             "agreement permits moving them before using this.")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    out_dir = args.out_dir.resolve()
    if out_dir.exists() and any(out_dir.iterdir()) and not args.overwrite:
        raise SystemExit(f"Refusing to overwrite a non-empty directory: {out_dir}")
    out_dir.mkdir(parents=True, exist_ok=True)

    config = json.loads(args.config.read_text(encoding="utf-8"))
    summary = json.loads(args.split_summary.read_text(encoding="utf-8"))
    panel = json.loads(args.panel_summary.read_text(encoding="utf-8"))

    shutil.copy2(args.config, out_dir / "CONFIG.json")
    write_split_summary(summary, out_dir / "split_manifest.csv")
    write_exclusions(summary, config, out_dir / "exclusions.json")

    verdicts: dict[str, object] = {}
    prediction_sources: list[Path] = []
    for trait in args.trait:
        label = TRAIT_LABELS.get(trait, trait)
        curve = args.export_dir / f"{trait}_validation_curve.csv"
        if not curve.exists():
            raise SystemExit(f"Missing curve for {trait}: {curve}")
        shutil.copy2(curve, out_dir / f"{label}_validation_curve.csv")
        counts = args.export_dir / f"{trait}_counts.json"
        if counts.exists():
            verdicts[trait] = json.loads(counts.read_text(encoding="utf-8"))
        predictions = args.export_dir / f"validation_predictions_{trait}.csv"
        if predictions.exists():
            prediction_sources.append(predictions)

    prediction_rows = 0
    if args.include_predictions and prediction_sources:
        prediction_rows = merge_predictions(
            prediction_sources, out_dir / "validation_predictions.csv"
        )

    (out_dir / "trait_counts.json").write_text(
        json.dumps(verdicts, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )

    panel_note = {
        "shape": panel.get("shape"),
        "chromosomes_included": config.get("panel", {}).get("chromosomes_included"),
        "chromosomes_excluded": config.get("panel", {}).get("chromosomes_excluded"),
        "dosage_definition": panel.get("dosage_definition"),
        "statistics_fitted_on": panel.get("statistics_fitted_on"),
        "filters": panel.get("filters"),
    }
    (out_dir / "panel_summary.json").write_text(
        json.dumps(panel_note, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )

    environment = {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
    }
    try:
        import numpy  # type: ignore

        environment["numpy"] = numpy.__version__
    except Exception:  # noqa: BLE001
        environment["numpy"] = None
    (out_dir / "environment.json").write_text(
        json.dumps(environment, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )

    with (out_dir / "hashes.txt").open("w", encoding="utf-8") as handle:
        for path in sorted(out_dir.iterdir()):
            if path.name == "hashes.txt" or path.is_dir():
                continue
            handle.write(f"{sha256_file(path)}  {path.name}\n")
        handle.write("\n# upstream artifacts, not copied into this bundle\n")
        for label, path in [
            ("split_manifest.tsv", args.split_summary.parent / "tq_split_manifest.tsv"),
            ("panel.int8.npy", args.panel_summary.parent / "panel.int8.npy"),
            ("panel_variants.tsv", args.panel_summary.parent / "panel_variants.tsv"),
        ]:
            if path.exists():
                handle.write(f"{sha256_file(path)}  {label}\n")

    readme = [
        "# Task Gate audit bundle",
        "",
        f"Traits: {', '.join(TRAIT_LABELS.get(t, t) for t in args.trait)}",
        "",
        "## What each file is",
        "",
        "| file | contents |",
        "| --- | --- |",
        "| `CONFIG.json` | the frozen protocol. `firewall.test_opened` is false and "
        "`ready_to_open_test` records that every binding field was set before any test access. |",
        "| `split_manifest.csv` | participant counts per stratum and split. Identifiers are not "
        "included; the manifest's own sha256 is in `hashes.txt`. |",
        "| `exclusions.json` | the order in which filters were applied, their counts, and the "
        "state of the withdrawal list. |",
        "| `*_validation_curve.csv` | one row per p-value threshold: variants selected, "
        "`r2_A`, `r2_B`, `delta_r2`, a paired bootstrap 95% CI, and Pearson correlations. "
        "`selected` marks the threshold chosen on validation. |",
        "| `trait_counts.json` | per trait: participants per split, the selected threshold, and "
        "the verdict, which is `VALIDATION_ONLY_TEST_NOT_OPENED`. |",
        "| `panel_summary.json` | panel shape, chromosome coverage and the filters that built it. |",
        "| `environment.json` | interpreter and numpy versions. |",
        "| `hashes.txt` | sha256 of every bundle file, plus the upstream split manifest, panel "
        "matrix and variant table that are too large to include. |",
        "",
        "## What an auditor can check here",
        "",
        "- Every `delta_r2` in the curves is computed on the validation split. The test split "
        "was never read: `trait_counts.json` carries the verdict and `CONFIG.json` carries "
        "`firewall.test_opened = false`.",
        "- The chosen threshold is the one flagged `selected`, and it was chosen by the largest "
        "validation `delta_r2`, which makes those numbers optimistic by construction.",
        "- Out-of-sample R2 uses the train mean as its reference, so it can go negative; the "
        "loose-threshold rows show this.",
        "",
        "## What this bundle cannot settle",
        "",
        "A negative `delta_r2` at loose thresholds is consistent with honest out-of-sample "
        "evaluation, but it is not proof that the implementation is leak-free. Confirming that "
        "needs the code path itself checked: that the split, the allele frequencies and "
        "missingness fills, the covariate coefficients, the per-variant effects and the "
        "threshold choice all use train and validation only. The scripts that produced these "
        "files are in the same package as this one.",
    ]
    if args.include_predictions:
        readme += [
            "",
            "## `validation_predictions.csv`",
            "",
            f"{prediction_rows} rows. Columns: `surrogate_id`, `trait`, `split`, `y_true_z`, "
            "`pred_A`, and the arm B prediction at the selected threshold.",
            "",
            "`surrogate_id` is `sha256(split_manifest_sha256 | sample_id)` truncated. The same "
            "participant keeps one id across the three traits, so the tables join. This is "
            "**pseudonymisation, not anonymisation**: anyone holding the underlying identifier "
            "list can recompute the mapping, and `y_true_z` is a participant's own phenotype "
            "value on a train-fitted z scale. Handle under the terms that govern the cohort.",
        ]
    (out_dir / "README.md").write_text("\n".join(readme) + "\n", encoding="utf-8")

    print(json.dumps({
        "status": "OK",
        "bundle": str(out_dir),
        "files": sorted(path.name for path in out_dir.iterdir()),
        "predictions_included": bool(args.include_predictions and prediction_rows),
        "prediction_rows": prediction_rows,
        "test_opened": config.get("firewall", {}).get("test_opened"),
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
