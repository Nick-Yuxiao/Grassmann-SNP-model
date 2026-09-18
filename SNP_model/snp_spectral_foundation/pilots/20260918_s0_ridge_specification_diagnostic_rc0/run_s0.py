"""S0: does penalising the covariate block explain ``macro R2(B) < macro R2(A)``?

Development individuals only.  ``task_gate`` and ``bridge_test`` outcomes are
never loaded, so this cannot consume a decision resource and cannot overturn a
frozen pilot.

The gene panel, the cis window, the uniform thinning, the QC, the genotype PCs
and the covariate block are all inherited from TQ-G1 by importing its loaders,
so the fit specification is the only thing that varies between arms.

    python -m unittest test_s0 -v
    python smoke_s0.py
    python validate_s0.py --results results_smoke
    python run_s0.py
"""

from __future__ import annotations

import json
import platform
import sys
from pathlib import Path

import numpy as np
import sklearn

PILOT_DIR = Path(__file__).resolve().parent
TQG1_DIR = PILOT_DIR.parent / "20260917_tqg1_gene_layer_readout_rc0"
if str(TQG1_DIR) not in sys.path:
    sys.path.insert(0, str(TQG1_DIR))

import run_tqg1 as T  # noqa: E402
from tqg1_core import (  # noqa: E402
    af_missingness_dosage,
    deterministic_gene_order,
    two_way_clustered_bootstrap,
    uniform_thin,
)

from s0_core import SPECIFICATIONS, cross_validated_predictions, group_kfold  # noqa: E402

CONFIG_PATH = PILOT_DIR / "CONFIG_S0.json"
BINDING_PATH = PILOT_DIR / "PRE_RUN_BINDING.json"
RESULT_DIR = PILOT_DIR / "results"


def resolve_development_rows(
    rows: list[dict[str, str]], config: dict[str, object]
) -> tuple[np.ndarray, list[str]]:
    """Bind development rows and refuse every forbidden role."""
    used = set(config["roles_used"])
    forbidden = set(config["forbidden_roles"])
    dev_idx = np.asarray(
        [i for i, r in enumerate(rows) if r["development_role"] in used], dtype=np.int64
    )
    if len(dev_idx) == 0:
        raise SystemExit("no development individuals in the manifest")
    leaked = sorted({rows[i]["primary_role"] for i in dev_idx} & forbidden)
    if leaked:
        raise SystemExit(f"forbidden roles reached the development set: {leaked}")
    train_idx = np.asarray(
        [i for i, r in enumerate(rows) if r["development_role"] == "dev_train"], dtype=np.int64
    )
    if len(train_idx) == 0:
        raise SystemExit("empty dev_train")
    families = [rows[i]["family"] for i in dev_idx]
    print(
        f"development bound: {len(dev_idx)} individuals, {len(set(families))} families "
        f"(dev_train={len(train_idx)}); sealed roles untouched: {sorted(forbidden)}",
        flush=True,
    )
    return dev_idx, families


def analyse(
    config: dict[str, object],
    rows: list[dict[str, str]],
    positions: np.ndarray,
    genotype: np.ndarray,
    load_records,
    dev_idx: np.ndarray,
    families: list[str],
) -> dict[str, object]:
    dev_rows = [rows[i] for i in dev_idx]
    local_train = np.flatnonzero(
        np.asarray([r["development_role"] == "dev_train" for r in dev_rows])
    )
    dev_genotype = genotype[dev_idx]

    af, missingness, dosage_all = af_missingness_dosage(dev_genotype, local_train)
    maf = np.minimum(af, 1.0 - af)
    qc = np.flatnonzero(
        (maf >= float(config["maf_min"]))
        & (missingness <= float(config["missingness_max"]))
        & np.isfinite(af)
    )
    print(f"QC variants on dev_train: {len(qc)} of {len(positions)}", flush=True)

    pc_panel = uniform_thin(qc, min(int(config["pc_panel_snps"]), len(qc)))
    pcs, pc_ratio = T.genotype_pcs(
        dosage_all[:, pc_panel], local_train, int(config["genotype_pc_count"])
    )
    covariates = np.column_stack([T.demographic_covariates(dev_rows), pcs]).astype(np.float64)

    records = load_records(dev_idx)
    radius = int(config["cis_radius_bp"])
    snps_per_gene = int(config["snps_per_gene"])
    eligible: dict[str, dict[str, object]] = {}
    for record in records:
        values = np.asarray(record["values"], dtype=np.float64)
        if not np.isfinite(values).all() or float(np.var(values)) <= 0:
            continue
        tss = int(record["tss_1based"])
        window = qc[(positions[qc] >= tss - radius) & (positions[qc] <= tss + radius)]
        if len(window) < snps_per_gene:
            continue
        eligible[str(record["gene_id"])] = {
            "tss_1based": tss,
            "variant_indices": uniform_thin(window, snps_per_gene),
            "cis_qc_variant_count": int(len(window)),
            "values": values,
        }
    ordered = deterministic_gene_order(sorted(eligible), str(config["gene_selection_salt"]))
    selected = sorted(ordered[: int(config["max_genes"])])
    if len(selected) < int(config["min_genes_required"]):
        raise SystemExit(
            f"DESIGN-INELIGIBLE: only {len(selected)} eligible genes, "
            f"{config['min_genes_required']} required"
        )
    print(f"genes: {len(eligible)} eligible, {len(selected)} analysed", flush=True)

    folds = group_kfold(families, int(config["cv_folds"]), int(config["cv_fold_seed"]))
    alphas = [float(a) for a in config["ridge_alphas"]]
    inner_seed = int(config["inner_split_seed"])
    counts = [int(c) for c in config["snp_counts"]]
    n_dev = len(dev_idx)

    truth = np.zeros((len(selected), n_dev))
    means = np.zeros(len(selected))
    arm_keys = [f"A@{s}" for s in SPECIFICATIONS]
    arm_keys += [f"B@{s}:{c}" for s in SPECIFICATIONS for c in counts]
    predictions = {key: np.zeros((len(selected), n_dev)) for key in arm_keys}
    alpha_log: dict[str, list[float]] = {key: [] for key in arm_keys}

    for g, gene in enumerate(selected):
        item = eligible[gene]
        y = np.asarray(item["values"], dtype=np.float64)
        truth[g] = y
        means[g] = float(y.mean())
        panel = np.asarray(item["variant_indices"], dtype=np.int64)
        for specification in SPECIFICATIONS:
            key = f"A@{specification}"
            fitted, chosen = cross_validated_predictions(
                specification, covariates, None, y, folds, families, alphas, inner_seed
            )
            predictions[key][g] = fitted
            alpha_log[key].extend(chosen)
            for count in counts:
                subset = uniform_thin(panel, count)
                dosage = dosage_all[:, subset].astype(np.float64)
                key = f"B@{specification}:{count}"
                fitted, chosen = cross_validated_predictions(
                    specification, covariates, dosage, y, folds, families, alphas, inner_seed
                )
                predictions[key][g] = fitted
                alpha_log[key].extend(chosen)
        if (g + 1) % 20 == 0 or g + 1 == len(selected):
            print(f"  gene {g + 1}/{len(selected)} done", flush=True)

    macro = {
        key: float(
            np.mean(
                1.0
                - np.sum((truth - value) ** 2, axis=1)
                / np.sum((truth - means[:, None]) ** 2, axis=1)
            )
        )
        for key, value in predictions.items()
    }
    contrasts: dict[str, object] = {}
    for specification in SPECIFICATIONS:
        for count in counts:
            name = f"B-A@{specification}:{count}"
            contrasts[name] = two_way_clustered_bootstrap(
                truth,
                {
                    "A": predictions[f"A@{specification}"],
                    "B": predictions[f"B@{specification}:{count}"],
                },
                means,
                ("B", "A"),
                int(config["bootstrap_seed"]),
                int(config["bootstrap_replicates"]),
            )

    return {
        "analysis_id": config["analysis_id"],
        "pilot_id": config["pilot_id"],
        "claim_class": config["claim_class"],
        "decides_nothing_about": config["decides_nothing_about"],
        "development": {
            "n_individuals": int(n_dev),
            "n_families": int(len(set(families))),
            "n_dev_train": int(len(local_train)),
            "cv_folds": int(config["cv_folds"]),
            "fold_sizes": [int((folds == f).sum()) for f in sorted(set(folds.tolist()))],
            "roles_used": config["roles_used"],
            "sealed_roles_never_read": config["forbidden_roles"],
        },
        "qc": {"scope": "dev_train_only", "variants_passing": int(len(qc))},
        "covariates": {
            "columns": int(covariates.shape[1]),
            "genotype_pc_count": int(config["genotype_pc_count"]),
            "pc_explained_variance_ratio": pc_ratio,
        },
        "genes": {
            "eligible": int(len(eligible)),
            "analysed": int(len(selected)),
            "selection_uses_phenotype_signal": False,
            "gene_ids": selected,
        },
        "snp_counts": counts,
        "specifications": list(SPECIFICATIONS),
        "macro_r2": macro,
        "contrasts": contrasts,
        "chosen_alpha_median": {k: float(np.median(v)) for k, v in alpha_log.items()},
        "interpretation_limits": config["interpretation_limits"],
    }


def verdict(results: dict[str, object], config: dict[str, object]) -> dict[str, object]:
    primary = int(config["snps_per_gene"])
    rows = {}
    for specification in SPECIFICATIONS:
        entry = results["contrasts"][f"B-A@{specification}:{primary}"]
        low, high = entry["paired_bootstrap_ci95"]
        rows[specification] = {
            "delta": float(entry["macro_delta_r2"]),
            "ci95": [float(low), float(high)],
            "clears_zero": bool(entry["pass"]),
        }
    artefact = (
        rows["S2"]["clears_zero"]
        and rows["S3"]["clears_zero"]
        and not rows["S1"]["clears_zero"]
    )
    return {
        "analysis_id": config["analysis_id"],
        "primary_snp_count": primary,
        "per_specification": rows,
        "status": "SPECIFICATION-ARTEFACT-SUPPORTED" if artefact else "SPECIFICATION-ARTEFACT-NOT-SUPPORTED",
        "overturns_tqg1": False,
        "sealed_roles_never_read": config["forbidden_roles"],
    }


def main() -> None:
    config = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
    binding = json.loads(BINDING_PATH.read_text(encoding="utf-8"))
    if (RESULT_DIR / "FINAL_STATUS_S0.json").exists():
        raise SystemExit("this diagnostic already has a final status; move results/ aside to rerun")

    paths = T.resolve_assets(binding)
    audit = T.audit_assets(paths, binding)
    manifest_path = (PILOT_DIR / str(binding["sample_manifest_path"])).resolve()
    manifest_sha = T.sha256_file(manifest_path)
    if manifest_sha != binding["sample_manifest_expected_sha256"]:
        raise SystemExit("SAMPLE_MANIFEST.tsv hash mismatch: the frozen split changed")
    rows = T.load_tsv(manifest_path)

    T.set_reproducibility(int(config["seed"]))
    positions, _ = T.read_pvar(
        paths["pvar"], str(config["chromosome"]), int(binding["expected_variant_count"])
    )
    sample_ids = [r["sample_id"] for r in rows]
    genotype = T.read_pgen(paths["pgen"], len(sample_ids), len(positions))
    if genotype.shape[0] != len(sample_ids):
        raise SystemExit("manifest/PSAM sample order mismatch")

    dev_idx, families = resolve_development_rows(rows, config)
    dev_sample_ids = [sample_ids[i] for i in dev_idx]

    def load_records(idx: np.ndarray):
        if not np.array_equal(np.asarray(idx), dev_idx):
            raise SystemExit("expression may only be read for the development rows")
        return T.load_expression(paths["expression"], str(config["chromosome"]), dev_sample_ids)

    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    results = analyse(config, rows, positions, genotype, load_records, dev_idx, families)
    status = verdict(results, config)

    T.dump_json(RESULT_DIR / "RESULTS_S0.json", results)
    T.dump_json(RESULT_DIR / "FINAL_STATUS_S0.json", status)
    T.dump_json(
        RESULT_DIR / "RUN_BINDING_S0.json",
        {
            "assets": audit,
            "sample_manifest_sha256": manifest_sha,
            "config_sha256": T.sha256_file(CONFIG_PATH),
            "core_sha256": T.sha256_file(PILOT_DIR / "s0_core.py"),
            "runtime": {
                "python": platform.python_version(),
                "platform": platform.platform(),
                "numpy": np.__version__,
                "sklearn": sklearn.__version__,
            },
        },
    )
    print("\n--- S0 summary ---")
    for specification, row in status["per_specification"].items():
        print(
            "  %-3s B-A=%+.6f CI=[%+.6f, %+.6f] clears_zero=%s"
            % (specification, row["delta"], row["ci95"][0], row["ci95"][1], row["clears_zero"])
        )
    print("  status:", status["status"])
    print("artifacts written to  :", RESULT_DIR)


if __name__ == "__main__":
    main()
