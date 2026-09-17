"""TQ-B1: measure the additive gene-effect bar. CPU only, numpy only.

Two phases, in this order and no other:

  1. Detectability gate. On the VALIDATION split, does adding a gene's cis
     dosage to the covariates raise held-out R2, averaged over genes? If not,
     the run stops and the evaluation split is never opened.
  2. The bar. On the EVALUATION split, perturb each gene's cis window, take the
     spread of the prediction change as that gene's effect score, and rank it
     against the independent WES burden effect.

The resulting Spearman is the number a foundation model has to beat later.

    python run_bar.py --binding BINDING.json --out results
"""

from __future__ import annotations

import argparse
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from tqb1_core import (  # noqa: E402
    PrimalRidge,
    af_missingness_dosage,
    clustered_macro_bootstrap_shared,
    decile_strata,
    deterministic_order,
    gene_bootstrap_spearman,
    matched_null_spearman,
    perturbation_score,
    r2_against_train_mean,
    read_bed_variants,
    spearman,
    uniform_thin,
)
from tqb1_io import (  # noqa: E402
    assign_splits,
    check_bed_size,
    dump_json,
    load_json,
    numeric_table,
    read_bim,
    read_burden,
    read_fam,
    read_genes,
    read_samples,
    sha256_file,
)


def full(value: float) -> str:
    """Full round-trip precision. Incident 01 of TQ-G1 rc0 was 6-decimal output."""
    return repr(float(value))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binding", default="BINDING.json")
    parser.add_argument("--config", default="CONFIG.json")
    parser.add_argument("--out", default="results")
    parser.add_argument("--gate-only", action="store_true", help="Stop after the detectability gate.")
    parser.add_argument("--max-genes", type=int, default=None, help="Override for a quick trial run.")
    args = parser.parse_args()

    config = load_json(HERE / args.config)
    binding = load_json(Path(args.binding))
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    started = time.time()

    # ---- bind -------------------------------------------------------------
    geno = binding["genotype"]
    bed, bim_path, fam_path = Path(geno["bed"]), Path(geno["bim"]), Path(geno["fam"])
    sample_ids_geno = read_fam(fam_path)
    chroms, positions, variant_ids = read_bim(bim_path)
    check_bed_size(bed, len(sample_ids_geno), len(positions))
    print(f"genotype: {len(sample_ids_geno)} samples x {len(positions)} variants", flush=True)

    ids, groups, supplied_splits = read_samples(Path(binding["samples_path"]))
    split_labels = assign_splits(
        ids, groups, supplied_splits, config["split_fractions"], int(config["seed"])
    )
    geno_row = {s: i for i, s in enumerate(sample_ids_geno)}
    missing = [s for s in ids if s not in geno_row]
    if missing:
        raise SystemExit(f"{len(missing)} analysed samples are absent from the .fam, e.g. {missing[:3]}")
    rows = np.asarray([geno_row[s] for s in ids], dtype=np.int64)

    cov_ids, cov_names, cov_values = numeric_table(Path(binding["covariates_path"]))
    phe_ids, trait_names, phe_values = numeric_table(Path(binding["phenotypes_path"]))
    if cov_ids != ids or phe_ids != ids:
        raise SystemExit("covariates and phenotypes must have exactly the sample order of the samples file")
    traits = list(binding["traits"]) or trait_names
    for trait in traits:
        if trait not in trait_names:
            raise SystemExit(f"trait '{trait}' is not a column of the phenotype file")

    train = np.flatnonzero(split_labels == "train")
    validation = np.flatnonzero(split_labels == "validation")
    evaluation = np.flatnonzero(split_labels == "evaluation")
    for name, idx in (("train", train), ("validation", validation), ("evaluation", evaluation)):
        if len(idx) < 2:
            raise SystemExit(f"split '{name}' has {len(idx)} individuals")
    group_of = np.asarray(groups)
    if set(group_of[train]) & set(group_of[validation]) or set(group_of[train]) & set(group_of[evaluation]):
        raise SystemExit("group leakage across splits")
    if set(group_of[validation]) & set(group_of[evaluation]):
        raise SystemExit("group leakage between validation and evaluation")
    print(f"splits: train={len(train)} validation={len(validation)} evaluation={len(evaluation)}", flush=True)

    # ---- eligible genes ---------------------------------------------------
    genes = read_genes(Path(binding["genes_path"]))
    wanted_chroms = set(binding.get("chromosomes") or [])
    if wanted_chroms:
        genes = [g for g in genes if g["chrom"] in wanted_chroms]
    chrom_array = np.asarray(chroms)
    by_chrom = {c: np.flatnonzero(chrom_array == c) for c in {g["chrom"] for g in genes}}
    radius = int(config["cis_radius_bp"])
    snps_per_gene = int(config["snps_per_gene"])

    eligible: dict[str, dict[str, object]] = {}
    for gene in genes:
        pool = by_chrom.get(gene["chrom"])
        if pool is None or len(pool) == 0:
            continue
        window = pool[np.abs(positions[pool] - int(gene["tss"])) <= radius]
        if len(window) < snps_per_gene:
            continue
        eligible[str(gene["gene_id"])] = {**gene, "window": window}
    ordered = deterministic_order(sorted(eligible), str(config["gene_selection_salt"]))
    cap = args.max_genes or int(config["max_genes"])
    selected = sorted(ordered[:cap])
    minimum = int(config["min_genes_required"]) if args.max_genes is None else 1
    if len(selected) < minimum:
        raise SystemExit(f"DESIGN-INELIGIBLE: {len(selected)} eligible genes, {minimum} required")
    print(f"genes: {len(eligible)} eligible, {len(selected)} analysed", flush=True)

    alphas = [float(a) for a in config["ridge_alphas"]]
    maf_min, missing_max = float(config["maf_min"]), float(config["missingness_max"])

    covariates = cov_values
    n_cov = covariates.shape[1]
    # Both arms treat covariates identically and leave them unpenalised, so the
    # only thing alpha shrinks is the genetic block.
    mask_a = np.zeros(n_cov)
    mask_b = np.concatenate([np.zeros(n_cov), np.ones(snps_per_gene)])
    a_models: dict[str, tuple[PrimalRidge, float]] = {}
    for trait in traits:
        y = phe_values[:, trait_names.index(trait)]
        model = PrimalRidge(covariates[train], y[train], mask_a)
        alpha, _ = model.select_alpha(covariates[validation], y[validation], alphas)
        a_models[trait] = (model, alpha)
    print(f"arm A fitted for {len(traits)} traits on {len(cov_names)} covariates", flush=True)

    # ---- per-gene pass ----------------------------------------------------
    n_selected = len(selected)
    pred_val = {t: {"A": np.zeros((n_selected, len(validation))), "B": np.zeros((n_selected, len(validation)))} for t in traits}
    pred_eval = {t: {"A": np.zeros((n_selected, len(evaluation))), "B": np.zeros((n_selected, len(evaluation)))} for t in traits}
    scores = {t: np.zeros(n_selected) for t in traits}
    gene_meta: list[dict[str, object]] = []

    for g, gene_id in enumerate(selected):
        item = eligible[gene_id]
        window = np.asarray(item["window"], dtype=np.int64)
        raw = read_bed_variants(bed, len(sample_ids_geno), window, len(positions))[rows]
        af, miss, _ = af_missingness_dosage(raw, train)
        maf = np.minimum(af, 1.0 - af)
        keep = np.flatnonzero((maf >= maf_min) & (miss <= missing_max) & np.isfinite(af))
        if len(keep) < snps_per_gene:
            gene_meta.append({"gene_id": gene_id, "skipped": "insufficient_qc_variants"})
            continue
        chosen = uniform_thin(keep, snps_per_gene)
        af_c, _, dosage = af_missingness_dosage(raw[:, chosen], train)
        design = np.column_stack([covariates, dosage])
        perturbed = np.column_stack(
            [covariates, np.broadcast_to(dosage[train].mean(axis=0), dosage.shape)]
        )
        for trait in traits:
            y = phe_values[:, trait_names.index(trait)]
            a_model, a_alpha = a_models[trait]
            b_model = PrimalRidge(design[train], y[train], mask_b)
            b_alpha, _ = b_model.select_alpha(design[validation], y[validation], alphas)
            pred_val[trait]["A"][g] = a_model.predict(covariates[validation], a_alpha)
            pred_val[trait]["B"][g] = b_model.predict(design[validation], b_alpha)
            pred_eval[trait]["A"][g] = a_model.predict(covariates[evaluation], a_alpha)
            pred_eval[trait]["B"][g] = b_model.predict(design[evaluation], b_alpha)
            _, scores[trait][g] = perturbation_score(
                b_model, design[evaluation], perturbed[evaluation], b_alpha
            )
        gene_meta.append(
            {
                "gene_id": gene_id,
                "chrom": item["chrom"],
                "tss": int(item["tss"]),
                "cis_qc_variant_count": int(len(keep)),
                "mean_maf": float(np.mean(np.minimum(af_c, 1.0 - af_c))),
            }
        )
        if (g + 1) % 50 == 0 or g + 1 == n_selected:
            rate = (g + 1) / max(1e-9, time.time() - started)
            print(f"  gene {g + 1}/{n_selected}  ({rate:.1f}/s)", flush=True)

    usable = np.asarray([i for i, m in enumerate(gene_meta) if "skipped" not in m], dtype=np.int64)
    if len(usable) < minimum:
        raise SystemExit(f"only {len(usable)} genes survived QC, {minimum} required")

    # ---- phase 1: detectability gate on VALIDATION -------------------------
    gate: dict[str, object] = {}
    for trait in traits:
        y = phe_values[:, trait_names.index(trait)]
        gate[trait] = clustered_macro_bootstrap_shared(
            y[validation],
            {arm: pred_val[trait][arm][usable] for arm in ("A", "B")},
            float(y[train].mean()),
            ("B", "A"),
            int(config["bootstrap_seed"]),
            int(config["bootstrap_replicates"]),
        )
        gate[trait].pop("per_gene_r2", None)
        d = gate[trait]
        print(
            f"gate {trait}: B-A delta={d['macro_delta_r2']:+.6f} "
            f"CI={[round(v, 6) for v in d['paired_bootstrap_ci95']]} "
            f"{'PASS' if d['pass'] else 'FAIL'}",
            flush=True,
        )
    passing = [t for t in traits if gate[t]["pass"]]
    dump_json(out_dir / "GATE.json", {"rule": config["detectability_gate"]["rule"], "per_trait": gate,
                                      "passing_traits": passing})

    if not passing or args.gate_only:
        dump_json(
            out_dir / "FINAL_STATUS.json",
            {
                "analysis_id": config["analysis_id"],
                "status": "GATE-ONLY" if args.gate_only else "TASK-INELIGIBLE",
                "passing_traits": passing,
                "evaluation_split_opened": False,
                "n_genes": int(len(usable)),
                "elapsed_seconds": round(time.time() - started, 1),
            },
        )
        print("\nevaluation split NOT opened." if not passing else "\ngate-only run complete.", flush=True)
        return 0 if passing else 2

    # ---- phase 2: the bar on EVALUATION ------------------------------------
    burden_path = binding.get("burden_path")
    bar: dict[str, object] = {}
    held_out: dict[str, object] = {}
    burden = read_burden(Path(burden_path)) if burden_path else {}
    eval_r2: dict[str, dict[str, np.ndarray]] = {}

    for trait in passing:
        y = phe_values[:, trait_names.index(trait)]
        summary = clustered_macro_bootstrap_shared(
            y[evaluation],
            {arm: pred_eval[trait][arm][usable] for arm in ("A", "B")},
            float(y[train].mean()),
            ("B", "A"),
            int(config["bootstrap_seed"]) + 1,
            int(config["bootstrap_replicates"]),
        )
        per_gene_r2 = summary.pop("per_gene_r2")
        eval_r2[trait] = per_gene_r2
        held_out[trait] = summary
        if not burden:
            continue
        have = [i for i in range(len(usable)) if (gene_meta[usable[i]]["gene_id"], trait) in burden]
        if len(have) < int(config["min_genes_required"]):
            bar[trait] = {"status": "INSUFFICIENT_BURDEN_OVERLAP", "n_genes_with_burden": len(have)}
            continue
        idx = np.asarray(have, dtype=np.int64)
        s = scores[trait][usable][idx]
        beta = np.asarray([abs(burden[(gene_meta[usable[i]]["gene_id"], trait)][0]) for i in have])
        se = np.asarray([burden[(gene_meta[usable[i]]["gene_id"], trait)][1] for i in have])
        strata_inputs = np.column_stack(
            [
                decile_strata(np.asarray([gene_meta[usable[i]]["cis_qc_variant_count"] for i in have], dtype=float)),
                decile_strata(np.asarray([gene_meta[usable[i]]["mean_maf"] for i in have], dtype=float)),
                decile_strata(se),
            ]
        )
        _, strata = np.unique(strata_inputs, axis=0, return_inverse=True)
        bar[trait] = {
            "n_genes_with_burden": int(len(have)),
            "matched_null": matched_null_spearman(
                s, beta, strata.astype(np.int64),
                int(config["permutation_seed"]), int(config["permutation_replicates"])
            ),
            "gene_bootstrap": gene_bootstrap_spearman(
                s, beta, int(config["bootstrap_seed"]) + 2, int(config["bootstrap_replicates"])
            ),
            "spearman_score_vs_heldout_cis_r2": spearman(s, np.asarray(per_gene_r2["B"])[idx]),
        }
        b = bar[trait]
        print(
            f"BAR {trait}: spearman={b['matched_null']['observed_spearman']:+.4f} "
            f"CI={[round(v, 4) for v in b['gene_bootstrap']['bootstrap_ci95']]} "
            f"matched-null p={b['matched_null']['p_one_sided_greater']:.4f} "
            f"({b['n_genes_with_burden']} genes)",
            flush=True,
        )

    # ---- artifacts ---------------------------------------------------------
    with (out_dir / "GENE_TABLE.tsv").open("w", encoding="utf-8", newline="") as handle:
        header = ["gene_id", "chrom", "tss", "cis_qc_variant_count", "mean_maf"]
        for trait in passing:
            header += [f"score_{trait}", f"r2_A_eval_{trait}", f"r2_B_eval_{trait}"]
        handle.write("\t".join(header) + "\n")
        for i, k in enumerate(usable):
            m = gene_meta[k]
            row = [str(m["gene_id"]), str(m["chrom"]), str(m["tss"]),
                   str(m["cis_qc_variant_count"]), full(m["mean_maf"])]
            for trait in passing:
                row += [
                    full(scores[trait][k]),
                    full(eval_r2[trait]["A"][i]),
                    full(eval_r2[trait]["B"][i]),
                ]
            handle.write("\t".join(row) + "\n")

    results = {
        "analysis_id": config["analysis_id"],
        "claim_class": config["claim_class"],
        "status": "VALID_COMPLETE",
        "splits": {"train": len(train), "validation": len(validation), "evaluation": len(evaluation),
                   "unit": config["split_unit"], "group_disjoint": True},
        "genes": {"eligible": len(eligible), "analysed": int(len(usable)),
                  "selection_rule": config["gene_selection_rule"],
                  "selection_uses_phenotype": False},
        "covariates": {"count": len(cov_names), "names": cov_names},
        "detectability_gate": gate,
        "passing_traits": passing,
        "held_out_macro": held_out,
        "bar": bar,
        "interpretation_limits": config["interpretation_limits"],
        "elapsed_seconds": round(time.time() - started, 1),
    }
    dump_json(out_dir / "RESULTS.json", results)
    dump_json(
        out_dir / "RUN_BINDING.json",
        {
            "inputs": {
                name: {"path": str(p), "sha256": sha256_file(Path(p))}
                for name, p in [
                    ("bed", geno["bed"]), ("bim", geno["bim"]), ("fam", geno["fam"]),
                    ("samples", binding["samples_path"]),
                    ("covariates", binding["covariates_path"]),
                    ("phenotypes", binding["phenotypes_path"]),
                    ("genes", binding["genes_path"]),
                ] + ([("burden", burden_path)] if burden_path else [])
            },
            "config_sha256": sha256_file(HERE / args.config),
            "core_sha256": sha256_file(HERE / "tqb1_core.py"),
            "runner_sha256": sha256_file(Path(__file__).resolve()),
            "runtime": {"python": platform.python_version(), "platform": platform.platform(),
                        "numpy": np.__version__, "gpu_used": False},
        },
    )
    dump_json(
        out_dir / "FINAL_STATUS.json",
        {
            "analysis_id": config["analysis_id"],
            "status": "VALID_COMPLETE",
            "passing_traits": passing,
            "evaluation_split_opened": True,
            "n_genes": int(len(usable)),
            "bar_spearman": {t: (bar[t]["matched_null"]["observed_spearman"]
                                 if "matched_null" in bar.get(t, {}) else None) for t in passing},
            "elapsed_seconds": round(time.time() - started, 1),
        },
    )
    print(f"\nTQ-B1 complete in {time.time() - started:.0f}s. Artifacts in {out_dir}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
