"""Synthetic end-to-end smoke run for TQ-G1.

No GEUVADIS asset is touched. It builds a small LD-structured chromosome with
real cis signal, then drives the identical ``run_analysis`` the frozen runner
uses, so the whole chain can be checked before any real outcome is opened.

    python smoke_tqg1.py

It writes to results_smoke/ and proves nothing about the science.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

PILOT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(PILOT_DIR))

from run_tqg1 import resolve_roles, run_analysis, set_reproducibility  # noqa: E402

SMOKE_DIR = PILOT_DIR / "results_smoke"

N_TRAIN, N_VALIDATION, N_EVAL = 60, 20, 40
N_VARIANTS = 12_000
SPAN_BP = 20_000_000
N_GENES = 24
POPULATIONS = ["CEU", "FIN", "GBR", "TSI", "YRI"]


def smoke_config() -> dict[str, object]:
    config = json.loads((PILOT_DIR / "CONFIG_FROZEN.json").read_text(encoding="utf-8"))
    config.update(
        {
            "snps_per_gene": 64,
            "pretrain_blocks": 16,
            "pretrain_epochs": 3,
            "pc_panel_snps": 2_000,
            "genotype_pc_count": 5,
            "max_genes": N_GENES,
            "min_genes_required": 5,
            "bootstrap_replicates": 200,
            "permutation_replicates": 200,
            "evaluation_role": "task_gate",
        }
    )
    return config


def synthetic_manifest(rng: np.random.Generator) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    total = N_TRAIN + N_VALIDATION + N_EVAL
    for i in range(total):
        if i < N_TRAIN:
            primary, development = "development", "dev_train"
        elif i < N_TRAIN + N_VALIDATION:
            primary, development = "development", "dev_validation"
        else:
            primary, development = "task_gate", ""
        rows.append(
            {
                "sample_id": f"SM{i:05d}",
                "family": f"FAM{i:05d}",
                "pop": POPULATIONS[i % len(POPULATIONS)],
                "super_pop": "EUR",
                "sex": "2" if rng.random() < 0.5 else "1",
                "primary_role": primary,
                "development_role": development,
            }
        )
    return rows


def synthetic_genotypes(rng: np.random.Generator, n: int) -> tuple[np.ndarray, np.ndarray]:
    """A Markov haplotype chain gives short-range LD for the encoder to learn."""
    positions = np.sort(rng.choice(SPAN_BP, size=N_VARIANTS, replace=False)).astype(np.int64)
    frequency = rng.uniform(0.08, 0.92, size=N_VARIANTS)
    switch = 0.12  # per-marker recombination-like probability
    haplotypes = np.zeros((2 * n, N_VARIANTS), dtype=np.int8)
    current = (rng.random(2 * n) < frequency[0]).astype(np.int8)
    haplotypes[:, 0] = current
    for j in range(1, N_VARIANTS):
        fresh = (rng.random(2 * n) < frequency[j]).astype(np.int8)
        keep = rng.random(2 * n) >= switch
        current = np.where(keep, current, fresh).astype(np.int8)
        haplotypes[:, j] = current
    genotype = (haplotypes[0::2] + haplotypes[1::2]).astype(np.int8)
    missing = rng.random(genotype.shape) < 0.002
    genotype[missing] = -1
    return positions, genotype


def synthetic_expression(
    rng: np.random.Generator, positions: np.ndarray, genotype: np.ndarray, rows: list[dict[str, str]]
):
    """Each gene gets a handful of causal cis SNPs plus sex and noise."""
    tss_points = np.linspace(1_500_000, SPAN_BP - 1_500_000, N_GENES).astype(np.int64)
    dosage = np.where(genotype < 0, 0.0, genotype).astype(np.float64)
    dosage = (dosage - dosage.mean(axis=0)) / np.maximum(dosage.std(axis=0), 1e-8)
    sex = np.asarray([1.0 if r["sex"] == "2" else 0.0 for r in rows])
    records = []
    for g, tss in enumerate(tss_points):
        window = np.flatnonzero(np.abs(positions - tss) <= 1_000_000)
        causal = rng.choice(window, size=6, replace=False)
        beta = rng.normal(scale=0.55, size=len(causal))
        signal = dosage[:, causal] @ beta
        values = signal + 0.4 * sex + rng.normal(scale=1.0, size=len(rows))
        records.append(
            {"gene_id": f"ENSGSMOKE{g:05d}", "tss_1based": int(tss), "values": values}
        )
    return records


def main() -> None:
    config = smoke_config()
    set_reproducibility(int(config["seed"]))
    rng = np.random.default_rng(int(config["seed"]))

    rows = synthetic_manifest(rng)
    positions, genotype = synthetic_genotypes(rng, len(rows))
    records = synthetic_expression(rng, positions, genotype, rows)
    print(f"synthetic chromosome: {genotype.shape[0]} people x {genotype.shape[1]} variants", flush=True)

    train_idx, val_idx, eval_idx = resolve_roles(rows, str(config["evaluation_role"]))

    def load_records(outcome_idx: np.ndarray) -> list[dict[str, object]]:
        return [
            {
                "gene_id": r["gene_id"],
                "tss_1based": r["tss_1based"],
                "values": np.asarray(r["values"])[outcome_idx],
            }
            for r in records
        ]

    SMOKE_DIR.mkdir(parents=True, exist_ok=True)
    results = run_analysis(
        config, rows, positions, genotype, load_records, train_idx, val_idx, eval_idx, SMOKE_DIR
    )

    primary = results["primary"]
    rehearsal = results["gene_effect_rehearsal"]
    print("\n--- smoke summary (synthetic; no scientific meaning) ---")
    print(f"genes analysed        : {results['genes']['analysed']}")
    print(f"macro R2              : " + ", ".join(
        f"{arm}={primary['macro_r2'][arm]:+.4f}" for arm in ("A", "B", "C_gene", "C_full")
    ))
    print(f"primary {primary['contrast']:<14}: delta={primary['macro_delta_r2']:+.4f} "
          f"CI95={[round(v, 4) for v in primary['paired_bootstrap_ci95']]} "
          f"{'PASS' if primary['pass'] else 'FAIL'}")
    for name, value in results["secondary"].items():
        print(f"secondary {name:<12}: delta={value['macro_delta_r2']:+.4f} "
              f"CI95={[round(v, 4) for v in value['paired_bootstrap_ci95']]}")
    print(f"rehearsal spearman    : model={rehearsal['model']['observed_spearman']:+.3f} "
          f"linear={rehearsal['linear']['observed_spearman']:+.3f} "
          f"difference={rehearsal['model_minus_linear']['observed_difference']:+.3f} "
          f"{'PASS' if rehearsal['model_minus_linear']['pass'] else 'FAIL'}")
    print(f"\nartifacts written to  : {SMOKE_DIR}")


if __name__ == "__main__":
    main()
