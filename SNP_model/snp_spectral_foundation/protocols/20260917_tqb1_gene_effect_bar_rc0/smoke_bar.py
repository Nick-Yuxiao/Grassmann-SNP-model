"""Synthetic end-to-end smoke run for TQ-B1.

Writes a small but genuine PLINK 1 .bed/.bim/.fam plus the five TSV inputs,
then runs the real binder and the real runner over them. Nothing here touches
a cohort. Use it to prove the chain works on a server before any real data is
bound, and to confirm the machine produces the expected numbers.

    python smoke_bar.py
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from tqb1_core import bed_bytes_per_variant  # noqa: E402

OUT = HERE / "smoke_workspace"
N_SAMPLES = 2500
N_GENES = 30
SNPS_PER_GENE_TRUE = 40
CHROMS = ["1", "2", "3"]
GENE_SPACING_BP = 500_000
TRAITS = ["TRAIT_A", "TRAIT_B"]

# dosage of A1 -> PLINK 1 two-bit code
_CODE = {2: 0b00, -1: 0b01, 1: 0b10, 0: 0b11}


def write_bed(path: Path, genotype: np.ndarray) -> None:
    """genotype is [n_samples, n_variants] of {0,1,2,-1}; writes SNP-major."""
    n_samples, n_variants = genotype.shape
    stride = bed_bytes_per_variant(n_samples)
    with path.open("wb") as handle:
        handle.write(b"\x6c\x1b\x01")
        for v in range(n_variants):
            codes = np.asarray([_CODE[int(x)] for x in genotype[:, v]], dtype=np.uint8)
            padded = np.zeros(stride * 4, dtype=np.uint8)
            padded[:n_samples] = codes
            bits = np.zeros((stride * 4, 2), dtype=np.uint8)
            bits[:, 0] = padded & 1
            bits[:, 1] = (padded >> 1) & 1
            handle.write(np.packbits(bits.reshape(-1), bitorder="little").tobytes())


def build() -> dict:
    rng = np.random.default_rng(20260917)
    OUT.mkdir(parents=True, exist_ok=True)

    genes = []
    chrom_list, position_list = [], []
    for c in CHROMS:
        for k in range(N_GENES // len(CHROMS)):
            tss = 1_000_000 + k * GENE_SPACING_BP
            genes.append({"gene_id": f"G{c}_{k:03d}", "chrom": c, "tss": tss})
            offsets = np.sort(rng.integers(-150_000, 150_000, size=SNPS_PER_GENE_TRUE))
            for o in offsets:
                chrom_list.append(c)
                position_list.append(int(tss + o))
    order = sorted(range(len(position_list)), key=lambda i: (CHROMS.index(chrom_list[i]), position_list[i]))
    chrom_list = [chrom_list[i] for i in order]
    position_list = [position_list[i] for i in order]
    n_variants = len(position_list)

    frequency = rng.uniform(0.05, 0.5, size=n_variants)
    genotype = rng.binomial(2, frequency, size=(N_SAMPLES, n_variants)).astype(np.int8)
    genotype[rng.random(genotype.shape) < 0.001] = -1

    sample_ids = [f"S{i:06d}" for i in range(N_SAMPLES)]
    with (OUT / "geno.fam").open("w") as h:
        for s in sample_ids:
            h.write(f"{s} {s} 0 0 0 -9\n")
    with (OUT / "geno.bim").open("w") as h:
        for i in range(n_variants):
            h.write(f"{chrom_list[i]}\tv{i:07d}\t0\t{position_list[i]}\tA\tG\n")
    write_bed(OUT / "geno.bed", genotype)

    dosage = np.where(genotype < 0, 2 * frequency, genotype).astype(np.float64)
    standardized = (dosage - dosage.mean(axis=0)) / np.maximum(dosage.std(axis=0), 1e-8)
    chrom_array, position_array = np.asarray(chrom_list), np.asarray(position_list)

    age = rng.normal(55, 8, size=N_SAMPLES)
    sex = rng.integers(0, 2, size=N_SAMPLES).astype(float)
    pcs = rng.normal(size=(N_SAMPLES, 4))
    covariates = np.column_stack([age, sex, pcs])

    # A third of genes carry a real cis effect; burden beta tracks it.
    tau = np.where(rng.random(len(genes)) < 0.34, rng.normal(0, 0.45, size=len(genes)), 0.0)
    signal = np.zeros(N_SAMPLES)
    burden_rows = []
    for gi, gene in enumerate(genes):
        window = np.flatnonzero(
            (chrom_array == gene["chrom"]) & (np.abs(position_array - gene["tss"]) <= 200_000)
        )
        causal = rng.choice(window, size=4, replace=False)
        signal += tau[gi] * standardized[:, causal].mean(axis=1)
        for trait in TRAITS:
            scale = 1.0 if trait == "TRAIT_A" else 0.35
            burden_rows.append(
                (gene["gene_id"], trait, tau[gi] * scale + rng.normal(0, 0.06), abs(rng.normal(0.09, 0.02)))
            )

    with (OUT / "phenotypes.tsv").open("w") as h:
        h.write("sample_id\t" + "\t".join(TRAITS) + "\n")
        y_a = 1.6 * signal + 0.02 * age + 0.3 * sex + rng.normal(0, 1.0, size=N_SAMPLES)
        y_b = 0.6 * signal + 0.01 * age + rng.normal(0, 1.0, size=N_SAMPLES)
        for i, s in enumerate(sample_ids):
            h.write(f"{s}\t{float(y_a[i])!r}\t{float(y_b[i])!r}\n")
    with (OUT / "covariates.tsv").open("w") as h:
        h.write("sample_id\tage\tsex\tPC1\tPC2\tPC3\tPC4\n")
        for i, s in enumerate(sample_ids):
            h.write(s + "\t" + "\t".join(repr(float(v)) for v in covariates[i]) + "\n")
    with (OUT / "samples.tsv").open("w") as h:
        h.write("sample_id\tgroup\n")
        for i, s in enumerate(sample_ids):
            h.write(f"{s}\tFAM{i // 2:06d}\n")  # deliberate sibling pairs
    with (OUT / "genes.tsv").open("w") as h:
        h.write("gene_id\tchrom\ttss\n")
        for gene in genes:
            h.write(f"{gene['gene_id']}\t{gene['chrom']}\t{gene['tss']}\n")
    with (OUT / "burden.tsv").open("w") as h:
        h.write("gene_id\ttrait\tbeta\tse\n")
        for gid, trait, beta, se in burden_rows:
            h.write(f"{gid}\t{trait}\t{float(beta)!r}\t{float(se)!r}\n")

    binding = {
        "genotype": {"bed": str(OUT / "geno.bed"), "bim": str(OUT / "geno.bim"),
                     "fam": str(OUT / "geno.fam"), "dosage_allele": "A1"},
        "samples_path": str(OUT / "samples.tsv"),
        "covariates_path": str(OUT / "covariates.tsv"),
        "phenotypes_path": str(OUT / "phenotypes.tsv"),
        "genes_path": str(OUT / "genes.tsv"),
        "burden_path": str(OUT / "burden.tsv"),
        "traits": TRAITS,
        "chromosomes": [],
        "expected_sha256": {},
    }
    (OUT / "BINDING.json").write_text(json.dumps(binding, indent=2) + "\n")

    config = json.loads((HERE / "CONFIG.json").read_text(encoding="utf-8-sig"))
    config.update({
        "cis_radius_bp": 200_000, "snps_per_gene": 32, "maf_min": 0.01,
        "max_genes": N_GENES, "min_genes_required": 5,
        "bootstrap_replicates": 300, "permutation_replicates": 300,
    })
    (OUT / "CONFIG.smoke.json").write_text(json.dumps(config, indent=2) + "\n")
    print(f"synthetic cohort: {N_SAMPLES} samples x {n_variants} variants, {len(genes)} genes", flush=True)
    return binding


def main() -> int:
    build()
    return subprocess.call([
        sys.executable, str(HERE / "run_bar.py"),
        "--binding", str(OUT / "BINDING.json"),
        "--config", str(OUT / "CONFIG.smoke.json"),
        "--out", str(OUT / "results"),
    ])


if __name__ == "__main__":
    raise SystemExit(main())
