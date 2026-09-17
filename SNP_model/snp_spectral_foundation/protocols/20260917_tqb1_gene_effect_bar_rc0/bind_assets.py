"""Read-only inventory and binding check for TQ-B1.

Run this before run_bar.py. It parses every input, checks that the dimensions
and identifiers agree, hashes each file, and writes the hashes back into the
binding. It never fits a model and never reports a phenotype summary, so it
cannot leak outcome information into a design decision.

    python bind_assets.py --binding BINDING.json
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

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


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binding", default="BINDING.json")
    parser.add_argument("--config", default="CONFIG.json")
    parser.add_argument("--write-hashes", action="store_true",
                        help="Store the observed sha256 values back into the binding file.")
    args = parser.parse_args()

    binding_path = Path(args.binding)
    binding = load_json(binding_path)
    config = load_json(HERE / args.config)
    problems: list[str] = []
    report: dict[str, object] = {}

    def need(key: str, value) -> Path | None:
        if not value:
            problems.append(f"{key} is not set in {binding_path.name}")
            return None
        path = Path(value)
        if not path.is_file():
            problems.append(f"{key} does not exist: {path}")
            return None
        return path

    geno = binding.get("genotype") or {}
    bed = need("genotype.bed", geno.get("bed"))
    bim = need("genotype.bim", geno.get("bim"))
    fam = need("genotype.fam", geno.get("fam"))
    samples_path = need("samples_path", binding.get("samples_path"))
    cov_path = need("covariates_path", binding.get("covariates_path"))
    phe_path = need("phenotypes_path", binding.get("phenotypes_path"))
    genes_path = need("genes_path", binding.get("genes_path"))
    burden_path = Path(binding["burden_path"]) if binding.get("burden_path") else None
    if burden_path is not None and not burden_path.is_file():
        problems.append(f"burden_path does not exist: {burden_path}")
        burden_path = None
    if problems:
        for p in problems:
            print(f"FAIL: {p}")
        return 1

    fam_ids = read_fam(fam)
    chroms, positions, variant_ids = read_bim(bim)
    check_bed_size(bed, len(fam_ids), len(positions))
    report["genotype"] = {"samples": len(fam_ids), "variants": len(positions),
                          "chromosomes": sorted(set(chroms))}
    print(f"genotype: {len(fam_ids)} samples, {len(positions)} variants, "
          f"{len(set(chroms))} chromosomes", flush=True)

    ids, groups, supplied = read_samples(samples_path)
    absent = [s for s in ids if s not in set(fam_ids)]
    if absent:
        problems.append(f"{len(absent)} analysed samples are missing from the .fam, e.g. {absent[:3]}")
    labels = assign_splits(ids, groups, supplied, config["split_fractions"], int(config["seed"]))
    counts = {k: int(np.sum(labels == k)) for k in ("train", "validation", "evaluation")}
    group_of = np.asarray(groups)
    for a, b in (("train", "validation"), ("train", "evaluation"), ("validation", "evaluation")):
        if set(group_of[labels == a]) & set(group_of[labels == b]):
            problems.append(f"group leakage between {a} and {b}")
    report["samples"] = {"analysed": len(ids), "groups": len(set(groups)),
                         "splits": counts, "split_source": "supplied" if any(supplied) else "derived"}
    print(f"samples: {len(ids)} in {len(set(groups))} groups -> {counts}", flush=True)

    cov_ids, cov_names, _ = numeric_table(cov_path)
    phe_ids, trait_names, _ = numeric_table(phe_path)
    if cov_ids != ids:
        problems.append("covariates file is not in the same sample order as the samples file")
    if phe_ids != ids:
        problems.append("phenotypes file is not in the same sample order as the samples file")
    traits = list(binding.get("traits") or trait_names)
    for trait in traits:
        if trait not in trait_names:
            problems.append(f"trait '{trait}' is not a column of the phenotype file")
    report["covariates"] = {"count": len(cov_names), "names": cov_names}
    report["traits"] = traits
    print(f"covariates: {len(cov_names)} -> {cov_names}", flush=True)
    print(f"traits: {traits}", flush=True)

    genes = read_genes(genes_path)
    wanted = set(binding.get("chromosomes") or [])
    if wanted:
        genes = [g for g in genes if g["chrom"] in wanted]
    unknown = sorted({g["chrom"] for g in genes} - set(chroms))
    if unknown:
        problems.append(f"genes name chromosomes absent from the .bim: {unknown[:5]} "
                        f"(the .bim uses {sorted(set(chroms))[:5]})")
    chrom_array = np.asarray(chroms)
    radius, per_gene = int(config["cis_radius_bp"]), int(config["snps_per_gene"])
    eligible = 0
    for gene in genes:
        pool = np.flatnonzero(chrom_array == gene["chrom"])
        if len(pool) and int(np.sum(np.abs(positions[pool] - int(gene["tss"])) <= radius)) >= per_gene:
            eligible += 1
    report["genes"] = {"supplied": len(genes), "eligible": eligible,
                       "min_required": int(config["min_genes_required"])}
    print(f"genes: {len(genes)} supplied, {eligible} with >= {per_gene} cis variants", flush=True)
    if eligible < int(config["min_genes_required"]):
        problems.append(f"only {eligible} eligible genes, {config['min_genes_required']} required")

    if burden_path is not None:
        burden = read_burden(burden_path)
        gene_ids = {g["gene_id"] for g in genes}
        overlap = {t: sum(1 for (g, tr) in burden if tr == t and g in gene_ids) for t in traits}
        report["burden"] = {"rows": len(burden), "genes_with_burden_per_trait": overlap}
        print(f"burden: {len(burden)} rows, overlap per trait {overlap}", flush=True)
        for trait, n in overlap.items():
            if n < int(config["min_genes_required"]):
                problems.append(f"trait '{trait}' has burden for only {n} analysed genes")
    else:
        report["burden"] = None
        print("burden: not bound (the run will stop after the detectability gate)", flush=True)

    hashes = {name: sha256_file(p) for name, p in [
        ("bed", bed), ("bim", bim), ("fam", fam), ("samples", samples_path),
        ("covariates", cov_path), ("phenotypes", phe_path), ("genes", genes_path),
    ] + ([("burden", burden_path)] if burden_path else [])}
    report["sha256"] = hashes
    expected = binding.get("expected_sha256") or {}
    for name, value in expected.items():
        if name in hashes and hashes[name] != value:
            problems.append(f"{name} hash mismatch: expected {value}, observed {hashes[name]}")

    dump_json(Path("INVENTORY.json"), {"status": "FAIL" if problems else "PASS",
                                       "problems": problems, "report": report})
    if args.write_hashes and not problems:
        binding["expected_sha256"] = hashes
        dump_json(binding_path, binding)
        print(f"wrote {len(hashes)} hashes into {binding_path.name}", flush=True)

    for p in problems:
        print(f"FAIL: {p}")
    print(f"\n{'FAIL' if problems else 'PASS'}: inventory written to INVENTORY.json")
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
