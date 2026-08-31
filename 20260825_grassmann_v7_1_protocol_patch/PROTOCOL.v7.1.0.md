# Grassmann V7 protocol v7.1.0

Status: **FROZEN FOR P0 RE-SIGN; A1 MUST NOT START BEFORE T01-T04 PASS**  
Supersedes for new runs: v7.0.1  
Preserves: all v7.0.1 files and completed T00 evidence

## 1. T01 branch decision

T01 Branch B is a **CONDITIONAL GO**.

- Use unrelated 1KGP haplotypes as the HAPNEST donor base.
- Generate synthetic training and validation corpora only after donor individuals
  (and relatedness groups) have been split. Synthetic train and validation
  corpora must never share donor individuals.
- Hold every HGDP individual out of generator fitting, donor copying, model
  selection and preprocessing fit. HGDP is real external evaluation data.
- HAPNEST is the primary reference-based generator. An msprime population-history
  simulation is a generator-sensitivity experiment, not a substitute for HAPNEST.
- The first scientific round is chr22 only. It can decide engineering feasibility
  and within-chromosome masked-genotype headroom, not genome-wide biological
  generality.

## 2. Site definition and sequence length

Freeze the site-selection algorithm, not a convenient tensor length:

1. GRCh38 chromosome 22, biallelic SNPs, PASS only.
2. Alleles must be harmonized and uniquely represented.
3. MAF is computed in the frozen unrelated 1KGP donor pool and must exceed 1%.
4. A site must be evaluable in the frozen HGDP panel.
5. Sites are coordinate sorted; duplicates and unresolved strand/allele conflicts
   are rejected.

`L` is the number of sites produced by those rules. `81,920` is an external
SNPBag-compatibility expectation only when the exact releases, samples and filters
match. No truncation, padding or MAF adjustment may be used to force that number.

## 3. Masking design

- Random masking at rates 0.50, 0.90 and 0.99 is diagnostic. It demonstrates how
  apparent headroom changes as local observations become sparse.
- The primary confirmatory cells are `ld_block` and
  `within_chrom_longrange` at overall target mask rate 0.90.
- The same structured masks at 0.99 are preregistered stress sensitivities. They
  are reported but cannot rescue a failed 0.90 primary gate.
- A 0.50 random-mask tie is an expected negative-control pattern and is not by
  itself a scientific NO-GO.

## 4. Comparator architecture and calibration

The primary fairness comparison contains three approximately 8M-parameter arms:

- `local_attn_8m_w256`
- `local_attn_gpc_8m_w256`
- `grassmann_full_8m_w256`

The local attention contract is at most 256 key tokens per query per layer. The P0
profiler uses a block-centered implementation of that contract and records the
exact implementation name. It is not labelled an exact SNPBag reproduction.

`snpbag_chr22_reproduction_38m` is a separate external-integrity calibration arm.
It must match the published site list, evaluation subset and architecture closely
enough for a declared comparison. It is not included in the 8M matched-parameter
or matched-compute superiority gate. A calibration result must be one of PASS,
FAIL or NOT_COMPARABLE, with reasons.

## 5. HGDP evaluation and uncertainty

- The published 216-individual, 54-population subset is reserved for the SNPBag
  calibration arm. Substitution is not called an exact reproduction.
- All other eligible unrelated HGDP individuals form the primary scientific
  holdout. The two sets are disjoint.
- The primary estimator gives each HGDP population equal weight; individual-weighted
  results are secondary.
- Training uncertainty and evaluation-population uncertainty are separate:
  - seed CI: paired data-seed differences after averaging two initialization seeds;
  - population CI: 10,000 seeded block-bootstrap replicates resampling HGDP
    populations, after averaging the frozen training-seed estimates.
- Both lower bounds must clear the frozen practical margin in every required
  fairness regime.

For the two primary confirmatory mask families, familywise alpha is controlled at
0.05 with Bonferroni simultaneous intervals. Diagnostic and stress-sensitivity
cells are not promoted to confirmatory after outcomes are seen.

## 6. 1KGP-neighbour sensitivity

Near/far status is derived before model outcomes from donor-only LD-pruned PCs.
Fit PCs on 1KGP donor-train individuals, project HGDP, and compute each HGDP
individual's distance to the nearest 1KGP population centroid. A population is
`near` when its median distance is within the frozen 95th percentile of 1KGP
leave-one-out own-population distances; otherwise it is `far`. Retain the
continuous distance as the primary effect-modification variable. Near/far results
are prespecified sensitivity analyses and do not replace the all-HGDP primary
estimate.

## 7. Sequential gates and claim ceiling

The order is A1 -> A2 -> A3. A failed or inapplicable upstream gate prevents the
downstream confirmatory claim. Exploratory results remain reportable as such.

Branch B can support claims about comparative masked-genotype prediction,
computational efficiency and architecture headroom under the frozen public-panel
and generator distribution. It cannot establish label efficiency, phenotype
utility, discovery beyond LD/MAF/ancestry, real-biobank performance or genome-wide
generality from chr22 alone.

Any later change to donors, holdout membership, site rules, mask role/rate,
comparator semantics, primary estimator, confidence procedure or gate logic
requires a new protocol version and manifest before outcomes are inspected.
