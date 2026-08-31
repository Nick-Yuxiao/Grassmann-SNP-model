# Grassmann V7 protocol addendum v7.1.3

Status: **FINAL CHR22 SITE RULE FROZEN BEFORE MODEL OUTCOMES**  
Amends: v7.1.2 HGDP-evaluability placeholder and final sequence length  
Preserves: all v7.1.0-v7.1.2 sample, model, masking, metric, GPU and source decisions

## Stage-2 evidence

After donor GT MAF >0.01 and rejection of every record at 368 duplicated
positions, 154,850 exact CHROM/POS/REF/ALT keys remained. All 154,850 keys were
present in the frozen HGDP evaluation extraction. For every key, all 768 HGDP
individuals had called diploid GT (`AN=1536`, `F_MISSING=0`).

## Final site contract

A chr22 site is included exactly when it satisfies all of the following:

1. it belongs to the source BCF pinned in v7.1.2;
2. source FILTER is the unset value `.`, truthfully labelled
   `SOURCE_PREFILTERED_FILTER_UNSET` rather than PASS;
3. it is a biallelic SNP;
4. MAF recomputed from GT over all 2,496 frozen 1KGP donors is strictly >0.01;
5. its genomic position occurs exactly once among records passing steps 1-4;
6. its exact CHROM/POS/REF/ALT key is present in the HGDP extraction;
7. every frozen HGDP individual has a called diploid GT at the key
   (`AN=1536` and `F_MISSING=0`).

The resulting sequence length is **L=154,850**. This is a data-derived result,
not a padded, truncated or MAF-tuned length. It differs from the 81,920-site
SNPBag compatibility expectation and reinforces the existing NOT_COMPARABLE
status for exact SNPBag calibration.

## Holdout-frequency non-selection

HGDP alternate-allele count is not an inclusion criterion. The 104 final sites
with `HGDP AC=0` remain in the primary panel because removing them after reading
holdout allele states would fit preprocessing to the external test population.
Metrics must additionally report a sensitivity restricted to the 154,746 sites
polymorphic in HGDP, but that sensitivity cannot replace the all-site primary
result.
