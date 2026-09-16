# rc2 append-only amendments

## Initial freeze — 2026-09-15

- Representation-specific phenotype results accessed: no.
- Primary contrast: `C−B` only.
- Removed from rc1: arbitrary `0.005` SESOI, P1 prerequisite, Holm family, Grassmann non-inferiority margin, external-anchor and a-priori-power execution gates.
- Data fallback: public Bioconductor GEUVADIS chr19 paired genotype-expression package because no authorized UKB phenotype is available and the EBI 87 MB matrix endpoint was not practically downloadable in the current session.

## Pre-run implementation binding — 2026-09-15

This clarification was written before extracting `H`, fitting any phenotype decoder, or opening test outcomes. It changes no estimand, contrast, support rule, seed, or alpha grid.

- Split is family-disjoint where 1000 Genomes pedigree metadata are available; an individual with no pedigree row forms a singleton group. Approximate 70/15/15 fractions are produced by two deterministic group splits.
- Expression is `log1p(expected transcript count)`. Candidate traits are restricted to the package-published `gene_id_subset.txt`, must be finite and nonzero in at least 80% of train individuals, and are ranked by train-only variance with transcript-ID tie breaking.
- Fixed sample covariates are sex plus a sex-missing indicator. Sex is train-mean imputed. Ten train-fitted genotype PCs are included in every arm.
- Candidate SNPs are restricted to package-published `snp_id_subset.txt`. Train-only filters are missingness <= 0.10 and MAF >= 0.05. Eligible SNPs are coordinate-sorted; at most 1,280 are selected at deterministic evenly spaced indices, then grouped into consecutive 32-SNP blocks. At least 256 SNPs are required.
- TSV SNP coordinates inherit the source VCF convention: GRCh37 chromosome 19, 1-based variant positions; `start == end` is required. No BED conversion or liftover occurs.
- Local encoder is fixed at width 16, two Transformer layers, four heads, 40 genotype-only epochs, mask probability 0.20, batch size 16, AdamW learning rate 0.002 and weight decay 0.0001.
- `D` is a rank-4 train-fitted PCA of the hidden channel axis, retaining aligned per-SNP scores. `E` is the rank-4 per-person/per-block hidden covariance subspace projector upper triangle plus normalized `log1p` eigenvalues.
- Test opening remains single-shot. Any implementation failure makes this Pilot version invalid; repair requires a new Pilot version rather than overwriting a valid result.
