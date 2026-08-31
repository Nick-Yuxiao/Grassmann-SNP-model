# v7.1.4 panel materialization addendum

Status: frozen for T02 materialization.

This addendum does not change the v7.1.3 sample or site selection. It defines
the lossless materialization of those frozen selections from the audited source
BCF.

## Frozen rules

- Source BCF SHA-256 must be
  `09e50c698522d19bbc2c43a2ea2d29b290d63cf8c1d6983c5198dcb4289ff3fa`.
- Variants are intersected by exact `CHROM,POS,REF,ALT` identity. Position-only
  matching is forbidden.
- The expected sequence length is 154,850 and is not forced or subsampled.
- All source `INFO` annotations and the unused malformed `FORMAT/PP` declaration
  and values are removed. This prevents source-wide AC/AF-like annotations from
  entering a donor-only artifact. Phased `FORMAT/GT` is retained unchanged; no
  genotypes are imputed or rewritten.
- Five indexed BCF artifacts are written: all source samples, all frozen release
  samples, donor train, donor validation, and HGDP primary.
- Every artifact must contain the same 154,850 variants in frozen order and the
  exact expected sample set for its role.
- The 3,264-sample joint release artifact must have zero sites with missing GT
  and zero sites containing an unphased `/` genotype.
- Every materialized header must contain no `INFO` IDs and exactly one `FORMAT`
  ID, `GT`.
- Source data and source index are read-only. The source-index mtime warning is
  recorded but does not authorize replacing the source CSI.
- Materialization is CPU/I/O-only. No GPU is selected or used.

## Transaction and non-interruption rule

The runner holds a non-blocking file lock, refuses an existing final directory,
and refuses to start while a matching source-panel `bcftools` task owned by the
same user is active. It never sends signals to another process. Outputs are
built under a unique temporary directory and renamed to the final directory
only after all audits and hashes pass.
