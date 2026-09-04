# SUPERSEDED

_This package (rc1) is superseded by `20260904_grassmann_phenotype_signal_compression_decision_rc2` and must NOT be authorized to run._

---

## Why

Review found rc1's primary test arm (RBF kernel-PCA on raw genotype) is neither the pretrained encoder nor Grassmann, yet rc1 wrote Stage-1/2 and Grassmann-specificity kills onto it — an over-claim. rc2 corrects the scope to **Gate 0A (nonlinear compression headroom)**, deletes the Grassmann kill, and adds the frozen DGP regimes, MAF standardization, total-budget accounting, the cost axis, `R²_genetic`, replicate-based inference, and a detectability gate.

rc1 is retained only as design provenance. Use rc2.
