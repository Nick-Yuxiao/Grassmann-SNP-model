# Gate 0A sigma-pilot protocol

_Prospectively frozen. Estimates the replicate-level SD of the paired `Δ R²_genetic` on the exact bound frozen contract, then maps it to the Gate 0A replicate count `R`. Variance calibration only. Date: 2026-09-04._

---

## 📋 Question

For the Gate 0A decision procedure, what is the replicate-level SD `σ` of the paired difference `Δ_r = R²_genetic(A_test) − R²_genetic(B_pca_z)` on the frozen contract, and therefore how many replicates `R` are needed to reach ≥0.90 power at the fixed margin 0.005?

## 🔒 Order (never reversed)

Bind panel hash → bind LD-block version hash → freeze/verify preprocessing (`BINDING.json → BOUND`) → export panel → run pilot → estimate `σ` → map to `R` → freeze `R_formal`. The run script refuses unless the contract is `BOUND`, so a `σ` measured against an unfrozen SNP set / block version / split can never feed `R`.

## 🎯 Scope: variance only

The pilot estimates `SD(Δ R²_genetic)` per cell. It never declares a winner, and its results may set **only** the replicate count and compute feasibility — not the margin, not whether a DGP regime is kept, not `k`, not the arms, not the success threshold, not which regime is primary. Those are frozen upstream (rc2).

## 🧮 Why no h² grid

The primary metric `R²_genetic` fits the downstream to the **noise-free** genetic value `g` and evaluates `R²(g, ĝ)` on held-out folds. Since `g` carries no phenotype noise, the replicate SD of `Δ R²_genetic` is **heritability-independent**; the pilot therefore needs no h² grid. (h² enters only the secondary `R²_pheno` that Gate 0A reports, not the primary sizing here.)

## 🧪 Design

- Regimes (representative): `spectral_tail_adversarial` (must — the primary gate), `major_LD_aligned` (easy reference), `between_block_interaction` (high-variance reference, bilinear head).
- Budget: per-block `k ∈ {8, 16}`; `R_pilot = 30`; outer 5-fold CV; ridge λ grid `{0.01,0.1,1,10}`; 3 causal directions.
- Arms fit on the training fold only; MAF-z standardization with training-fold stats; primary arm KPCA-RBF (median-heuristic bandwidth on the training fold), null `B_pca_z`.
- Additive regimes read a bounded set of blocks (≤40, recorded) so the pilot is affordable; interaction regimes use the matched bilinear cross-product head on the involved block pair.

## 📐 σ → R and R_formal

Each cell's `σ̂` is mapped to the smallest `R` in `{20,50,100,200,400,800}` reaching ≥0.90 power at the fixed 0.005 margin, using the exact detectability decision procedure (percentile-bootstrap CI lower bound > 0). Then

```
R_formal = max over primary-relevant cells of required R    (never the average)
```

so `spectral_tail` cannot be under-powered because an easy regime needed fewer replicates. If no candidate `R` suffices, increase `R` — never lower the margin.

## 🧬 Reference implementation

`src/gate0a` implements the arms (`pca_z`, `kpca_rbf`), the per-regime DGP (`g` in the block spectral tail / leading directions / cross-block product), the matched bilinear head, and `R²_genetic` via nested CV. It is the reference Gate 0A will reuse. It is numpy-only and unit-tested on synthetic panels (`tests/test_gate0a.py`); the real `σ` comes from running it on the bound panel on the server.

## 🔒 Firewall

`VARIANCE_CALIBRATION_ONLY`. Runs only when bound; refuses a length-mismatched panel; outputs set `R` and expose feasibility, nothing else.
