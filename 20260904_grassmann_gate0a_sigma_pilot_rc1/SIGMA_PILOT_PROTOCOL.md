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

## 🧪 Design (contract-faithful)

- **Training population only.** The pilot runs on the 2247 donor-train rows; the 249 donor-validation rows are never touched.
- **Folds drawn once.** The outer 5-fold split is drawn a single time and shared byte-identically across all regimes, budgets, arms, and replicates.
- **A replicate is one fixed `g`.** `g_r` is drawn once per replicate (named component + diffuse polygenic background, `frac=0.5`) and predicted across all five outer folds — `g` never changes per fold.
- **Real nested CV.** Penalty selection is a nested inner 5-fold CV (not a single 80/20 split).
- **Caching.** Arm representations do not depend on `g`, so they are computed once per (arm, fold, block) at `k=16` and sliced to `k=8`; this turns the run from tens of thousands of KPCA eigendecompositions into `2 arms × 5 folds × |blocks|`.
- Regimes (representative): `spectral_tail_adversarial` (must — the primary gate), `major_LD_aligned` (easy reference), `between_block_interaction` (high-variance reference, bilinear head). Per-block `k ∈ {8, 16}`; `R_pilot = 30`; ridge λ `{0.01,0.1,1,10}`; 3 causal directions.
- Arms fit on the fold-train only; MAF-z standardization with fold-train stats; primary arm KPCA-RBF (median-heuristic bandwidth on fold-train), null `B_pca_z`. Additive regimes read a bounded block set (≤40, recorded); interaction regimes use the matched bilinear cross-product head on the involved block pair.

## 🔗 Provenance & isolation

The run refuses unless `BINDING.json` is `BOUND` **and** an `EXPORT_MANIFEST.json` (written by `write_export_manifest.py`) ties the four exported arrays to the bound panel/block SHA-256 and confirms L=154850 / train=2247 / val=249. So the pilot cannot silently run on an export that was not derived from the bound artifacts, nor on a length-mismatched panel. `BINDING.json` is excluded from `MANIFEST.sha256` (mutable state), so binding does not invalidate the frozen-design manifest.

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
