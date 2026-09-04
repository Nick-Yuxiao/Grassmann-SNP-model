# Gate 0A protocol — nonlinear compression headroom

_Prospectively frozen, closed-form, no-Transformer decision experiment. Design only; run execution is NOT authorized. Supersedes rc1. Date: 2026-09-04._

---

## 📋 The one proposition Gate 0A decides

> On the shared genotype panel, per pre-registered DGP regime, is there any cheap nonlinear low-dimensional map that preserves more **genetic** signal than the matched linear reconstruction null (block PCA), at fixed representation budget?

Gate 0A is the only truly no-Transformer gate. It **does not** test the pretrained encoder and **does not** test Grassmann geometry — those are Gate 0B and Gate 1 (see `GATE_LADDER_ROADMAP.md`). rc1's error was writing "Stage 1 + Stage 2" kill criteria onto a kernel-PCA experiment; that scope is corrected here.

## 🎯 Why this gate is worth running (the mechanistic risk)

Unsupervised reconstruction optimizes **genotype variance**. Trait-causal directions have no reason to align with top genotype variance — leading genotype PCs are dominated by ancestry and LD block structure, not by what is causal for a given trait. So there is a real prior that unsupervised spectral compression discards phenotype signal unless its nonlinearity specifically reorganizes it. Gate 0A tests exactly that, cheaply, and a clean negative in the adversarial regime would kill the unsupervised-compression premise for the whole programme.

## 🎯 Why block PCA is the null (and its exact status)

The objective under test is unsupervised reconstruction. PCA is the linear **reconstruction optimum**, so it is the matched linear comparator — but it is **not** an upper bound over all linear phenotype-retention methods (a supervised PLS would be a different, and here disallowed, comparator). The primary null `B_pca_z` runs on the same MAF-standardized input the kernel encoder sees, so preprocessing cannot confound the comparison; `B_pca_raw` is a preprocessing-sensitivity secondary.

## 🧪 Arms (matched total budget, identical inputs/splits/traits)

| Arm | Representation | Role |
| --- | --- | --- |
| `A_test` (primary) | RBF kernel-PCA top-`k` | primary test arm |
| `A_test` (secondary) | shallow per-block autoencoder | robustness only |
| `B_pca_z` | block PCA top-`k` on MAF-standardized input | **PRIMARY LINEAR NULL** |
| `B_pca_raw` | block PCA top-`k` on raw-centered input | preprocessing-sensitivity null |
| `B_rand_bilinear` | same-dim random bilinear features | **exploratory only, no kill** |
| `B_rand_proj` | Gaussian random projection | trivial reference |
| `B_ldprune` | LD-pruned SNP subset, budget-matched | naive selection reference |
| `B_blockmean` | block mean / haplotype dosage | floor |

`B_rand_bilinear` carries no Grassmann kill here — Gate 0A has no Grassmann arm. It is the matched comparator only at Gate 1.

## 📐 Budget, not "k"

Per-block `k` yields total dimension `D_total(k) = Σ_b min(k, m_b)`, not `k`. The decision axis is **total representation budget** and compression ratio; per-block `k` is reported alongside. Primary budget points are per-block `k ∈ {8, 16}` (fixed in advance to avoid `argmax_k` selection bias); the grid `{1,2,4,8,16,32,64}` is exploratory.

## 🧬 DGP regimes (decided separately)

Six pre-registered regimes (`config/DGP_REGIMES.json`), quantitative traits only, `h² ∈ {0.05,0.1,0.2,0.4}`, each with a polygenic background:

1. `major_LD_aligned` — additive, causal on leading-PC/strong-LD SNPs (PCA-favorable).
2. `spectral_tail` — additive, causal in the low-variance tail (**adversarial; the primary decision regime**).
3. `low_MAF` — additive, causal on low-MAF variants.
4. `within_block_interaction` — `γ·G_i·G_j` inside a block; matched interaction head.
5. `between_block_additive` — `β_A G_A + β_B G_B` across distant blocks; **distributed polygenic, NOT long-range interaction**.
6. `between_block_interaction` — `γ·G_A·G_B` across distant blocks; **the true long-range test**; matched interaction head.

For interaction regimes the downstream is ridge on the degree-2 polynomial expansion of the representation, applied **identically to every arm** — a linear ridge cannot read an interaction out of any representation, so an interaction-capable head is mandatory, and being identical across arms it cannot advantage one. A representation that discarded the interacting factors fails here, which is the intended signal.

## 📊 Metrics

- **Primary decision metric:** `R²_genetic(k) = R²(g, ĝ)` where `ĝ = downstream(Z_k)` fit to the **known** genetic value `g`. This separates representation loss from phenotype noise and avoids the `p ≫ n` full-genotype-ceiling problem.
- **Secondary realistic metric:** `R²_pheno(k) = R²(y, ŷ)`.
- Binary/AUC is dropped from round 1 (no penetrance/liability/prevalence/scale defined).

## 🔬 Inference (replicate-based, not fold-based)

The inferential unit is the **simulation replicate (seed)**, not the CV fold. Per replicate `r`, `Δ_r = R²_genetic(A_test, r) − R²_genetic(B_pca_z, r)`; the 95% CI is a percentile bootstrap over replicates (10,000 resamples). CV folds are used only for penalty selection and fitting.

## 🔪 Kill criterion (Gate 0A only)

Within **each** pre-registered DGP regime **separately**, at the primary budget points, if `A_test`'s paired-mean `Δ R²_genetic` vs `B_pca_z` has a replicate-bootstrap 95% CI that does **not** exclude zero, then unsupervised nonlinear compression has no headroom **in that regime**.

> The paired-mean criterion is evaluated separately within each pre-registered DGP regime; no pooled cross-regime average is used for the primary decision.

Program implication: if `A_test` fails even the `spectral_tail` and `between_block_interaction` regimes, the unsupervised spectral-compression premise is wrong and the unsupervised path to Stage 1/2 stops; regimes where it passes define the headroom carried to Gate 0B. There is **no** Grassmann kill and **no** blanket "Stage 1+2 don't exist" kill in this gate.

## 💰 Cost axis (restored)

Record bytes/sample, peak memory, encode time, downstream fit/predict time, and complexity vs #SNPs. Equal retention at lower budget/cost is a win; a tiny retention gain at large cost is not. This restores the two-dimensional (retention × cost) judgement that was a frozen conclusion.

## 🚧 Pre-run gates (all required)

1. Fix and verify **centering AND standardization** using training-fold stats only — an uncentered or unscaled PCA null is not a valid comparator.
2. Bind the data-contract hashes (`DATA_CONTRACT.json`).
3. Run the **detectability gate**: simulate known `Δ ∈ {0, 0.002, 0.005, 0.01}` and fix the margin and the replicate count so the design can actually separate them (20 seeds is likely underpowered at `h²=0.05, N≈2247`).
4. Run a **compute smoke** (3 blocks × 2 `k` × 1 outer fold) to validate the one-machine/one-week claim, or pre-register the randomized-SVD / Nyström / Lanczos eigensolver fallback rather than switching algorithms mid-run.
5. Obtain explicit project-lead run authorization; then write the run code.

## 🔒 Evidence firewall

`DESIGN_ONLY_RUN_NOT_AUTHORIZED`: no execution, no Transformer, no Gate 0B/1/2, no phenotype gating, no GPU/A1-R, no binary/real-phenotype claim, no architecture GO/NO-GO, no HGDP. When later authorized, results are per-regime headroom verdicts feeding the gate ladder — nothing more.
