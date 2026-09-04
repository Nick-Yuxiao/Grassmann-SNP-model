# Phenotype-signal compression decision protocol

_Prospectively frozen, closed-form decision experiment. Design only; run execution is NOT authorized by this package. Date: 2026-09-04._

---

## 📋 The one question this experiment answers

Before any further Stage-1 (spectral compression) or Stage-2 (Grassmann mixing) investment, decide the single upstream question the whole architecture rests on:

> For an LD block, does the top-`k` of a learned nonlinear reconstruction representation preserve more phenotype signal than the information-matched best linear baseline at the same `k`?

If the answer is no, Stage 1 + Stage 2 have no reason to exist and the GPU / A1-R programme stays stopped. This is the foundation-load test that the four-stage design skipped.

This experiment trains **no Transformer** and uses **no GPU**. It is closed-form / CPU and is intended to finish on one machine within about a week.

## 🎯 Why block PCA is the correct null (and PLS is not)

The compression objective under test is **unsupervised reconstruction**. Under a variance-maximizing (reconstruction) objective, the top-`k` subspace of a *linear* encoder **is** the block PCA subspace. Therefore a near-linear nonlinear encoder degenerates, in the exact mathematical sense, to block PCA, and block PCA is the precise linear ceiling the test arm must clear.

A phenotype-aware null (block-PLS / supervised projection) is deliberately **excluded**: it would let the test arm "win" by using label information the encoder never saw, converting the comparison from "nonlinear vs linear" into "used labels vs did not." The null must see exactly the information the encoder sees — for a reconstruction encoder that is block PCA.

## 🧪 Arms (all at matched `k`, identical inputs/splits/traits)

| Arm | Representation | Role |
| --- | --- | --- |
| `A_test` (primary) | RBF kernel-PCA top-`k` (closed-form) | the nonlinear reconstruction compressor under test |
| `A_test` (secondary) | shallow per-block autoencoder bottleneck `k` | robustness only, not the decision |
| `B_pca` | block PCA top-`k` | **PRIMARY LINEAR NULL** |
| `B_rand_bilinear` | random bilinear / same-dim quadratic features | **reinstated ablation**: geometry vs generic nonlinearity |
| `B_rand_proj` | Gaussian random projection to dim `k` | trivial lower reference |
| `B_ldprune` | LD-pruned SNP subset, count-matched to `k` | naive SNP-selection reference |
| `B_blockmean` | block mean / simple haplotype dosage | floor reference |

`B_rand_bilinear` is the ablation the latest design dropped. It is not redundant with `B_pca`: PCA answers "did nonlinearity beat the linear optimum," while random bilinear answers "did the Grassmann *geometry* do anything beyond adding a same-dimension quadratic feature." Both questions must be answered on the same curve.

## 📊 Downstream evaluation

Fixed downstream: ridge regression (continuous) / logistic regression (binary), penalty chosen by inner CV on training folds only, nested 5×5 CV. For each arm and each `k`, report R² / AUC as a curve in `k`, per trait and aggregated across the trait panel.

## 🧬 Phenotype simulator (unifies both readouts)

The frozen panel is 1KGP genotypes with no rich real phenotypes, so all phenotype-signal readouts are on **simulated** traits with a known genetic architecture: empirical LD from the frozen v7 chr22 panel, a polygenic background, known causal positions, and an h² grid `{0.05, 0.1, 0.2, 0.4}`. Causal effects are injected in two regimes — **within-block** and **between-block** — reported separately.

The R²/AUC-vs-`k` curve and the causal-retention-vs-compression curve are two readouts of the **same** simulated runs. `causal_signal_retention` is defined as the fraction of injected causal variance recoverable by the frozen downstream model from the `k`-dim representation, relative to the full-genotype ceiling. The **between-block** retention curve is the only direct evidence of LD-preserving / long-range compression — the actual Grassmann pitch.

## 🔪 Pre-registered kill criteria (written now, frozen before any run)

1. **Stage-1/2 kill.** If `A_test` does not beat `B_pca` at matched `k` across the trait panel in **both** within-block and between-block regimes, by a margin whose CV-based uncertainty interval excludes zero, then Stage 1 + Stage 2 have no reason to exist.
2. **Geometry-specificity kill.** If `A_test` beats `B_pca` but does not beat `B_rand_bilinear` by the same standard, the win is a generic quadratic feature and the Grassmann-specific claim fails even if nonlinearity helps.

The margin threshold and the CV-fold uncertainty rule are frozen before any run and may not be adjusted after seeing arm results.

## 🚧 Pre-run gates (all required before run authorization)

1. Fix and verify the previously flagged **centering defect** — every reconstruction/PCA arm depends on correct per-block centering; an uncentered PCA null is not a valid ceiling.
2. Bind the data contract to the frozen v7 chr22 panel manifest hash and the block-version hash (see `DATA_CONTRACT.json`).
3. Freeze the margin and uncertainty rule.
4. Obtain explicit project-lead run authorization; then, and only then, write the run code.

## 🔒 Evidence firewall

This package is `DESIGN_ONLY_RUN_NOT_AUTHORIZED`. It authorizes no execution, no Transformer training, no phenotype gating, no GPU or A1-R work, and no reading of the HGDP holdout. When later authorized and run, results decide only the stated upstream question; results on simulated phenotypes may not be promoted to a real-phenotype prediction claim, nor to an architecture GO/NO-GO beyond the stated kill logic. Any later formal use requires binding the pre-run gates above.
