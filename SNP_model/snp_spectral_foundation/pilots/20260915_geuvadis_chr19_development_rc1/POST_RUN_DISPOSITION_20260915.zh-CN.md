# RC2 post-run disposition

_0526 Stage 2 bridge · 2026-09-15 · append-only decision record written after the frozen run_

---

## 🎯 Frozen disposition

> **0526 Stage 2 bridge, RC2: NO-GO in this setting.**

This wording is final for Pilot `20260915_geuvadis_chr19_development_rc1`. It must not be changed by rerunning the 15 test individuals, selecting another transcript, changing the seed, tuning the rank, replacing the decoder, or adding a new pass threshold.

| Frozen item | Value |
| --- | ---: |
| Test individuals | 15 |
| Dosage `B` R² | 0.022519 |
| Full-H `C` R² | 0.015716 |
| Primary `C-B` | -0.006803 |
| Paired 95% CI | [-0.059185, 0.044373] |
| Decision | `NO-GO in this setting` |

The immutable result is [`results/RESULTS.json`](results/RESULTS.json), SHA-256 `10e0624fba7820e93bb03c0740ad5e6eb0a2695d82b6462149845b070df879ff`. Its original 13-entry [`results/RUN_MANIFEST.sha256`](results/RUN_MANIFEST.sha256) remains unchanged; this post-run interpretation file is intentionally not represented as a pre-run artifact.

## 📋 What this result establishes

- Genotype-only pretraining learned local genotype structure in this run
- Full aligned `H` did not improve held-out predictive R² over same-panel dosage
- Grassmann arm `E` did not outperform dosage or aligned low-rank `D`
- No phenotype-conditioned weighting or biology-prior stage is authorized from RC2

## ⚠️ What this result does not establish

The run does not distinguish between an absent generic phenotype bridge and inadequate bridge sensitivity caused by a weakly predictable task plus a 15-person test set. It therefore does not, by itself, terminate the complete 0526 research direction.

The positive point estimate of `D-B` is registered only as a secondary hypothesis for an independent Stage 2B. It does not authorize any further analysis of the RC2 test individuals.

## ✍️ Next authorized question

Only a separately versioned Stage 2B may ask whether frozen pretrained `H` improves phenotype prediction after the phenotype task has demonstrated genotype predictability on data disjoint from the final bridge test.

No new foundation architecture, CropNet/CropARNet-style weighting, rank search, or Grassmann optimization is part of that question.
