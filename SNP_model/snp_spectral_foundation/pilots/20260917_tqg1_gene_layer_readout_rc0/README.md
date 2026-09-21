# TQ-G1: annotation-free gene-layer readout Pilot (chr18, rc0)

Tests one premise before any UKB or WES work is planned: **does pooling a frozen
shared genotype encoder into a gene unit beat the same-window additive dosage
ridge, when the gene is the unit of replication?**

It is **not** a burden test. No WES, no LoF, no rare variants, no functional
annotation. SNP-to-gene routing is a deterministic function of TSS distance
only, so a future comparison against real burden statistics stays uncontaminated.

Read `FROZEN_PROTOCOL.zh-CN.md` before running anything.

## What you have to supply

One line. Set `asset_dir` in `PRE_RUN_BINDING.json` to the directory that already
holds the six chr18 GEUVADIS assets used by
`../20260915_geuvadis_chr18_stage2b_rc1`:

```
GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.chr18.{pgen,pvar,psam}
GEUVADIS.445_samples.expression.bed.gz
integrated_call_samples_v3.20130502.ALL.panel
integrated_call_samples_v3.20250704.ALL.ped
```

Their SHA-256 values are already pinned from the rc1 run and are enforced, so a
wrong or re-downloaded file aborts the run instead of silently changing results.
The sample split is reused from rc1's frozen `SAMPLE_MANIFEST.tsv` — it is never
regenerated here.

Nothing else needs binding. No annotation file, no gene model, no eQTL table.

## Order of operations

```bash
python -m unittest test_tqg1 -v     # 26 asset-free contract tests
python smoke_tqg1.py                # synthetic end-to-end run, ~10 s, writes results_smoke/
python validate_tqg1.py --results results_smoke
# only once the three above are clean, and asset_dir is bound:
python run_tqg1.py
python validate_tqg1.py
```

`run_tqg1.py` refuses to run twice once `results/FINAL_STATUS.json` exists.

## Held-out policy

| Role | n | Use here |
| --- | --- | --- |
| `dev_train` | 120 | QC, AF, PCs, encoder pretraining, ridge fitting |
| `dev_validation` | 30 | ridge alpha, encoder validation, the rehearsal target `T_g` |
| `task_gate` | 50 | **evaluation set for this Pilot** |
| `bridge_test` | 230 | **sealed** |

`bridge_test` stays sealed. Pointing `evaluation_role` at it aborts unless you
pass `--open-sealed-test` *and* write an `OPEN_AUTHORIZATION.json`. Both guards
are tested.

`task_gate` was already opened once by Stage 2B rc1 for its `B-A` gate, so it is
not a virgin test set. The estimand here is different, but that limitation is
recorded in `CONFIG_FROZEN.json` and echoed into every result file.

## Arms

| Arm | Input |
| --- | --- |
| A | sex + population + 10 `dev_train`-fitted genotype PCs |
| B | A + 256 additive cis dosages — the baseline that matters |
| C_gene | A + annotation-free gene-pooled readout of the frozen encoder |
| C_full | A + 256 x 16 coordinate-aligned hidden state (the rc1 arm, for reference) |

Primary: macro `R2(C_gene) - R2(B)` over genes, two-way clustered bootstrap over
individuals **and** genes.

## Two deliberate differences from Stage 2B rc1

1. **One shared encoder** with `relative_continuous` positions and no block
   lookup table, instead of one encoder per trait. A per-trait encoder cannot
   support a "shared local encoder" claim.
2. **Cis SNPs are chosen by uniform genomic thinning**, not by top development
   `|r|`. Phenotype-driven SNP selection inflates B and makes `C-B` hard to read.

Genes are capped by a fixed hash of the gene ID, never by signal.

## The burden-test rehearsal

Per gene, `s_g` is the spread across evaluation individuals of
`prediction - prediction with the gene's cis genotypes removed`, computed for the
model arm and for the linear arm. Both are correlated against `T_g`, the
`dev_validation` cis predictability of that gene.

The number that decides anything is the difference
`spearman(s_model, T) - spearman(s_linear, T)`, with a matched-null permutation
that holds cis SNP count, MAF spectrum and expression variance fixed.

If the model score does not rank genes better than the linear score does, the
foundation encoder adds nothing at the gene level and the UKB + WES plan should
not start. Migrating to a real burden test later changes exactly one thing:
`T_g` becomes the independent WES `beta_g`. Estimator, null and decision rule
stay as written.

## Files

| File | Role |
| --- | --- |
| `FROZEN_PROTOCOL.zh-CN.md` | pre-registration; freeze before running |
| `CONFIG_FROZEN.json` | every hyperparameter and decision rule |
| `PRE_RUN_BINDING.json` | the one file you edit |
| `tqg1_core.py` | asset-free statistics, pooling and selection |
| `run_tqg1.py` | asset binding and the single real run |
| `smoke_tqg1.py` | synthetic end-to-end rehearsal |
| `validate_tqg1.py` | independent re-derivation of every headline number |
| `test_tqg1.py` | contract tests, no assets required |
