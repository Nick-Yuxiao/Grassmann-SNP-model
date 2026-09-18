# S0: ridge-specification diagnostic (chr18, development only)

Tests one thing: **is `macro R2(B) < macro R2(A)` an artefact of penalising the
covariate block?** Arm B's design matrix contains arm A's, so B should not be
able to do worse — unless the fit prevents it.

This is **not** a decision pilot. It reads `dev_train` and `dev_validation`
only, it cannot consume a decision resource, and it cannot overturn a frozen
verdict. `task_gate` and `bridge_test` are refused by `resolve_development_rows`
and the refusal is tested.

## Why it exists

`TQ-G1 rc0` returned `NOT-SUPPORTED` with `C_gene - B = -0.031` (CI excluding
zero), but its own secondary contrast was `B - A = -0.014`: adding 256 cis
dosages did not beat sex + population + 10 genotype PCs. Stratifying TQ-G1's
200 genes by `target_T_g` (measured on `dev_validation`, so not circular) shows
`B <= A` in every stratum, including the 14 genes with the highest cis
predictability. The same pattern is on record in `Stage 2B` (`B-A = +0.011`,
`TASK-INELIGIBLE`) and `TQ1` (`B-A = -0.017` on 230 individuals).

All four pilots fit arms the same way — `run_stage2b.fit_ridge_predictions` and
`tqg1_core.fit_ridge` both standardize every column of
`[covariates | genotype features]` and solve one dual ridge with a single alpha
over the whole design matrix. After standardization each column contributes
equally to the Gram matrix:

| arm | features | covariate share of the kernel |
| --- | ---: | ---: |
| A | 12 | 100% |
| B | 12 + 256 | 4.5% |
| C_full | 12 + 4096 | 0.3% |

`(K + alpha I)^-1` shrinks every direction together. A small alpha leaves 256
dosage directions essentially unpenalized on 120 training rows; a large alpha
crushes the covariate block along with them. No alpha recovers arm A's fit from
arm B's design matrix. That is a hypothesis with a clean fingerprint, and this
package is what tests it.

## The three fits

Nothing varies between them except how the covariate block is treated.

| | fit |
| --- | --- |
| `S1` | every column standardized, one alpha over `[covariates \| dosage]` — the programme's current specification, reproduced |
| `S2` | covariates unpenalized, dosage penalized, fitted jointly (partial ridge via Frisch-Waugh-Lovell) |
| `S3` | covariates fitted by OLS, dosage ridged on the residual (standard two-stage eQTL recipe) |

Arm A is the covariate-only model under the same treatment: penalized ridge in
`S1`, OLS in `S2` and `S3`. `S2` and `S3` are identical for arm A, which the
validator checks for free.

`test_s0.TestS1Fidelity` asserts `S1` reproduces `tqg1_core.fit_ridge`
prediction-for-prediction, so if `S1` and `S2`/`S3` diverge on real data, the
fit is the only thing that changed.

## Reading the result

| outcome | meaning |
| --- | --- |
| `B-A > 0` under `S2` and `S3` with CI excluding zero, while `S1` stays at or below zero | the artefact is real. Every `B-A` in the programme was measured through a fit that could not let B win, and the first rung of the hierarchy has never been tested fairly. |
| `B-A` clears zero under no specification | the artefact is refuted. Suspicion moves to the uniform thinning (256 SNPs across 2 Mb is one per ~7.8 kb, so lead cis-eQTLs are probably absent from the panel), to 120 training rows, or to the upstream expression normalisation. |

The `snp_counts` axis is the sharper read. Under `S1` the damage should grow
with the size of the dosage block; under `S2`/`S3` it should not.

## What you have to supply

One line. Set `asset_dir` in `PRE_RUN_BINDING.json` to the directory TQ-G1 used.
The six chr18 GEUVADIS sha256 values are pinned from the `rc1` run and enforced,
so a wrong or re-downloaded file aborts the run. The sample split is reused from
the frozen `Stage 2B rc1` `SAMPLE_MANIFEST.tsv` and never regenerated.

## Order of operations

```bash
python -m unittest test_s0 -v            # 19 asset-free contract tests
python smoke_s0.py                       # synthetic end-to-end, ~7 s
python validate_s0.py --results results_smoke
# only once the three above are clean, and asset_dir is bound:
python run_s0.py
python validate_s0.py
```

No encoder is trained and no GPU is touched — arms A and B are the only arms, so
this is minutes of CPU, not hours.

## Known limits

- Cross-validated within development, so this is not a held-out replication.
- AF, QC and the genotype PCs are fitted on `dev_train` exactly as TQ-G1 fits
  them, so that the panel is identical. Those 120 rows also appear in CV test
  folds, which leaks a little — equally into every arm, including A, so the
  `B - A` contrast is common-mode and the leak cannot manufacture the effect
  under test.
- Uniform thinning is kept on purpose. A null here does not rule out that the
  lead cis-eQTLs are simply not in the panel; that is a separate diagnostic.
- No encoder arm is evaluated. This says nothing about whether the
  representation has value.
