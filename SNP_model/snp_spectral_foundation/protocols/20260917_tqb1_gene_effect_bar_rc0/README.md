# TQ-B1: the additive gene-effect bar

Measures how well a **purely additive** cis model already predicts an independent
WES burden effect. That number is the bar any foundation model has to clear.
Knowing it costs days; discovering it after building a model costs months.

**No GPU. No PyTorch. No compiled genetics library. numpy only.**

## Why this exists

TQ-G1 rc0 spent a run comparing representations on a regression where additive
dosage did not beat covariates. Nothing could be concluded. This package puts
the detectability question first and makes it a hard gate: if adding cis dosage
does not raise held-out R2 on the validation split, the evaluation split is
never opened.

## Install

```bash
python -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt      # numpy, nothing else
```

Python 3.9+.

## Prove it works before binding anything

```bash
python -m unittest test_tqb1 -v      # 20 tests, no data needed
python smoke_bar.py                  # builds a synthetic cohort, runs the real pipeline
python validate_bar.py --results smoke_workspace/results \
                       --config smoke_workspace/CONFIG.smoke.json \
                       --burden smoke_workspace/burden.tsv
```

Expected: `20 tests OK`; the smoke gate passes `TRAIT_A` and stops `TRAIT_B`; the
validator reports `18/18`. The smoke run takes about a second and proves the
whole chain including the PLINK reader, the gate and the bar.

## Inputs

Copy `BINDING_TEMPLATE.json` to `BINDING.json` and fill in seven paths. Six
inputs are plain TSVs; only the genotype is binary.

| Key | What |
| --- | --- |
| `genotype.{bed,bim,fam}` | PLINK 1 **SNP-major** .bed. `plink2 --make-bed` produces it. |
| `samples_path` | `sample_id` + optional `group` + optional `split` |
| `covariates_path` | `sample_id` + numeric columns (age, sex, PCs, centre, batch) |
| `phenotypes_path` | `sample_id` + one column per trait |
| `genes_path` | `gene_id`, `chrom`, `tss` (1-based; `chrom` must match the .bim) |
| `burden_path` | `gene_id`, `trait`, `beta`, `se` from an **independent** WES analysis |

`group` is the kinship or family component. Splits are made by group, so
relatives never straddle a boundary. **If you leave `group` out, every sample is
treated as unrelated** and a related cohort will leak across splits.

Leave `burden_path` null to run the gate alone.

## Run

```bash
python bind_assets.py --binding BINDING.json --write-hashes   # read-only inventory
python run_bar.py     --binding BINDING.json --out results
python validate_bar.py --results results --burden <burden.tsv>
```

`bind_assets.py` parses and hashes every input and checks that identifiers and
dimensions agree. It fits nothing and prints no phenotype summary, so it cannot
leak outcome information into a design decision. Fix everything it reports
before running the analysis.

Try a subset first: `python run_bar.py --binding BINDING.json --max-genes 50`,
and `--gate-only` to stop after the gate.

## What comes out

| File | Contents |
| --- | --- |
| `GATE.json` | detectability gate per trait; which traits may proceed |
| `RESULTS.json` | gate, held-out macro R2, the bar, limits |
| `GENE_TABLE.tsv` | per gene: perturbation score and held-out R2, full precision |
| `FINAL_STATUS.json` | one-screen summary |
| `RUN_BINDING.json` | input hashes, code hashes, runtime |
| `VALIDATION.json` | independent re-derivation |

The number you want is `bar.<trait>.matched_null.observed_spearman`.

## Reading the bar

`s_g` is the spread, across evaluation individuals, of how much the prediction
moves when gene `g`'s cis genotypes are replaced by the training mean. It is a
**magnitude**, not a signed effect, so it is ranked against `|beta_g|`.

- **Bar near 0.8** — additive cis prediction already recovers the burden ranking.
  A foundation model has very little room and probably is not worth building for
  this purpose.
- **Bar near 0.3** — additive methods resolve genes poorly. There is real room,
  and a model that beats it would be saying something.
- **Matched-null p not small** — the ranking is explained by cis SNP count, MAF
  spectrum or burden precision alone. Then the bar is not measuring biology and
  the strata need rethinking before anything is built on it.

## Scale

Ridge is solved in the primal, so the normal equations are `p x p` with `p` a few
hundred, whatever the sample count. 500,000 individuals cost the same linear
algebra as 500. Cost is dominated by reading cis windows; each gene is one seek
per variant and nothing else in the .bed is touched.

Expect a few genes per second per trait on one core. Use `--max-genes` to time a
subset before committing.

## What this does not do

- It does not test a foundation model. It measures the baseline it must beat.
- It does not compute burden statistics. You supply them, from an analysis that
  never saw the array data used here.
- It makes no causal claim.
- It handles relatedness only through the `group` column you supply.
