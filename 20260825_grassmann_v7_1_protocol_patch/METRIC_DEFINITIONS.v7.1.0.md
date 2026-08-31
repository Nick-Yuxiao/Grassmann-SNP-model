# V7 metric definitions v7.1.0

All v7.0.1 numeric thresholds are retained. The changes below freeze the
evaluation population, masking roles and two distinct uncertainty sources.

| Name | Definition | Decision use |
|---|---|---|
| `masked_nll` | Natural-log cross entropy summed over evaluated masked genotype tokens and divided by their count. Padding and unmasked tokens are excluded. | Lower is better. |
| `delta` | `masked_nll_comparator - masked_nll_candidate`, paired within data seed, init seed, evaluation split, mask/rate and fairness regime. | Positive means candidate is better. |
| `delta_min` | 0.010 nats/masked-token. | A primary GO requires both adjusted lower bounds to be strictly above it. |
| `delta_NI` | 0.010 nats/masked-token for `candidate - comparator`. | A3 non-inferiority requires the upper bound below it. |
| `delta_LD` | 0.010 nats/masked-token. | Only for the closed-form LD reference. |
| `dosage_r2` | Squared Pearson correlation between expected dosage `P(g=1)+2P(g=2)` and observed dosage over evaluated masked tokens. | Secondary; zero-variance cells are undefined. |
| `population_equal_delta` | Compute delta within every HGDP population, then take the unweighted mean across populations. | Primary HGDP estimator. |
| `individual_equal_delta` | Pool/equally weight HGDP individuals. | Secondary sensitivity only. |

## Seed confidence interval

Average the two initialization seeds within each of the five data seeds. Compute
five paired `population_equal_delta` values. For each of the two primary mask
families, report a two-sided 97.5% Student-t interval (`df=4`), which is the frozen
Bonferroni simultaneous interval for familywise alpha 0.05. Record all five values,
the quantile implementation and software version. Initialization seeds are not ten
independent replicates.

## Population confidence interval

After averaging the frozen training-seed estimates, resample whole HGDP populations
with replacement 10,000 times using seed 71001. Recompute the equal-population mean
for each replicate and report the two-sided 97.5% percentile interval with explicit
quantile convention. Do not resample tokens or treat individuals from one population
as independent population replicates.

The seed CI measures training/data-generation variability. The population block
bootstrap measures sensitivity to which HGDP populations were observed. Neither is
a substitute for the other. Both must pass the primary gate in both fairness
regimes.

## Missingness and failures

OOM, non-finite values, empty masks, missing population cells, missing paired seeds,
hash mismatches and unplanned retries remain explicit failed runs. No failed row is
deleted or converted to a numeric metric. Stress-sensitivity success cannot replace
a failed or missing primary 0.90 cell.

## Calibration metric

The 38M SNPBag reproduction reports the published-compatible imputation metric on
the exact 216-person subset when possible. Its verdict is PASS, FAIL or
NOT_COMPARABLE. It calibrates the data/evaluation harness and is not pooled into the
8M architecture superiority estimate.
