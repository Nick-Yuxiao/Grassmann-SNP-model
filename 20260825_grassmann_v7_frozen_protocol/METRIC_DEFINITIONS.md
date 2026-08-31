# V7 metric definitions (v7.0.1)

This file is the implementation-facing dictionary for the constants frozen in
`DECISIONS.v7.0.1.tsv`.  Values and directions must not be tuned after looking
at A1/A2/A3 outcomes.

| Name | Definition | Unit | Direction / decision use |
|---|---|---|---|
| `masked_nll` | Sum of natural-log cross entropy over masked genotype tokens divided by the number of evaluated masked tokens. Unmasked and padding tokens are excluded. | nats/masked-token | Lower is better. |
| `delta` | `masked_nll_comparator - masked_nll_candidate`, paired within data seed, init seed, split, mask family and fairness regime. | nats/masked-token | Positive means the candidate is better. |
| `delta_min` | Minimum practically important superiority margin. | 0.010 nats/masked-token | A superiority GO requires the lower confidence bound for `delta` to be strictly greater than this value in every required fairness regime. |
| `delta_NI` | Non-inferiority loss margin for A3. | 0.010 nats/masked-token | With `loss_diff = masked_nll_candidate - masked_nll_comparator`, non-inferiority requires the upper confidence bound to be strictly less than this value. |
| `delta_LD` | Minimum superiority margin against the closed-form LD reference. It has an independent name even though its numeric value equals `delta_min`. | 0.010 nats/masked-token | Use only for the A1 classical LD reference comparison. |
| `overfit_thr` | Maximum tiny-set train NLL compatible with a successful overfit check. | 0.010 nats/masked-token | Engineering gate only; it is not an effect-size threshold. |
| `pc_control_thr` | Minimum advantage of true global PCs over shuffled PCs. | 0.005 nats/masked-token | `NLL_shuffled_PC - NLL_true_PC` must meet the threshold before the PC arm is informative. |
| `leak_accuracy_tol` | Allowed excess of leakage-probe accuracy over its empirical chance baseline. | 0.010 accuracy | Leakage check passes only when excess accuracy is no larger than the tolerance. |
| `resume_param_tol` | Maximum absolute parameter difference between uninterrupted and checkpoint-resumed deterministic runs. | 1e-6 | Engineering pass requires the maximum difference not to exceed the tolerance. |

## Pairing and confidence intervals

The independent inferential unit is the data seed (`n=5` in the confirmatory
tranche). The two initialization seeds are averaged within each data seed
before paired differences are formed. A confidence interval implementation
must record the interval method, confidence level, quantile convention and all
five paired values. Initialization seeds must never be treated as ten
independent replicates.

## Secondary metrics

- `dosage_r2`: squared Pearson correlation between expected dosage
  `P(g=1)+2P(g=2)` and observed dosage over evaluated masked tokens.
- `r2_by_maf_bin`: the same statistic within preregistered MAF bins; bins with
  zero variance are reported as undefined, never as zero.
- `tokens_per_s`: evaluated masked tokens divided by synchronized training or
  inference wall time, with scope explicitly named.
- `peak_vram_gb`: CUDA peak allocated and reserved bytes, reported separately,
  using `2^30` bytes per GiB.
- `wallclock_h`: monotonic elapsed wall time divided by 3600.

## Missing values and failures

OOM, non-finite loss/gradients, empty masks and missing paired cells are failed
runs. They are retained in the run ledger and are not silently retried or
converted to numeric metric values.
