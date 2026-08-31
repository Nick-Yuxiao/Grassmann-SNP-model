# v7.1.7 real-1KGP A1-R addendum

Status: frozen before any A1-R model-comparison outcome is inspected.

This addendum supersedes the unexecuted synthetic-data semantics in v7.1.6 for
the first A1 run. It does not replace the external-HGDP confirmatory experiment.
It defines a real-panel preliminary gate (`A1-R`) using the already frozen 1KGP
donor partitions and chr22 sites.

## 1. Scope and data isolation

- Train only on `DONOR_TRAIN_2247.bcf` and evaluate model selection and A1-R
  metrics only on `DONOR_VALIDATION_249.bcf`.
- Sequence length is the frozen data-derived `L=154850`.
- HGDP access is forbidden throughout A1-R training, convergence selection,
  evaluation, inference and troubleshooting. HGDP remains sealed for a later
  external experiment.
- The result is limited to comparative masked-genotype prediction on the real
  1KGP chr22 panel. It is not a population-independent or biobank claim.

## 2. Repetition semantics

The frozen panel is one biological dataset. The five former `data_seed` values
are renamed `mask_seed`. They generate shared masks and minibatch order for the
three paired model arms; they are algorithmic repeats, not five independent
biological datasets. Two initialization seeds are averaged within mask seed.

The donor-validation population bootstrap remains a separate uncertainty view.
It resamples whole 1KGP populations and never tokens or masks. Neither uncertainty
view may be described as replication across independent cohorts.

## 3. Nested donor-size diagnostic

Before outcomes are inspected, build deterministic population-stratified nested
subsets from the 2,247 donor-training individuals:

- 25%: 562 individuals;
- 50%: 1,124 individuals;
- 100%: all 2,247 individuals.

Within each population, samples are ranked by a frozen SHA-256 key. Hamilton
largest-remainder allocation gives the exact global target, and the 25% set must
be a subset of the 50% set, which must be a subset of the 100% set. All IDs,
population counts and file hashes are frozen before training.

The 100% primary grid remains 120 runs. The 25% and 50% grids are diagnostic:
each uses the three models, two primary masks, `matched_compute`, mask seeds 1/2
and init seed 1, for 12 runs per fraction. The corresponding 100% cells already
in the primary grid are the size-curve anchors. Total post-pilot runs are 144.

## 4. Convergence

Stage C0 chooses one common fixed step budget K from 4k/6k/8k/10k using only
within-curve donor-validation NLL. No between-model delta may be inspected while
choosing K. Per-run early stopping is forbidden.

In addition, every one of the 120 primary runs must pass the same frozen tail
convergence audit at K. A negative point estimate from a non-converged or missing
primary run cannot support NO-GO.

## 5. Three-state A1-R decision

The decision comparison is `grassmann_full_8m_w256` versus the frozen
`local_attn_gpc_8m_w256` comparator. `local_attn_8m_w256` remains a mechanistic
control. Delta is comparator NLL minus Grassmann NLL; positive favors Grassmann.
The minimum meaningful advantage remains `delta_min=0.010` nats/masked-token.

For both structured masks and both fairness regimes, compute the frozen 97.5%
simultaneous mask-seed t interval and the 97.5% population-block-bootstrap
interval on donor validation.

- `A1R_PRELIMINARY_GO`: every required lower bound is strictly greater than
  `+0.010`.
- `NO_GO_EQUIVALENT_OR_WORSE`: all 120 primary runs are complete and converged;
  every required upper bound is at most `+0.010`; all donor-size diagnostic cells
  are complete; and no prespecified positive donor-size trend is present.
- `INCONCLUSIVE_SAMPLE_LIMITED_HAPNEST`: every other valid non-GO result,
  including wide intervals, missing convergence, incomplete cells, or a positive
  donor-size trend. This result activates the HAPNEST data-expansion branch and
  does not terminate the project.

The diagnostic positive-trend flag is raised when, for either primary mask, the
paired Grassmann-minus-local+PC advantage increases from 25% to 50% and again
from 50% to 100%, with total increase at least 0.002 nats/masked-token. This is a
conservative branch trigger, not confirmatory evidence of a scaling law.

HAPNEST results are a new data-distribution experiment. They cannot retrospectively
rescue, relabel or overwrite the real-panel A1-R result without a new protocol.

## 6. GPU and execution policy

- Allowed physical GPUs are 1,2,3,4,5,6. Physical GPU 0 is forbidden.
- Every launch requires a read-only idle audit and a non-blocking project lock.
- Existing processes are never signalled, killed or preempted.
- The three model arms in a paired block use the same physical GPU and frozen
  inputs. Diagnostic 25/50/100 blocks sharing a size-curve ID use the same GPU.

