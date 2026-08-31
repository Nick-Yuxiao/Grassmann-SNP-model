# v7.1.6 primary-first multi-GPU addendum

Status: frozen before any pilot model-comparison outcome is inspected.

This addendum supersedes the unexecuted v7.1.5 pilot execution grid. It does not
change the confirmatory estimand, margins, masks, fairness definitions, seed
replication, or inference procedures. It changes execution order and capacity
allocation so the primary decision is obtained before non-decisive diagnostics.

## Stage C0: convergence pilot

- Use `L=154850`, batch size 1, the three frozen models, the two structured
  0.90 masks, pilot data seeds 91001/91002, and pilot init seed 91011.
- Evaluate donor-validation synthetic examples every 250 steps through step
  10,000. HGDP access is forbidden.
- Candidate common budgets are 4,000, 6,000, 8,000 and 10,000 steps.
- At a candidate budget K, calculate each curve's mean validation NLL in the
  final five evaluations at K and at K-2,000. K is acceptable only when the NLL
  improvement over that interval is at most 0.002 nats/masked-token for every
  model, mask and pilot data seed, with finite losses throughout.
- Select the smallest acceptable K. The rule may inspect within-curve convergence
  only; it must not inspect between-model deltas when selecting K.
- All Primary A1 runs use the same selected K. Per-run or per-model early stopping
  is forbidden. If no candidate is acceptable, Primary A1 is blocked for re-plan.

The C0 comparison at the selected K may be reported as preliminary signal only;
two pilot data seeds do not support confirmatory inference.

## Stage P1: Primary A1 first

Run exactly:

`2 structured masks x 2 fairness regimes x 3 models x 5 data seeds x 2 init seeds = 120 runs`.

The sequence length remains 154,850. Both `matched_parameter` and
`matched_compute` remain required. The primary GO/NO-GO/INCONCLUSIVE rules and
simultaneous confidence procedures remain those frozen in v7.1.0.

Diagnostic random/contiguous cells and 0.99 stress sensitivities are deferred.
They cannot rescue a failed or inconclusive Primary A1 gate. After a definitive
Primary result they may be run for mechanism and completeness, but are not
required to obtain the Primary decision. Any confirmatory extension after an
inconclusive result requires a new signed protocol.

## Multi-GPU execution contract

- Allowed physical GPUs: 1, 2, 3, 4, 5 and 6. Physical GPU 0 is forbidden.
- Every GPU must pass a read-only idle check immediately before a run and must
  have a non-blocking project lock. Busy GPUs are never preempted or signalled.
- The paired block is `(mask, fairness, data_seed, init_seed)`. The three model
  runs in a block use the same physical GPU and shared frozen data/mask inputs.
- Blocks are deterministically shuffled and balanced across GPUs. Model order is
  rotated within blocks so architecture is not confounded with run order.
- If fewer than six GPUs are idle, work waits; it is not silently reassigned in a
  way that changes the frozen schedule.
