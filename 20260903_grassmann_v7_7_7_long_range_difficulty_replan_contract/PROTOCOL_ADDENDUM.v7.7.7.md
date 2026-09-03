# v7.7.7 long-range task-difficulty re-plan contract

## Why this release exists

The v7.7.6 CPU execution returned `LONG_RANGE_TASK_DIFFICULTY_UNRESOLVED`
(`selected_k = null`): for every K in {4,6,8} the conventional global baseline
solved the parity task to essentially zero NLL, so the headroom
`H = NLL(baseline) - NLL(oracle)` was ~1e-5 to ~1e-4, far below the 0.020 gate.
The frozen `next_if_unresolved` transition is `STOP_NO_AUTOMATIC_TASK_EXPANSION_OR_GPU`;
the reserved factorial-implementation slot is therefore not unlocked and no GPU,
Grassmann, or architecture work is authorized.

This release is the honest successor to that stop: a CPU-only, planning-only
contract that re-plans task difficulty so a future audit can decide whether a fair
headroom corridor exists at all. It trains nothing, allocates no GPU, and fits no
Grassmann. The v7.6.4 A1-R decision and every prior v7.7.x artifact remain
immutable.

## Feasibility findings that shape the contract

A CPU due-diligence sweep (documented in `FEASIBILITY_FINDINGS.v7.7.7.md`, not part
of the signed manifest) tested whether a fair, adequately trained conventional
global baseline can land in the corridor "beats local AND leaves >= 0.020 headroom
below the oracle." Four regimes were probed:

1. positions handed + parity (the v7.7.6 design): baseline solves it, headroom ~0.
2. positions randomized + graded majority label: a fair attention baseline still
   solves it, headroom ~0.
3. positions randomized + parity, small K: baseline is partial but seed-bimodal and
   the target-shuffled control fails (overfitting), so no stable, valid corridor.
4. positions randomized + parity, larger K: baseline cannot beat local, so LR1 fails.

The structural lesson is that a fairly trained conventional baseline is bimodal on
these distal tasks: it either solves the task (no headroom for any added mechanism)
or fails to beat local (not a valid test), with only an unstable, control-failing
sliver between. A legitimate corridor is not obtained by simply raising parity order
or by handing the baseline the source positions.

## Frozen re-plan design

The re-planned difficulty family must satisfy all of:

- the conventional global baseline is NOT handed the informative source positions;
  it must retrieve them over the full (reduced-proxy) sequence;
- the informative bit is carried in the marker token itself, so the baseline can in
  principle read it;
- the conventional baseline is fairly and adequately trained under a frozen,
  architecture-neutral compute/data budget; headroom produced by an under-trained
  baseline is illegitimate and forbidden;
- the label function and retrieval geometry are chosen only from a pre-declared grid.

## Corridor definition (unchanged thresholds, now seed-stable)

For the selected configuration, across the independent truth seeds:

- `LR0`: `delta_local = NLL(marginal) - NLL(local)` paired two-sided 90% CI wholly
  inside `[-0.010, +0.010]`;
- `LR1`: `delta_global = NLL(local) - NLL(conventional global)` paired one-sided 95%
  lower bound strictly greater than `0.010`;
- `headroom`: `H = NLL(conventional global) - NLL(oracle)` with per-truth-seed
  minimum `>= 0.020` and paired one-sided 95% lower bound strictly greater than
  `0.010`;
- `target-shuffled`: `delta_shuffled` paired two-sided 90% CI wholly inside
  `[-0.010, +0.010]`.

A configuration is corridor-eligible only if all four hold on the same seeds.

## Fairness, blinding, and anti-shopping

- The difficulty grid, label-function family, retrieval geometry, baseline budget,
  and selection rule are pre-declared and frozen before any run.
- Selection uses only the conventional baseline and the oracle. No Grassmann branch
  is implemented, fitted, or consulted at any point in this line.
- `matched-parameter` and `matched-compute` labels remain forbidden until a
  non-identity audit and a realized-compute audit pass.
- All re-plan/pilot data is barred from any later formal analysis.

## Stop and transition rules

- The only permitted transition is
  `IMPLEMENT_V7_7_8_LONG_RANGE_DIFFICULTY_REPLAN_HARNESS_NO_LAUNCH`, which will
  implement the pre-declared grid search under the fairness and control rules and
  select at most one corridor-eligible configuration.
- If the pre-declared grid yields no fair, seed-stable, control-passing corridor
  configuration, this is not a task-engineering shortfall to be patched by more
  searching. It triggers the registered conclusion that, on fair synthetic distal
  tasks, there is no measurable `>= 0.010` window for an incremental Grassmann
  benefit, and the v7 Grassmann-primary route is closed rather than expanded.
- No GPU, Grassmann training, A3, HAPNEST, HGDP, phenotype, or biological
  long-range work is authorized by this release.
