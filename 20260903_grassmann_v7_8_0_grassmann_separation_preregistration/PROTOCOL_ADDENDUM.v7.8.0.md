# v7.8.0 long-range Grassmann separation pre-registration (CPU-first, GPU-gated)

## Purpose and standing

This release opens a new, PI-initiated affirmative track parallel to the
conventional-only difficulty track (v7.7.8). It is the first stage that formally
admits a Grassmann arm. It pre-registers a separation experiment before any Grassmann
result is seen. It trains nothing, allocates no GPU, and issues no verdict; it freezes
the estimand, controls, fairness audits, and the gate that alone could later authorize
GPU-scale training.

Despite the minor-version bump, this release authorizes only a CPU-first separation
pre-registration and its GPU gate. It does not authorize A3 scaling, HAPNEST, HGDP,
phenotype work, or any biological long-range claim. The v7.6.4 `A1R_LD_REGIME_DEPENDENT`
decision and all v7.7.x artifacts remain immutable.

## Why a separate estimand from v7.7.0

The v7.7.0 primary contrast is the *incremental* Grassmann benefit conditional on a
shared global router that already works, so it requires LR1 (a conventional router
that beats local). The v7.7.6 execution plus the v7.7.7 feasibility findings showed a
fairly trained conventional baseline is bimodal on these distal tasks: it either
solves the task (no incremental room) or cannot beat local (LR1 fails). The
incremental estimand cannot capture the case the PI raised:

> a distal signal that a fair, adequately resourced conventional architecture cannot
> extract, but Grassmann can.

That case is a *separation*, not an increment. This release pre-registers it. A
positive separation, if it survives fair matching, is the affirmative A1 result (not
A1-R): it is the logical complement of A1-R, which found conventional dominant or
equivalent under matched compute.

## The A1-R trap this pre-registration must not re-create

It is trivial to design a synthetic task only Grassmann solves; that is task-shopping
and is scientifically empty. A separation counts only if all of the following hold:

- the conventional comparator is not a strawman: a suite of well-tuned conventional
  long-range families (attention, state-space, dilated-convolution, mlp-mixer), each
  trained to convergence at a parameter and realized-compute budget greater than or
  equal to the Grassmann arm;
- a convergence / compute-sufficiency audit shows each conventional arm has plateaued,
  so its failure is intrinsic to the architecture, not under-training;
- the task family is motivated a priori by the genotype / LD low-rank haplotype
  subspace geometry that motivated Grassmann modelling in the first place, and is
  frozen before any Grassmann run; it is never reverse-engineered from Grassmann's
  own inductive bias;
- oracle solvability, local insufficiency, and target-shuffled controls all pass;
- results replicate across truth seeds with paired confidence intervals, the Grassmann
  arm means are blinded until validity and fairness are locked, and no
  result-dependent task, architecture, seed, or stopping-rule selection is used.

## Pre-registered estimand and gate

Per truth seed, using held-out NLL in nats per target:

- `NLL_oracle` (analytic floor), `NLL_local` (local probe), `NLL_marginal = ln 2`;
- `NLL_conv_best` = minimum over the fair conventional suite;
- `NLL_grass` = Grassmann arm;
- separation gap `S = NLL_conv_best - NLL_grass` (positive favors Grassmann).

Task validity (required, LR1 deliberately NOT required):

- `LR0`: `NLL_marginal - NLL_local` paired two-sided 90% CI inside `[-0.010, +0.010]`;
- oracle-solvable: `NLL_local - NLL_oracle` paired one-sided 95% LCB `> 0.050`;
- target-shuffled: paired two-sided 90% CI inside `[-0.010, +0.010]`.

Primary separation GO (the strong "only-Grassmann" form):

1. conventional genuinely fails under a fair budget: `NLL_conv_best - NLL_local`
   paired two-sided 90% CI inside `[-0.010, +0.010]` (the best fair conventional arm
   is practically equivalent to local, i.e. extracts no distal signal); and
2. Grassmann solves it: `NLL_local - NLL_grass` paired one-sided 95% LCB `> 0.050`
   and per-seed `NLL_grass - NLL_oracle <= 0.050`; and
3. all fairness audits pass (matched-parameter, matched-compute, realized-compute,
   conventional convergence/compute-sufficiency, non-identity), all controls pass,
   and the effect is seed-stable.

Secondary separation signal (weaker, cannot substitute for the primary): `S` paired
one-sided 95% LCB `> 0.050` even when the conventional suite partially works.

## GPU gate and stop rules

- GPU-scale Grassmann training is authorized only by a later, separately signed stage,
  and only if the CPU-first primary separation GO passes with all fairness and control
  gates. A separation absent at CPU proxy is not rescued by scale.
- If, under a fair strong conventional suite that has demonstrably converged, no
  separation appears at CPU proxy, this is the A1-R pattern again: close the v7
  Grassmann-primary route rather than expand the task search or escalate to GPU.
- `matched-parameter` and `matched-compute` labels remain forbidden until the audits
  pass.

The only permitted transition is
`IMPLEMENT_V7_8_1_CPU_GRASSMANN_SEPARATION_HARNESS_NO_LAUNCH`, which requires a
CPU-runnable Grassmann module and the fair conventional suite, and still does not
authorize GPU.
