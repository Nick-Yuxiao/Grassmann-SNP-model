# v7.7.3 long-range architecture pilot contract (no GPU)

## Scope and what the v7.7.2 PASS does and does not mean

This CPU-only release drafts the prospective long-range architecture pilot that the
signed v7.7.0 factorial contract deferred. It is authorized only because the v7.7.2
execution returned `LONG_RANGE_TASK_VALIDITY_PASS`: the frozen synthetic task
separates "local is insufficient" from "a conventional global route can solve it."
The local linear probe was practically equivalent to the marginal predictor, the
explicit two-source oracle succeeded, the conventional global MLP cleared its
positive control, and the target-shuffled control stayed equivalent.

Task validity is a mechanism and harness property only. It is not a Grassmann
success, not a biological long-range signal, and not population generalization. This
contract does not train Grassmann, does not allocate or query a GPU, and issues no
architecture verdict. It only freezes the design under which a later, separately
authorized comparison could be run, and authorizes drafting the CPU blinded variance
pilot that must precede any formal launch.

The immutable sources are the v7.6.4 A1-R audit, the v7.7.0 Gate contract, the
v7.7.1 harness, and the v7.7.2 task-validity execution. None of them is modified,
rescued, or reclassified here. The registered `A1R_LD_REGIME_DEPENDENT` decision
stands.

## Registered factorial architecture

The pilot fixes the complete 2x2 factorial defined in v7.7.0. Two factors are
crossed: the Grassmann component (absent or present) and a shared global router
(absent or present).

- `LR_A00` absent/absent is the local-only cell.
- `LR_A10` present/absent is the Grassmann-only cell.
- `LR_A01` absent/present is the conventional global router and the primary
  comparator.
- `LR_A11` present/present is Grassmann added to the identical global router and the
  primary candidate.

The shared global router is a conventional non-Grassmann long-range sequence model.
Unlike the v7.7.2 task-validity positive control, the router is not handed the two
source positions; its receptive field is the full 8,192-token sequence and its input
is the full token sequence. Its implementation, receptive field, inputs, optimization
schedule, and training data are identical between the two router-present arms
`LR_A01` and `LR_A11`. This is the property that makes the incremental Grassmann
contrast interpretable.

## Estimands

The primary estimand is the router-conditional incremental Grassmann benefit:

`delta_primary = NLL(LR_A01) - NLL(LR_A11)`
= `NLL(no-Grassmann + router) - NLL(Grassmann + same router)`.

Positive values favor Grassmann. The prospective practical-benefit margin is 0.010
nats per target. A future GO requires a paired one-sided 95% lower confidence bound
strictly greater than 0.010 nats per target, together with task-validity PASS,
target-shuffled-control PASS, and fairness PASS. A numerical improvement or a nominal
p-value alone is insufficient.

The Grassmann main effect averaged over router levels and the Grassmann-by-router
interaction are secondary estimands. They cannot override or substitute for the
primary router-conditional contrast.

## Fairness, non-identity, and matched labels

Absent components require auditable parameter and compute compensation. The labels
`matched-parameter` and `matched-compute` are forbidden until both a non-identity
audit and a realized-compute audit pass. The non-identity audit must show that a
compensated absent-component arm is not numerically equivalent to a trivial
pass-through of its present-component counterpart. The realized-compute audit must
confirm that the compensation reflects realized, not nominal, compute, reusing the
v7.6 realized-compute fairness discipline.

## Replication and the blinded variance pilot

The independent algorithmic repeat unit is the synthetic truth seed. Multiple
initialization seeds within one truth seed are averaged before inference and are not
independent replicates. All arms share paired truth, masking, split, and
initialization schedules within a block.

The formal truth-seed count is not chosen here. This contract authorizes drafting a
CPU-only, blinded, planning-only variance pilot whose sole outputs are the dispersion
of `delta_primary` across truth seeds and the frozen formal sample size derived from
a power calculation targeting detection of the 0.010 margin at the registered
confidence. The pilot runs on explicit CPU tensors at a reduced proxy scale. Arm
means are withheld; only variance and the resulting seed count may be released. Any
data used by the pilot is barred from the formal analysis. No result-dependent task
family, distance, threshold, model, seed, or stopping-rule selection is permitted.

## Stop and transition rules

- If the future blinded variance pilot cannot reach adequate precision within a
  feasible truth-seed count and compute budget, stop with no GPU launch and no
  architecture verdict.
- If a non-identity or realized-compute audit fails, stop with no GPU launch; the
  matched labels remain forbidden.
- If the later preregistered paired comparison is adequately precise and Grassmann
  fails to show a practical accuracy, efficiency, or transfer benefit against the
  conventional global comparator, close the v7 Grassmann-primary route rather than
  search additional outcome-driven tasks.

The only permitted transition from this contract is
`IMPLEMENT_V7_7_4_CPU_BLINDED_VARIANCE_PILOT_NO_GPU`. It does not authorize Grassmann
training beyond a blinded CPU proxy, a GPU factorial, A3, HAPNEST, HGDP, phenotype
work, or any biological long-range claim.
