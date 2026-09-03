# v7.7.3 long-range architecture-pilot contract

## Purpose and authority

This release is a planning-only successor to the immutable v7.7.2 task-validity
execution. Its source status must be `LONG_RANGE_TASK_VALIDITY_PASS`. It defines
the architecture-pilot estimands, blocking, controls, and future fairness audit;
it neither implements nor launches a learned Grassmann arm. GPU use, architecture
ranking, A3 scaling, HAPNEST, HGDP, phenotype work, and biological claims remain
forbidden.

The v7.6 A1-R result remains `A1R_LD_REGIME_DEPENDENT` and is not reclassified.
The present synthetic task is a new mechanism experiment, not a rescue analysis.

## Frozen 2x2 pilot design

The two factors are `global_router` (absent/present) and `grassmann_component`
(absent/present). The four cells are:

| Cell | global router | Grassmann component | role |
| --- | --- | --- | --- |
| `R0_LOCAL` | absent | absent | local-negative diagnostic only |
| `R1_ROUTER` | present | absent | primary conventional comparator |
| `R2_GRASSMANN_ONLY` | absent | present | diagnostic only |
| `R3_ROUTER_GRASSMANN` | present | present | primary target arm |

The primary future contrast is exactly
`NLL(R1_ROUTER) - NLL(R3_ROUTER_GRASSMANN)`; positive values favor the added
Grassmann component conditional on an otherwise shared global router. Main effects
and interaction are secondary and cannot override the primary contrast.

For router-present cells, source positions, token embedding, router topology,
training split, target mask, optimizer family, learning-rate schedule, batch
construction, and evaluation code must be identical. The future implementation
must expose the Grassmann branch as an additive residual around that shared router;
no extra raw input, oracle feature, different target, or different validation set
may enter only the Grassmann arm.

## Pilot, replication, and blinding

This contract freezes a planning pilot of six new, disjoint synthetic truth seeds
`77301` through `77306`, with init seeds `87401` and `87402` nested within each
truth seed. Thus the future pilot has 48 executions but effective n=6. Seeds used
in v7.7.2 (`77201`–`77205`) are excluded from all pilot and formal analyses.

The pilot is used only to audit functional nonidentity, realized parameter counts,
realized CPU/GPU compute, and paired dispersion. Before any formal launch its
aggregate arm means, rankings, p-values, and GO/NO-GO labels must remain blinded.
Pilot data cannot be folded into the later formal inference.

## Required gates before any pilot launch

The next implementation-only release must provide all of the following without
launching a GPU job:

1. a deterministic four-cell harness with the exact factor labels above;
2. a static audit showing identical router-present input and router identities;
3. a nonidentity audit showing the Grassmann residual is computationally active
   and not an identity/no-op;
4. parameter and per-step-compute profiling sufficient to derive, rather than
   assume, future matched-parameter and matched-compute schedules;
5. target-shuffled versions of both router-present cells;
6. a blinding firewall that releases only pilot dispersion and engineering data.

The later formal sample size must be based on a predeclared practical margin of
0.010 nats/target and an independent truth-seed precision calculation. It must not
use an observed pilot mean effect to choose n.

## Stop rules

If v7.7.2 source hashes/status fail, the shared-router identity audit fails, the
Grassmann branch is a no-op, controls fail, or blinding cannot be enforced, stop
with `LONG_RANGE_ARCHITECTURE_PILOT_NOT_READY`. Do not launch a GPU pilot and do
not replace the task after seeing arm outcomes. Only a separately signed launch
authorization may permit a future pilot.

The sole transition permitted by this release is
`IMPLEMENT_V7_7_4_LONG_RANGE_PILOT_HARNESS_NO_LAUNCH`.
