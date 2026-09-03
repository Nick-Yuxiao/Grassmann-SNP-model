# v7.7.7 feasibility findings (CPU due-diligence, NOT part of the signed manifest)

These reduced-proxy CPU smokes (L=96, target/local radius mirrored from the frozen
task, fair attention baseline NOT handed positions, 3 truth seeds) were run to check
whether a legitimate headroom corridor is reachable before freezing the re-plan.
They are engineering evidence, not a formal task-validity run, and are barred from
any formal analysis.

## What "corridor" means

A configuration is usable only if a FAIRLY TRAINED conventional global baseline:
- beats local (LR1): `d_global = NLL(local) - NLL(global) > 0`, and
- leaves headroom (H): `NLL(global) - NLL(oracle) >= 0.020`, and
- passes the target-shuffled control, on the same seeds.

## Results

| regime | local | global | oracle | d_global (beat local) | headroom | verdict |
| --- | --- | --- | --- | --- | --- | --- |
| positions HANDED + parity K=4,6,8 (v7.7.6 canonical) | ~0.69 | ~1e-5 | ~1e-6 | large | ~0 | baseline SOLVES → ceiling |
| positions RANDOM + majority K=3,5,7 | ~0.70 | ~1e-4 | ~0 | large | ~0 | baseline SOLVES → ceiling |
| positions RANDOM + parity K=3 | 0.706 | 0.476 | ~0 | +0.23 | 0.01–0.7 (min 0.0099) | PARTIAL but seed-bimodal; shuffled FAILS |
| positions RANDOM + parity K=5,7 | ~0.70 | 0.82–0.97 | ~0 | negative | n/a | baseline CANNOT beat local → LR1 fails |

(An intermediate smoke where the bit was stored at the marker's neighbor is excluded:
the pooling attention structurally could not read the neighbor, an artifact of the
proxy, corrected by carrying the bit inside the marker token.)

## Structural conclusion

A fairly trained conventional baseline is BIMODAL on these distal tasks:
- when retrieval is easy (positions handed) or the label gives strong graded signal
  (majority), it SOLVES the task → no headroom for any added mechanism;
- when retrieval is hard (random positions) and the label gives no partial-credit
  gradient (parity, larger K), it CANNOT beat local → not a valid test;
- only a narrow, seed-unstable, control-failing sliver sits between.

Consequently the >= 0.010 "practically meaningful incremental benefit" window is
elusive by construction. v7.7.7 therefore freezes a pre-declared grid search for a
fair, seed-stable corridor and, crucially, routes a corridor-empty outcome to closing
the v7 Grassmann-primary route rather than to further outcome-driven task search.
