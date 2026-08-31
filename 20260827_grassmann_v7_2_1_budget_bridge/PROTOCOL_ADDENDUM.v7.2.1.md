# v7.2.1 budget-bridge addendum

Status: frozen after v7.2.0 selected shared peak LR 4e-4 and before any
budget-bridge optimization step.

## Purpose

The v7.2.0 pilot selected a shared peak LR but did not determine the duration of
the stable training phase or the total duration of a later cosine schedule. This
bridge measures time-to-target and terminal dynamics on a common evaluator.

Old v7.1.13 C0 absolute NLL and v7.2.0 pilot absolute NLL are not matched because
their validation individuals, mask seeds, initialization seeds and warmup differ.
They may not be used to compute an acceleration ratio.

## Frozen paired design

- LRs: 1e-4 control and selected 4e-4.
- Models: local attention, local attention plus global PCs, Grassmann full.
- Masks: ld_block_0p90 and within_chrom_longrange_0p90.
- Shared mask seed 92001 and initialization seed 82001.
- Same 2,247 training individuals, deterministic sample stream, mask bank and
  historically viewed LR-tuning validation firewall for both LRs.
- Start all 12 runs from step zero. Resume from v7.1 checkpoints is forbidden.
- 500-step linear warmup from 10% to peak, followed by constant LR.
- First stage ends at step 20,000. Evaluation interval is 250 steps.
- Physical GPUs 1,3,4,5,6,7 each block one model/mask cell and run both LRs.
  LR run order is counterbalanced across the six GPU blocks.
- Save complete model and optimizer checkpoints at step 20,000.
- Record train NLL, validation NLL, LR, gradient norm, throughput and memory.

No decision-holdout individual and no HGDP individual may be read. Architecture
ranking and GO/NO-GO are forbidden.

## Reproducibility and acceleration

At step 4,000, tail-five validation NLL from each bridge run must reproduce the
matching v7.2.0 pilot run within 0.001 nats/masked-token.

Within each model/mask cell, monotonically smoothed 1e-4 and 4e-4 curves are
compared at five evenly spaced interior NLL targets in their common reachable
range. Acceleration is time-to-target at 1e-4 divided by time-to-target at 4e-4.
Cell-specific ratios and their target range are primary; no cross-cell absolute
NLL comparison is permitted.

## Step-20k decision

For every 4e-4 curve, calculate tail-five-mean validation-NLL changes over the
three non-overlapping intervals 14k-16k, 16k-18k and 18k-20k.

- `BUDGET_BRIDGE_20K_ADEQUATE`: every absolute interval change is at most 0.002,
  no terminal tail mean exceeds the best historical tail mean by more than 0.002,
  and all reproducibility checks pass. The proposed formal schedule is then
  warmup 500, stable peak through 20k, cosine decay to 0.1 peak through 30k.
- `BUDGET_BRIDGE_EXTEND_ALL_TO_30K`: reproducibility and stability pass, but at
  least one 4e-4 curve has an absolute interval change above 0.002. Any extension
  must resume all 12 runs, retain constant LR, and use the frozen step-20k
  checkpoints. Selective extension is forbidden.
- `BUDGET_BRIDGE_REPLAN_INSTABILITY`: reproducibility fails, any run fails, or a
  terminal tail mean is more than 0.002 worse than its best historical tail mean.

The bridge does not itself launch an extension. This prevents GPU reacquisition
from occurring without a fresh non-interference audit.

