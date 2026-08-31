# v7.2.4 final optimization-horizon diagnostic

Status: prospectively frozen after the immutable v7.2.3 run returned
`BUDGET_EXTENSION_30K_NOT_ADEQUATE_REPLAN` and before any step above 30,000.

## Purpose and evidence boundary

The 30k audit completed all 12 continuations without failure or lineage error. Four of
six selected-LR cells passed the terminal-change Gate. The two failures were
`grassmann_full_8m_w256|within_chrom_longrange_0p90` and
`local_attn_8m_w256|within_chrom_longrange_0p90`; their training and validation NLLs
both continued to decrease. Because the local baseline had the largest terminal drop,
reading between-model differences at 30k could favor Grassmann through differential
undertraining.

This release performs one final, bounded optimization-horizon diagnostic. It cannot
read architecture deltas, decision holdout data, HGDP, formal A1-R outcomes or any
phenotype label. It cannot produce an architecture GO/NO-GO decision.

## Frozen continuation design

- Resume all six selected-LR (`4e-4`) model-by-mask cells from their exact immutable
  v7.2.3 `CHECKPOINT_STEP30000.pt` files.
- The complete `3 models x 2 masks` factorial is retained. Continuing only the two failed cells is forbidden because it would make training horizon depend on observed
  model identity and mask outcome.
- The six `1e-4` cells are excluded prospectively. Their sole role was descriptive LR
  calibration and common-support acceleration estimation, both completed in v7.2.3;
  they are not part of the selected formal optimization protocol.
- Preserve optimizer state, model, mask, mask seed 92001, initialization seed 82001,
  deterministic sample stream, mask-bank phase and validation firewall.
- Continue constant LR `4e-4` from step 30,001 through step 40,000, evaluate every
  250 steps and save complete model/optimizer/RNG checkpoints at 40k.
- Use physical GPUs 1,3,4,5,6,7, with one frozen model-mask cell per GPU. GPUs 0 and 2
  are forbidden. Every GPU requires a fresh read-only audit and advisory lock.
- This is the final automatic horizon diagnostic. No run may continue beyond 40k
  without a new protocol, compute contract and independently justified estimand.

## Step-40k audit

For every cell, calculate tail-five-mean validation-NLL changes over 34k-36k,
36k-38k and 38k-40k. The primary operational adequacy rule remains:

- every interval has absolute change at most 0.002 nats/masked-token;
- best-to-terminal degradation is at most 0.002;
- all six continuations and exact lineage checks pass.

The prospective non-primary shape classification is unchanged:

- `STABLE`: the primary interval criterion passes and either not all drops are
  positive or the last/first absolute-change ratio is at most 1.5;
- `STABLE_BUT_ACCELERATING`: all three drops are positive and individually at most
  0.002, but the last/first ratio is greater than 1.5;
- `NOT_STABLE`: at least one interval exceeds 0.002.

`STABLE_BUT_ACCELERATING` blocks automatic schedule freeze even though the primary
absolute-change rule passes. `NOT_STABLE` is a primary failure. These are operational
training-adequacy rules, not statistical proof that a global optimum was reached.

## Frozen decision branches

- `FINAL_BUDGET_40K_ADEQUATE`: all six cells pass, no shape flag and lineage is
  intact. This authorizes a later formal schedule contract with warmup 500, constant
  `4e-4` through 40k and cosine decay to `4e-5` through 50k. It does not itself start
  that schedule or authorize architecture comparison.
- `FINAL_BUDGET_40K_PRIMARY_ADEQUATE_SHAPE_REVIEW`: all primary interval criteria
  pass but at least one selected cell is `STABLE_BUT_ACCELERATING`. Stop for review;
  no formal schedule is frozen automatically.
- `FINAL_BUDGET_40K_NOT_ADEQUATE_STOP`: lineage and stability checks pass, but one or
  more cells fail the terminal-change criterion. No further constant-LR extension is
  authorized. Replan options are a separately signed tail-slope/uncertainty protocol
  or a descriptive-only real-1KGP A1-R scope.
- `FINAL_BUDGET_REPLAN_INSTABILITY`: run, sequence, lineage or degradation failure.
  Stop and diagnose.

No branch may convert a sample-limited negative result into architecture NO-GO.
