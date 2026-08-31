# v7.2.2 budget-extension addendum

Status: frozen after the immutable v7.2.1 bridge returned
`BUDGET_BRIDGE_EXTEND_ALL_TO_30K` and before any extension optimization step.

## Purpose

The v7.2.1 bridge reproduced all 12 pilot cells, had no run failure and no
instability, but only two of six selected-LR model/mask cells satisfied the
frozen step-20k terminal-change criterion. This extension determines whether
the common optimization horizon is adequate at step 30,000.

## Frozen continuation design

- Resume all 12 v7.2.1 runs from their exact `CHECKPOINT_STEP20000.pt` files.
- Selective continuation, restart from step zero and optimizer reset are forbidden.
- Preserve model, mask, mask seed 92001, initialization seed 82001, learning
  rate, deterministic sample stream, mask-bank phase and validation firewall.
- Continue the constant learning rate from step 20,001 through step 30,000.
- Evaluate every 250 steps and save complete model/optimizer checkpoints at 30k.
- Physical GPUs 1,3,4,5,6,7 retain their original model/mask blocks and paired
  LR order. GPUs 0 and 2 are forbidden.
- Read only the historically viewed LR-tuning validation partition. Decision
  holdout and HGDP access remain forbidden.

The extension is an optimization-budget diagnostic. Architecture ranking,
between-model deltas and GO/NO-GO decisions are forbidden.

## Step-30k audit

For every 4e-4 curve, append the audited extension to its immutable 0-20k
source curve and calculate tail-five-mean validation-NLL changes over 24k-26k,
26k-28k and 28k-30k.

- `BUDGET_EXTENSION_30K_ADEQUATE`: all six cells have absolute change at most
  0.002 in all three intervals, no terminal degradation exceeds 0.002, and all
  12 continuations pass lineage checks. The proposed formal schedule becomes
  warmup 500, stable 4e-4 through 30k, then cosine decay to 0.1 peak through 40k.
- `BUDGET_EXTENSION_30K_NOT_ADEQUATE_REPLAN`: lineage and stability pass, but
  at least one terminal interval exceeds 0.002. No further automatic extension
  is authorized.
- `BUDGET_EXTENSION_REPLAN_INSTABILITY`: any continuation fails, lineage is
  broken, or terminal degradation exceeds 0.002.

No branch launched by this addendum may read the decision holdout or formal
A1-R outcomes.
