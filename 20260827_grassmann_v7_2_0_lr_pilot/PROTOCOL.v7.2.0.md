# v7.2.0 shared learning-rate pilot protocol

Status: frozen before any v7.2.0 optimization run.

## Purpose and boundary

The v7.1.13 extension and v7.1.14 audit established continuing validation-NLL
improvement under constant LR 1e-4, but did not establish an architecture effect
or a staged-learning mechanism. Formal A1-R remains blocked.

This pilot selects one shared peak learning rate for a subsequent C0 schedule
confirmation. It is optimization-only, starts every run from step zero, uses no
HGDP data, and may not support GO/NO-GO or architecture ranking.

## Validation firewall

The historical C0 evaluation rule selected 32 donor-validation individuals for
each of mask seeds 91001 and 91002. Their union is frozen as the LR-tuning set.
Every other donor-validation individual is frozen as the untouched decision
holdout. The LR pilot may read only the tuning indices. Decision-holdout metrics
are forbidden until the formal protocol authorizes them.

Population membership is audited from the immutable v7.1.1 donor-validation TSV.
No historically viewed individual may enter the decision holdout.

## Frozen pilot design

- Sequence length: 154,850.
- Training data: all 2,247 real 1KGP donor-train individuals.
- Models: local attention, local attention plus global PCs, and Grassmann full.
- Masks: ld_block_0p90 and within_chrom_longrange_0p90.
- Shared pilot mask seed: 92001; shared initialization seed: 82001.
- Peak LR grid: 1e-4, 2e-4, 4e-4 and 8e-4.
- Optimizer: AdamW; 500-step linear warmup from 10% to 100% of peak, then
  constant peak through step 4,000.
- One run per LR/model/mask cell: 24 runs, 96,000 run-steps.
- Each physical GPU is a block for one model/mask cell and runs all four LRs;
  GPU-local run order is prespecified and varied across the six blocks.
- Evaluation every 250 steps. Per interval record training masked NLL,
  validation masked NLL, learning rate, gradient norm, throughput and memory.
- Per-run early stopping is forbidden. Non-finite loss/gradient or training NLL
  above 5.0 after warmup is a hard failure retained for audit.

LR 8e-4 is a prespecified stress boundary, not a presumed optimum.

## Shared-LR selection

For each model/mask cell, compare the tail-five mean validation NLL with its own
1e-4 control. Each of the six cells receives equal weight.

- A candidate is ineligible if any run failed or any cell is worse than its
  1e-4 control by more than 0.002 nats/masked-token.
- Candidate score is the mean of the six within-cell gains.
- Among eligible candidates within 0.0005 of the best score, select the lowest LR.
- If 1e-4 is selected, return LR_PILOT_NO_ACCELERATION_REPLAN; do not silently
  claim that the LR hypothesis succeeded.
- Models may not receive different primary LRs. Before any architecture NO-GO,
  an equal-search-budget per-model LR sensitivity check remains required.

## Next boundary

A selected shared peak LR authorizes only a new all-12-run C0 schedule
confirmation from step zero. It does not authorize the 144 formal A1-R runs.
Training adequacy and the n=5 paired-variance pilot must be signed before A1-R.
