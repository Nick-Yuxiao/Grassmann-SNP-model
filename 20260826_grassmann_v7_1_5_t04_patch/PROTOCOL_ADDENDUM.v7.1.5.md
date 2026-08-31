# v7.1.5 T04 two-tier compute addendum

Status: frozen before pilot outcomes are inspected.

This addendum records the completed T02/T03 engineering gates and separates a
rapid exploratory pilot from the unchanged confirmatory A1 design.

## Completed gates

- T02 materialization: GO at data-derived `L=154850`; exact allele keys, exact
  sample sets, phased complete GT, no retained INFO annotations.
- T03 profiling: GO on physical GPU 1; all six 100-step cells passed and the
  maximum observed reserved memory was below 14 GiB.

## Rapid pilot

The pilot is a paired randomized-block screening experiment, not a confirmatory
analysis. The block is `(pilot_data_seed, mask_cell)`. Within every block all
three models receive the same generated examples and mask realization. Model
execution order is deterministically rotated from the pilot seed so run order is
not confounded with architecture.

- Models: local attention, local attention plus ancestry PC, and Grassmann.
- Mask cells: random 0.50, random 0.90, LD-block 0.90, and within-chromosome
  long-range 0.90.
- Pilot data seeds: 91001 and 91002; pilot initialization seed: 91011.
- Training budget: 2,000 steps, batch size 1, `L=154850`.
- Evaluation: donor-validation-derived synthetic examples only. HGDP is forbidden
  during pilot training, tuning, stopping, and model selection.
- Random 0.50 is a negative control and cannot create a positive signal verdict.

Preliminary signal continuation requires both structured 0.90 mask cells to have
positive mean paired delta versus the best comparator, with no pilot seed worse
than `-0.010 nats/masked-token`. Preliminary futility requires both structured
cells to have mean delta at most `-0.010` and both pilot seeds negative in each.
All other patterns are `PILOT_INCONCLUSIVE`. No confidence interval or scientific
superiority claim is made from two pilot data seeds.

## Confirmatory A1 remains unchanged

The pilot seeds and pilot outputs are excluded from confirmatory inference. A1
still requires five confirmatory data seeds, two initialization seeds, both
fairness regimes where frozen, and the prespecified population and seed confidence
procedures. Full A1 capacity is not treated as reserved merely because GPUs were
idle during T03.
