# Informal design-feasibility protocol

_Prospectively frozen synthetic sensitivity analysis; permanently excluded from the formal v6.1 evidence chain._

---

## 📋 Question

Before investing further in calibration, determine whether the frozen rank-2 estimand is plausibly estimable and detectably separated from null geometry under planning ranges for sample size, MAF, population eigengap, and true maximum principal angle.

This is a sensitivity analysis of a proposed design, not post-hoc observed power. It uses no real phenotype, formal candidate list, formal outcome p-value, GC seed, or T14/T16 configuration.

## 🎯 Synthetic design

For each independent replicate, target dosage follows `Binomial(2, MAF)`. Region features have a frozen AR(1)-like conditional-LD covariance that varies modestly across dosage groups. Five outcomes follow the frozen joint conditional model

```text
Y = C alpha + G a + X B + (G * X) Gamma + E,
M_g = B + g Gamma.
```

`M_0` and `M_2` are constructed with a known rank-2 leading subspace angle and a third singular component that sets the population relative gap. `M_1` is the linear midpoint implied by the model. Residuals are heteroskedastic across dosage groups. The same frozen ridge estimator and outcome whitening used by v6.1 are fitted to every replicate.

The crossed factors are sample size, MAF, population relative gap, and true maximum principal angle. The design uses 12 independent family replicates per cell. Replicate indices are paired across angle levels within `(n, MAF, gap)` only to reduce Monte Carlo noise in descriptive contrasts; the family remains the independent unit.

## ⚙️ Estimability and detectability proxy

A replicate is jointly estimable only if:

1. every dosage group contains at least 50 subjects; and
2. the minimum fitted relative rank gap across `g=0,1,2` is at least 0.10.

The fitted direction score is set to zero for non-estimable replicates. Within each `(n, MAF, population gap)` stratum, the 95th percentile of the angle-zero gated scores is a Monte Carlo null reference. The exceedance rate at nonzero angles is reported as `null_exceedance_rate`.

This rate is not a p-value or a calibrated power estimate. It has only 12 null replicates, does not use synchronized wild bootstrap, does not control a candidate family, and cannot support a method comparison.

## 📊 Frozen planning classification

The global planning subset contains cells with population gap at least 0.10 and true angle at least 10 degrees. A cell is called workable when both `joint_estimable_rate >= 0.60` and `null_exceedance_rate >= 0.50`.

| Workable planning-cell fraction | Informal label | Permitted response |
| ---: | --- | --- |
| `< 0.25` | `STOP_REVIEW` | Reconsider the estimand or study scale before more calibration work |
| `0.25–<0.50` | `CAUTION` | Continue only with explicit feasibility limitations |
| `>= 0.50` | `PROVISIONALLY_FEASIBLE` | Continue T09 design; no formal method claim |

Equal weighting of planning cells is not an estimate of the fraction of real candidates. A future null-blind, manifest-bound candidate distribution is required before any population-level prevalence statement.

## 🔒 Evidence firewall

The design and all outputs are marked `EXPLORATORY_NON_EVIDENCE`. They may trigger stopping or broad study redesign. They may not select a formal sample subset, candidate family, seed, alpha, gap threshold, rank, effect size, baseline, or T14/T16 cell. No result from this package may be promoted by copying it into a formal directory; any later formal analysis requires a new prospective protocol and independent seed block.

