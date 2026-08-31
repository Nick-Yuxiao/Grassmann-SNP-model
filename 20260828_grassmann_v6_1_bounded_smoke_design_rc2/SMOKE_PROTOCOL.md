# Bounded shared-family smoke rc2 protocol

_Prospective logic-smoke protocol; frozen before any rc2 result exists._

---

## 📋 Purpose and evidence boundary

This smoke tests implementation contracts after the archived `GC-screen rc1` failure and the R1/R1.5 repairs. It does not estimate FWER or power. Three family replicates per cell cannot support either claim.

The source experimental unit is one independently generated shared-subject family. Candidate regions, genotype groups, traits, and multiplier resamples are nested objects and are not independent replicates.

## 🎯 Frozen cells

Four null cells cross conditional LD and heteroskedasticity orthogonally. Three controls test target selection, direction-negative amplification, and an obvious rotation signal.

| Cell class | Conditional LD | Heteroskedasticity | Role |
| --- | --- | --- | --- |
| Null 1 | Off | Off | Logic null |
| Null 2 | On | Off | Logic null |
| Null 3 | Off | On | Logic null |
| Null 4 | On | On | Logic null |
| Selected null | On | On | Selection contract |
| Pure amplification | On | On | Direction negative |
| High-SNR rotation | On | On | Direction positive |

Each cell uses three paired replicate indices, `n=360`, four aligned candidates, rank 2, ridge `0.01`, and 39 synchronized Rademacher multiplier resamples. The new seed namespaces are recorded in `config/BOUNDED_SMOKE_CONFIG.json` and do not overlap the R1, R1.5, or rc1 smoke blocks.

## ⚙️ Shared-family and D29 contracts

Every family contains one ordered `(subject_id, Y, G, C)` and four aligned `X_k` matrices. Candidate-specific common-subspace null fits remain separate, while each resample uses one multiplier vector keyed to canonical subject ID and shared across all four candidates.

For every observed candidate and every resample, the statistic is

```text
T_k(z) = D_k(z) * I[q_k(z) >= 0.10].
```

An observed ineligible candidate has statistic zero and candidate p-value one. It remains in the four-candidate family. Bootstrap-ineligible statistics are also zero before maxT aggregation.

The target-selection cell uses subject-disjoint selection and inference samples. The selected target column, not a recorded placeholder index, becomes inference `G`.

## 🔍 Machine checks

`BOUNDED_SMOKE_PASS` requires every prospective check below to be true:

1. Exactly 21 planned families are present, with no failed or duplicate run IDs
2. Every family retains four aligned candidates and 39 multiplier fingerprints
3. All family and candidate p-values are finite and lie on the frozen 1/40 grid
4. D29 is applied identically to observed and resampled statistics
5. Observed ineligible candidates are retained with statistic zero and p-value one
6. One canonical multiplier fingerprint is used across candidates within each resample
7. Selection and inference subjects are disjoint and the selected target is applied
8. Truth labels for null, pure amplification, and rotation controls are correct
9. Null controls are not all degenerate at the minimum p-value
10. Pure amplification is rejected in at most one of three family replicates
11. High-SNR rotation is rejected in at least two of three family replicates

Rank-ineligible candidates are not a runtime failure and are not deleted. Their frequency is reported descriptively because D29 defines the conservative no-call behavior.

## 📊 Run accounting and resource ceiling

The frozen accounting is `7 cells × 3 families × 39 resamples = 819` family-level resamples. With four candidates, this implies 84 observed and 3,276 bootstrap candidate-statistic evaluations.

The run is CPU-only, single-process, and limited to one BLAS thread. Based on the bounded dimensions and local contract-test timing, the prospective expected CPU time is 2–10 minutes. The planning ceiling is 30 CPU-minutes, 1 GiB peak RSS, and 100 MiB output. Exceeding a ceiling produces `RESOURCE_REVIEW_REQUIRED`; it does not authorize silent parallelization or GPU use.

## 🚫 Authorization boundary

Validation and unit tests are authorized locally. Running the smoke requires a new project-owner approval record bound to this package manifest. Even `BOUNDED_SMOKE_PASS` authorizes only drafting a prospective `GC-screen rc2` protocol. It does not authorize running that protocol, formal calibration, power, real data, GPU, or v7 work.
