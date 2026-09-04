# Gate 0A detectability — planning proxy

_No real data. `sigma` = replicate-to-replicate SD of the paired R2_genetic difference; it is a planning axis a real-panel pilot must pin. Decision is one-sided (CI lower bound > 0), so FPR@0 for a calibrated 95% CI is ~0.025._

Target: **P(CI lower>0 | Δ=0.005) ≥ 0.9**  (M=3000, B=1500, seed=20260904)

## sigma = 0.005

| Replicates R | FPR @ 0 | Power @ .002 | Power @ .005 | Power @ .010 |
| ---: | ---: | ---: | ---: | ---: |
| 20 | 0.033 | 0.482 | 0.992 | 1.000 |
| 50 | 0.036 | 0.807 | 1.000 | 1.000 |
| 100 | 0.030 | 0.977 | 1.000 | 1.000 |
| 200 | 0.029 | 1.000 | 1.000 | 1.000 |

**Min R for ≥0.90 power at Δ=0.005:** 20

## sigma = 0.01

| Replicates R | FPR @ 0 | Power @ .002 | Power @ .005 | Power @ .010 |
| ---: | ---: | ---: | ---: | ---: |
| 20 | 0.033 | 0.177 | 0.637 | 0.992 |
| 50 | 0.036 | 0.316 | 0.945 | 1.000 |
| 100 | 0.030 | 0.536 | 0.998 | 1.000 |
| 200 | 0.029 | 0.810 | 1.000 | 1.000 |

**Min R for ≥0.90 power at Δ=0.005:** 50

## sigma = 0.02

| Replicates R | FPR @ 0 | Power @ .002 | Power @ .005 | Power @ .010 |
| ---: | ---: | ---: | ---: | ---: |
| 20 | 0.033 | 0.089 | 0.241 | 0.637 |
| 50 | 0.036 | 0.122 | 0.439 | 0.945 |
| 100 | 0.030 | 0.182 | 0.720 | 0.998 |
| 200 | 0.029 | 0.290 | 0.937 | 1.000 |

**Min R for ≥0.90 power at Δ=0.005:** 200

## sigma = 0.03

| Replicates R | FPR @ 0 | Power @ .002 | Power @ .005 | Power @ .010 |
| ---: | ---: | ---: | ---: | ---: |
| 20 | 0.033 | 0.064 | 0.143 | 0.373 |
| 50 | 0.036 | 0.086 | 0.233 | 0.664 |
| 100 | 0.030 | 0.107 | 0.410 | 0.914 |
| 200 | 0.029 | 0.150 | 0.655 | 0.996 |

**Min R for ≥0.90 power at Δ=0.005:** > 200 (none in grid; increase R)

## sigma = 0.05

| Replicates R | FPR @ 0 | Power @ .002 | Power @ .005 | Power @ .010 |
| ---: | ---: | ---: | ---: | ---: |
| 20 | 0.033 | 0.053 | 0.089 | 0.177 |
| 50 | 0.036 | 0.063 | 0.122 | 0.316 |
| 100 | 0.030 | 0.064 | 0.182 | 0.536 |
| 200 | 0.029 | 0.078 | 0.290 | 0.810 |

**Min R for ≥0.90 power at Δ=0.005:** > 200 (none in grid; increase R)
