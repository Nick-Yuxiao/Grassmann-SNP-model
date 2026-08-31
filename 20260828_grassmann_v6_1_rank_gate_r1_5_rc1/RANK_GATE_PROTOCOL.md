# D29 rank-gated statistic

For candidate k and dataset z, let `D_k(z)` be the frozen rank-2 direction score and
let `q_k(z)` be the minimum relative singular-value gap across g=0,1,2. With the
approved fixed threshold delta=0.10, define

`T_k(z) = D_k(z) * I[q_k(z) >= delta]`.

The same deterministic mapping is applied to the observed data and every synchronized
multiplier bootstrap resample. If the observed candidate is ineligible, its observed
statistic is zero and its candidate p-value is exactly one because all bootstrap
statistics are nonnegative. The candidate is not deleted, so family membership and
multiplicity scope do not change.

This modifies the primary testing procedure but not the `M_g=B+g Gamma` estimand.
The threshold was approved after the bounded-smoke design exposed weak rank; it may
not be retuned against subsequent p-values. R1.5 PASS authorizes only redesigning the
bounded smoke with a new seed block.

