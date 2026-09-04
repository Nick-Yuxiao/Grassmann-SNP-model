# Gate 0A detectability protocol (planning proxy)

_Prospectively frozen. Uses no real data. Sizes the Gate 0A replicate count `R` at the fixed primary margin 0.005. Date: 2026-09-04._

---

## 📋 Question

For the Gate 0A decision procedure, how many simulation replicates `R` are needed to reliably detect the **fixed** primary margin `Δ R²_genetic = 0.005`, and what is the procedure's false-positive rate at `Δ = 0`?

The margin is fixed a priori. This gate sizes `R`; it never selects the margin from results (that would be silent goal drift).

## 🎯 Procedure fidelity

The simulation mirrors Gate 0A exactly: the inferential unit is the simulation replicate; the statistic is the paired difference `Δ_r = R²_genetic(A) − R²_genetic(B)` per replicate; uncertainty is a percentile bootstrap CI over the `R` replicates; a regime is declared to have positive headroom iff the CI lower bound exceeds 0 (one-sided).

## 🧮 Model and the shift/scale identity

Replicate differences are modelled as `Δ_r ~ Normal(Δ_true, σ²)`, where `σ` is the replicate-to-replicate SD of the paired `R²_genetic` difference (trait-architecture plus estimation noise). Because for `X_r = Δ + σ·ε_r` the bootstrap CI lower bound satisfies `L(Δ,σ) = Δ + σ·L₀` with `L₀` computed on the centered unit sample, the decision `L(Δ,σ) > 0` is exactly `L₀ > −Δ/σ`. So `L₀` is simulated once per `R` and every `(Δ, σ)` cell follows in closed form.

## 📊 Frozen grids

- `Δ ∈ {0, 0.002, 0.005, 0.010}` — FPR, sensitivity, **primary power**, strong-effect power.
- `R ∈ {20, 50, 100, 200}`.
- `σ ∈ {0.005, 0.01, 0.02, 0.03, 0.05}` — a **planning axis**; the real value is pinned by a real-panel pilot.
- Monte Carlo: `M = 3000` outer draws, `B = 1500` bootstrap resamples, seed `20260904`.

## ✅ Criterion

Fix the smallest `R` with `P(CI lower > 0 | Δ = 0.005) ≥ 0.90` at the pinned `σ`, provided `FPR@0` stays near the one-sided nominal `~0.025`. If no gridded `R` suffices, increase `R` — never lower the margin.

## 📈 Result (this run; planning proxy)

`results/detectability_table.md` / `.json`. Minimum `R` for ≥0.90 power at `Δ=0.005`, by `σ`:

| σ | 0.005 | 0.01 | 0.02 | 0.03 | 0.05 |
| --- | ---: | ---: | ---: | ---: | ---: |
| min R | 20 | 50 | 200 | > 200 | > 200 |

`FPR@0 ≈ 0.03` across `R` (one-sided, as expected). **Feasibility warning:** if the paired-`R²_genetic` replicate SD is `σ ≳ 0.03`, even 200 replicates cannot reach 90% power at 0.005 — Gate 0A would then need either more replicates or a reconsidered `N`/`h²`. This is exactly why the real-panel `σ` pilot is a pre-run gate.

## 🔒 Firewall

`PLANNING_PROXY`. Running is permitted (no real data, no biological claim), but `R` is fixed for the real Gate 0A run only after a real-panel pilot pins `σ`. Results are never formal Gate 0A evidence.
