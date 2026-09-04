#!/usr/bin/env python3
"""Run the Gate 0A detectability planning simulation and write the frozen-structure
power/FPR tables. Planning proxy only: no real data, no biological claim.

Usage:
    python scripts/run_detectability.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

PKG = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PKG / "src"))

from detectability import power_table, min_replicates  # noqa: E402

CFG = json.loads((PKG / "config" / "DETECTABILITY_CONFIG.json").read_text())
G = CFG["grids"]
MC = CFG["monte_carlo"]
DELTAS = G["delta_grid"]
SIGMAS = G["sigma_grid"]
RS = G["replicate_grid_R"]
TARGET_DELTA = CFG["margins"]["primary_practical_margin"]
TARGET_POWER = 0.90


def main() -> None:
    out_json = PKG / "results" / "detectability_table.json"
    out_md = PKG / "results" / "detectability_table.md"

    table, _ = power_table(DELTAS, SIGMAS, RS, M=MC["outer_draws_M"],
                           B=MC["bootstrap_resamples_B"], ci=0.95, seed=MC["seed"])

    # structured records
    records = []
    summary = {}
    for sigma in SIGMAS:
        for R in RS:
            rec = {
                "sigma": sigma,
                "R": R,
                "FPR@0": table[(0.0, sigma, R)],
                "power@0.002": table[(0.002, sigma, R)],
                "power@0.005": table[(0.005, sigma, R)],
                "power@0.010": table[(0.010, sigma, R)],
            }
            records.append(rec)
        summary[str(sigma)] = min_replicates(table, sigma, RS, TARGET_DELTA, TARGET_POWER)

    payload = {
        "classification": "PLANNING_PROXY",
        "target": f"P(CI lower>0 | Delta={TARGET_DELTA}) >= {TARGET_POWER}",
        "monte_carlo": MC,
        "records": records,
        "min_R_for_0.90_power_at_0.005_per_sigma": summary,
        "reminder": "sigma is a planning axis; a real-panel pilot must pin sigma before R is fixed for Gate 0A",
    }
    out_json.write_text(json.dumps(payload, indent=2))

    lines = [
        "# Gate 0A detectability — planning proxy",
        "",
        "_No real data. `sigma` = replicate-to-replicate SD of the paired R2_genetic "
        "difference; it is a planning axis a real-panel pilot must pin. Decision is "
        "one-sided (CI lower bound > 0), so FPR@0 for a calibrated 95% CI is ~0.025._",
        "",
        f"Target: **P(CI lower>0 | Δ={TARGET_DELTA}) ≥ {TARGET_POWER}**  "
        f"(M={MC['outer_draws_M']}, B={MC['bootstrap_resamples_B']}, seed={MC['seed']})",
        "",
    ]
    for sigma in SIGMAS:
        lines.append(f"## sigma = {sigma}")
        lines.append("")
        lines.append("| Replicates R | FPR @ 0 | Power @ .002 | Power @ .005 | Power @ .010 |")
        lines.append("| ---: | ---: | ---: | ---: | ---: |")
        for R in RS:
            r = next(x for x in records if x["sigma"] == sigma and x["R"] == R)
            lines.append(
                f"| {R} | {r['FPR@0']:.3f} | {r['power@0.002']:.3f} | "
                f"{r['power@0.005']:.3f} | {r['power@0.010']:.3f} |"
            )
        mr = summary[str(sigma)]
        lines.append("")
        lines.append(f"**Min R for ≥0.90 power at Δ=0.005:** "
                     f"{mr if mr is not None else '> 200 (none in grid; increase R)'}")
        lines.append("")
    out_md.write_text("\n".join(lines))

    print("wrote:", out_json.name, out_md.name)
    print("min_R_for_0.90_power_at_0.005_per_sigma:", summary)


if __name__ == "__main__":
    main()
