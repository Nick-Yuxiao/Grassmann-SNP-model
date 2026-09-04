#!/usr/bin/env python3
"""Run the Gate 0A sigma pilot on the bound frozen contract (server-side).

Refuses to run unless BINDING.json status == BOUND. Estimates SD(Delta R2_genetic)
per (regime, k) cell, maps each to a required R, and reports R_formal = max over
cells. Estimates variance ONLY.

Usage:
    python scripts/run_sigma_pilot.py \
        --binding /path/to/rc2/BINDING.json \
        --data-dir /path/to/exported_frozen_panel_npy \
        --out results/sigma_pilot_result.json
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

PKG = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PKG / "src"))

import numpy as np  # noqa: E402
from gate0a import estimator, load_frozen_panel  # noqa: E402
from gate0a.sigma_map import required_replicates  # noqa: E402

CFG = json.loads((PKG / "config" / "SIGMA_PILOT_CONFIG.json").read_text())


def _require_bound(binding_path):
    b = json.loads(Path(binding_path).read_text())
    if b.get("status") != "BOUND":
        sys.exit(f"RUN_BLOCKED: BINDING.json status={b.get('status')} (need BOUND). "
                 "Bind the panel/block hashes and verify preprocessing first.")
    return b


def _pick_blocks(blocks, regime, max_blocks):
    uniq = list(np.unique(blocks))
    sizes = {b: int(np.sum(blocks == b)) for b in uniq}
    big = [b for b in uniq if sizes[b] >= 6]
    if not big:
        sys.exit("RUN_BLOCKED: no block has >=6 sites; cannot form directions")
    if regime == "between_block_interaction":
        sig = [big[0], big[-1]]           # two distant blocks
    else:
        sig = [big[0]]
    if max_blocks is None:
        eval_ids = uniq
    else:
        # signal blocks + context, capped
        ctx = [b for b in big if b not in sig]
        eval_ids = list(dict.fromkeys(sig + ctx))[:max_blocks]
    return sig, eval_ids


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--binding", required=True)
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--out", default=str(PKG / "results" / "sigma_pilot_result.json"))
    args = ap.parse_args()

    _require_bound(args.binding)
    panel = load_frozen_panel(args.data_dir)
    G, blocks = panel["G"], panel["blocks"]

    regimes = list(CFG["representative_regimes"].keys())
    ks = CFG["budget_points_per_block_k"]
    R_pilot = CFG["R_pilot"]
    folds = CFG["outer_cv_folds"]
    max_blocks = CFG["eval_block_subset"]["max_blocks"]
    seed_base = CFG["seed_base"]
    R_grid = CFG["replicate_count_rule"]["R_candidate_grid"]

    cells, required = [], []
    for ri, regime in enumerate(regimes):
        sig, eval_ids = _pick_blocks(blocks, regime, max_blocks)
        for ki, k in enumerate(ks):
            seed = seed_base + 1000 * ri + k
            res = estimator.sigma_for_cell(G, blocks, regime, sig, k, R_pilot,
                                           seed, folds=folds, eval_block_ids=eval_ids)
            R_needed = required_replicates(max(res["sigma_hat"], 1e-6), R_grid=R_grid)
            cells.append({
                "regime": regime, "k": k,
                "sig_blocks": [int(x) for x in sig],
                "n_eval_blocks": len(eval_ids),
                "mean_delta": res["mean_delta"],
                "sigma_hat": res["sigma_hat"],
                "R_needed": R_needed,
            })
            required.append(R_needed)

    finite = [r for r in required if r is not None]
    R_formal = max(finite) if finite and len(finite) == len(required) else None
    payload = {
        "classification": "VARIANCE_CALIBRATION_ONLY",
        "estimand": CFG["estimand"],
        "R_pilot": R_pilot,
        "cells": cells,
        "R_formal": R_formal,
        "R_formal_rule": "max over primary-relevant cells; None means at least one cell "
                         "is underpowered on the candidate grid -> increase R, never lower margin",
        "reminder": "variance calibration only; does not decide winner/margin/DGP/k/arm/threshold/primary",
    }
    Path(args.out).write_text(json.dumps(payload, indent=2))
    print("wrote", args.out)
    for c in cells:
        print(f"  {c['regime']:<28} k={c['k']:<3} sigma_hat={c['sigma_hat']:.4f} "
              f"R_needed={c['R_needed']}")
    print("R_formal =", R_formal)


if __name__ == "__main__":
    main()
