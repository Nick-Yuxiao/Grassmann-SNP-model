#!/usr/bin/env python3
"""Run the Gate 0A sigma pilot on the bound frozen contract (server-side).

Refuses to run unless: BINDING.json status == BOUND, an EXPORT_MANIFEST.json is
present and consistent with BINDING and the exported arrays, and the panel matches
the contract. Runs on the TRAINING population only (val 249 never touched).
Estimates SD(Delta R2_genetic) per (regime, k), maps each to a required R, and
reports R_formal = max over cells. Variance calibration only.

Usage:
    python scripts/run_sigma_pilot.py \
        --binding /path/to/rc2/BINDING.json \
        --data-dir /path/to/exported_bound_panel \
        --out results/sigma_pilot_result.json
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

PKG = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PKG / "src"))

import numpy as np  # noqa: E402
from gate0a import estimator, load_frozen_panel  # noqa: E402
from gate0a.sigma_map import required_replicates  # noqa: E402

CFG = json.loads((PKG / "config" / "SIGMA_PILOT_CONFIG.json").read_text())


def _sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def _require_bound(binding_path):
    b = json.loads(Path(binding_path).read_text())
    if b.get("status") != "BOUND":
        sys.exit(f"RUN_BLOCKED: BINDING.json status={b.get('status')} (need BOUND).")
    return b


def _verify_export(binding, data_dir):
    em_path = Path(data_dir) / "EXPORT_MANIFEST.json"
    if not em_path.is_file():
        sys.exit("RUN_BLOCKED: EXPORT_MANIFEST.json missing; run write_export_manifest.py first.")
    em = json.loads(em_path.read_text())
    if em.get("source_panel_manifest_sha256") != binding["panel_manifest"]["sha256"]:
        sys.exit("RUN_BLOCKED: export panel SHA does not match BINDING panel SHA.")
    if em.get("source_block_version_sha256") != binding["block_version"]["sha256"]:
        sys.exit("RUN_BLOCKED: export block SHA does not match BINDING block SHA.")
    for f, h in em["arrays"].items():
        if _sha(Path(data_dir) / f) != h:
            sys.exit(f"RUN_BLOCKED: exported {f} hash differs from EXPORT_MANIFEST.")
    return em


def _pick_blocks(blocks, regime, max_blocks):
    uniq = list(np.unique(blocks))
    sizes = {b: int(np.sum(blocks == b)) for b in uniq}
    big = [b for b in uniq if sizes[b] >= 6]
    if not big:
        sys.exit("RUN_BLOCKED: no block has >=6 sites; cannot form directions")
    sig = [big[0], big[-1]] if regime == "between_block_interaction" else [big[0]]
    if max_blocks is None:
        eval_ids = list(big)
    else:
        ctx = [b for b in big if b not in sig]
        eval_ids = list(dict.fromkeys(sig + ctx))[:max_blocks]
    return sig, eval_ids


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--binding", required=True)
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--out", default=str(PKG / "results" / "sigma_pilot_result.json"))
    args = ap.parse_args()

    binding = _require_bound(args.binding)
    _verify_export(binding, args.data_dir)
    panel = load_frozen_panel(args.data_dir)

    # TRAINING POPULATION ONLY -- the 249 validation rows are never touched
    tr_idx = panel["train_idx"]
    G = np.asarray(panel["G"])[tr_idx]
    blocks = panel["blocks"]

    regimes = list(CFG["representative_regimes"].keys())
    ks = CFG["budget_points_per_block_k"]
    k_max = max(ks)
    R_pilot = CFG["R_pilot"]
    folds_n = CFG["outer_cv_folds"]
    inner_folds = CFG["inner_cv_folds"]
    lambdas = CFG["ridge_lambda_grid"]
    n_dir = CFG["n_causal_directions"]
    poly_frac = CFG["polygenic_background_fraction"]
    max_blocks = CFG["eval_block_subset"]["max_blocks"]
    seed_base = CFG["seed_base"]
    R_grid = CFG["replicate_count_rule"]["R_candidate_grid"]

    # per-regime block choices; union of all blocks we must precompute reps for
    regime_blocks = {r: _pick_blocks(blocks, r, max_blocks) for r in regimes}
    needed = sorted(set().union(*[set(sig) | set(ev) for sig, ev in regime_blocks.values()]))

    # folds drawn ONCE over the training population, shared everywhere
    folds = estimator.build_folds(G.shape[0], folds_n, seed=seed_base)
    # arm reps computed ONCE per (arm, fold, block) at k_max; sliced for k=8
    cache = estimator.precompute_reps(G, blocks, folds, k_max, needed)

    cells, required = [], []
    for ri, regime in enumerate(regimes):
        sig, eval_ids = regime_blocks[regime]
        for k in ks:
            res = estimator.sigma_for_cell(
                cache, G, blocks, folds, regime, sig, eval_ids, k,
                R_pilot, seed=seed_base + 1000 * ri + k,
                lambdas=lambdas, inner_folds=inner_folds, n_dir=n_dir, poly_frac=poly_frac)
            R_needed = required_replicates(max(res["sigma_hat"], 1e-6), R_grid=R_grid)
            cells.append({
                "regime": regime, "k": k, "sig_blocks": [int(x) for x in sig],
                "n_eval_blocks": len(eval_ids),
                "mean_delta_NON_EVIDENCE_DIAGNOSTIC_ONLY":
                    res["mean_delta_NON_EVIDENCE_DIAGNOSTIC_ONLY"],
                "sigma_hat": res["sigma_hat"], "R_needed": R_needed,
            })
            required.append(R_needed)

    R_formal = (max(r for r in required) if required and all(r is not None for r in required)
                else None)
    payload = {
        "classification": "VARIANCE_CALIBRATION_ONLY",
        "estimand": CFG["estimand"],
        "trained_on": "donor-train only (2247); validation 249 untouched",
        "outer_folds_drawn_once": True,
        "R_pilot": R_pilot,
        "cells": cells,
        "R_formal": R_formal,
        "R_formal_rule": "max over primary-relevant cells; None => a cell is underpowered "
                         "on the candidate grid -> increase R, never lower margin",
        "reminder": "variance calibration only; decides only R and compute feasibility",
    }
    Path(args.out).write_text(json.dumps(payload, indent=2))
    print("wrote", args.out)
    for c in cells:
        print(f"  {c['regime']:<28} k={c['k']:<3} sigma_hat={c['sigma_hat']:.4f} "
              f"R_needed={c['R_needed']}")
    print("R_formal =", R_formal)


if __name__ == "__main__":
    main()
