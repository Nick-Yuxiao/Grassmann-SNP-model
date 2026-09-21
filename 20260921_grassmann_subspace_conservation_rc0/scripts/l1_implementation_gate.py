#!/usr/bin/env python3
"""LEVEL 1 -- IMPLEMENTATION GATE.  NOT SCIENTIFIC EVIDENCE.

Under a GL(2) world where the phenotype is a function of the column space only,
the Grassmann arm transfers and the coordinate arm does not.  That is three
lines of algebra, not an experiment, and running it proves nothing about
genetics.  What it can establish -- and what this script is for -- is whether
the *implementation* is invariant where the algebra says it should be, and where
it stops being invariant once estimation noise is in the picture.

Four things get measured:

1. GL(2) invariance of the projector, against the condition number of M.
2. The O(2) / GL(2) separation.  Under right-O(2) the second moment E E.T is
   exactly invariant and strictly richer than the projector, so an O(2)-only
   world hands the spectrum arm a free win over Grassmann. This is the flaw
   that reshaped the design; it is measured here rather than asserted.
3. The identifiability curve: how fast the estimated subspace decays as
   sigma2/noise falls. This sets the y axis of the Level 2 surface.
4. The per-pair GL(2) alignment identity: the best per-pair GL(2) realignment
   residual *is* the chordal subspace distance, so any baseline fitting M on
   the test pair's own target data is the Grassmann arm plus leakage.
"""

from __future__ import annotations

import argparse
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from subspace_conservation.geometry import (  # noqa: E402
    gl2_alignment_residual,
    projector,
    random_gl2,
    random_orthogonal,
    subspace_distance,
)
from subspace_conservation.simulate import TruthSpec, draw_truth, observe  # noqa: E402

STATUS = "IMPLEMENTATION_GATE_NOT_SCIENTIFIC_EVIDENCE"


def invariance_vs_condition(rng, r, n_trials, conds):
    rows = []
    for cond in conds:
        devs = []
        for _ in range(n_trials):
            E = rng.standard_normal((r, 2))
            M = random_gl2(cond, rng)
            devs.append(subspace_distance(E, E @ M))
        rows.append(
            {
                "cond_M": float(cond),
                "max_subspace_distance": float(np.max(devs)),
                "median_subspace_distance": float(np.median(devs)),
                "invariant_to_1e-9": bool(np.max(devs) < 1e-9),
            }
        )
    return rows


def o2_versus_gl2(rng, r, n_trials):
    """Is the second moment invariant under right-O(2) but not right-GL(2)?"""
    o2_dev, gl2_dev = [], []
    for _ in range(n_trials):
        E = rng.standard_normal((r, 2))
        G = E @ E.T
        R = random_orthogonal(2, rng)
        M = random_gl2(8.0, rng)
        scale = np.linalg.norm(G)
        o2_dev.append(np.linalg.norm((E @ R) @ (E @ R).T - G) / scale)
        gl2_dev.append(np.linalg.norm((E @ M) @ (E @ M).T - G) / scale)
    return {
        "second_moment_under_O2_max_rel_dev": float(np.max(o2_dev)),
        "second_moment_under_GL2_median_rel_dev": float(np.median(gl2_dev)),
        "verdict": (
            "O2_ONLY_WORLD_IS_DEGENERATE"
            if np.max(o2_dev) < 1e-10 < np.median(gl2_dev)
            else "UNEXPECTED"
        ),
        "why_it_matters": (
            "E E.T survives right-O(2) untouched and carries the spectrum the "
            "projector throws away. A rotation-only positive control therefore "
            "rewards the spectrum arm, not the Grassmann arm. Only GL(2) makes "
            "the column space the maximal invariant."
        ),
    }


def identifiability_curve(rng, r, n_trials, ratios, noise=1.0, sigma1=8.0):
    rows = []
    for ratio in ratios:
        spec = TruthSpec(r=r, sigma1=sigma1, sigma2=ratio * noise, noise=noise)
        truth = draw_truth(n_trials, spec, rng)
        hat, _ = observe(truth, noise, rng)
        d_sub = [subspace_distance(truth[p], hat[p]) for p in range(n_trials)]
        # A uniformly random plane is the "no information left" reference.
        rand = [
            subspace_distance(truth[p], rng.standard_normal((r, 2))) for p in range(n_trials)
        ]
        rows.append(
            {
                "sigma2_over_noise": float(ratio),
                "sigma1_over_sigma2": float(sigma1 / (ratio * noise)),
                "mean_subspace_error": float(np.mean(d_sub)),
                "mean_random_plane_distance": float(np.mean(rand)),
                "fraction_of_random": float(np.mean(d_sub) / np.mean(rand)),
            }
        )
    return rows


def alignment_identity(rng, r, n_trials):
    gaps = []
    for _ in range(n_trials):
        A = rng.standard_normal((r, 2))
        B = rng.standard_normal((r, 2))
        gaps.append(abs(gl2_alignment_residual(A, B) - subspace_distance(A, B, "chordal")))
    return {
        "max_abs_gap": float(np.max(gaps)),
        "identity_holds_to_1e-10": bool(np.max(gaps) < 1e-10),
        "consequence": (
            "min_M ||Qb - Qa M||_F equals the chordal subspace distance, so a "
            "per-pair GL(2) baseline fitted on the test pair's own target data "
            "is not a baseline -- it is the Grassmann arm with the adaptation "
            "cost hidden. Transports must be fitted on calibration pairs."
        ),
    }


def rank_deficiency_behaviour(rng, r, n_trials):
    """What the projector does as the pair becomes collinear."""
    rows = []
    for eps in (1e-1, 1e-3, 1e-6, 1e-9, 0.0):
        finite, ranks = 0, []
        for _ in range(n_trials):
            v = rng.standard_normal(r)
            w = rng.standard_normal(r)
            E = np.column_stack([v, v + eps * w])
            P = projector(E)
            finite += int(np.all(np.isfinite(P)))
            ranks.append(int(round(np.trace(P))))
        rows.append(
            {
                "collinearity_eps": float(eps),
                "all_finite": bool(finite == n_trials),
                "median_effective_rank": float(np.median(ranks)),
            }
        )
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seed", type=int, default=20260921)
    ap.add_argument("--r", type=int, default=8, help="number of traits")
    ap.add_argument("--trials", type=int, default=400)
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parents[1] / "results")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    report = {
        "status": STATUS,
        "level": 1,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "seed": args.seed,
        "r_traits": args.r,
        "trials": args.trials,
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "platform": platform.platform(),
        },
        "gl2_invariance": invariance_vs_condition(
            rng, args.r, args.trials, [1.0, 2.0, 10.0, 1e2, 1e4, 1e6]
        ),
        "o2_versus_gl2": o2_versus_gl2(rng, args.r, args.trials),
        "identifiability": identifiability_curve(
            rng, args.r, args.trials, [0.25, 0.5, 1.0, 2.0, 4.0, 8.0]
        ),
        "gl2_alignment_identity": alignment_identity(rng, args.r, args.trials),
        "rank_deficiency": rank_deficiency_behaviour(rng, args.r, min(args.trials, 200)),
    }

    checks = {
        "projector_gl2_invariant_to_cond_1e4": all(
            row["invariant_to_1e-9"] for row in report["gl2_invariance"] if row["cond_M"] <= 1e4
        ),
        "o2_world_is_degenerate": report["o2_versus_gl2"]["verdict"] == "O2_ONLY_WORLD_IS_DEGENERATE",
        "alignment_identity_holds": report["gl2_alignment_identity"]["identity_holds_to_1e-10"],
        "no_nan_under_collinearity": all(row["all_finite"] for row in report["rank_deficiency"]),
    }
    report["checks"] = checks
    report["gate"] = "PASS" if all(checks.values()) else "FAIL"

    args.out.mkdir(parents=True, exist_ok=True)
    path = args.out / "L1_IMPLEMENTATION_GATE.json"
    path.write_text(json.dumps(report, indent=2))
    print(f"[L1] {report['gate']}  ->  {path}")
    for k, v in checks.items():
        print(f"  {'PASS' if v else 'FAIL'}  {k}")
    return 0 if report["gate"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
