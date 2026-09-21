#!/usr/bin/env python3
"""LEVEL 0 CALIBRATION -- estimator validation on known truth.  PLANNING PROXY.

This is not the reality gate.  It is the thing that has to pass before the
reality gate is worth opening, in the same sense as this project's existing
detectability package: it fixes how many pairs are needed and shows the decision
procedure does what it claims, using synthetic worlds where the answer is known.

Four worlds, and the fourth carries the weight:

  S  subspace-conserved   E_B = E_A M, M in GL(2)      the hypothesis
  C  fully conserved      E_B = E_A                    coordinates survive too
  N  unrelated            independent draw             nothing survives
  D  decorative           per-column rotations matched to S, joint part scrambled

Worlds S and D have *identical single-SNP marginals*.  Every per-SNP direction
change and every per-SNP magnitude is the same.  They differ only in whether the
two columns moved together.  If the deflation test cannot separate them, no
pair-level claim -- and therefore no Grassmann claim -- can be supported by this
estimand, because two individually stable non-collinear vectors span a stable
plane as a matter of arithmetic rather than biology.

Reported per world: I_sub, I_coord, their gap with a paired CI, and the excess
joint conservation test.  Plus a detectability sweep over the number of pairs.
"""

from __future__ import annotations

import argparse
import csv
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from subspace_conservation.estimators import (  # noqa: E402
    conservation_fraction,
    excess_joint_conservation,
    level0_verdict,
    observed_distances,
    paired_bootstrap_ci,
    permutation_reference,
    standardised_excess,
)
from subspace_conservation.nulls import parametric_null  # noqa: E402
from subspace_conservation.simulate import (  # noqa: E402
    TruthSpec,
    apply_distortion,
    draw_truth,
    make_decorative_matched,
    observe,
)

STATUS = "PLANNING_PROXY_NOT_REALITY_GATE_EVIDENCE"

EXPECTED = {
    "S": "PROCEED_TO_LEVEL_1",
    "C": "NO_CONTEXT_DISTORTION_ROUTE_UNMOTIVATED",
    "N": "NO_SUBSPACE_CONSERVATION_STOP",
    "D": "CONSERVED_BUT_DECORATIVE_STOP",
}


def build_world(kind, n_pairs, spec, distortion, rng):
    A = draw_truth(n_pairs, spec, rng)
    if kind == "S":
        B = apply_distortion(A, distortion, rng)
    elif kind == "C":
        B = A.copy()
    elif kind == "N":
        B = draw_truth(n_pairs, spec, rng)
    elif kind == "D":
        B = make_decorative_matched(A, apply_distortion(A, distortion, rng), rng)
    else:
        raise ValueError(kind)
    EA, SA = observe(A, spec.noise, rng)
    EB, SB = observe(B, spec.noise, rng)
    return EA, SA, EB, SB


def analyse(EA, SA, EB, SB, rng, n_null, n_surrogate, n_boot, metric="chordal"):
    null = parametric_null([EA, EB], [SA, SB], n_null, rng, metric=metric)
    d_sub, d_coord = observed_distances(EA, EB, metric=metric)

    res_sub = standardised_excess(d_sub, null.sub, "subspace")
    res_coord = standardised_excess(d_coord, null.coord, "coordinate")

    # Both metrics are put on one absolute [0, 1] scale by dividing out their
    # own noise floor and their own unrelated-pair ceiling. Raw distances on two
    # different spaces cannot be compared; these can.
    perm_sub, perm_coord = permutation_reference(EA, EB, rng, metric=metric)
    c_sub = conservation_fraction(d_sub, null.sub, perm_sub)
    c_coord = conservation_fraction(d_coord, null.coord, perm_coord)
    gap, gap_lo, gap_hi = paired_bootstrap_ci(c_sub, c_coord, n_boot, rng)

    joint = excess_joint_conservation(EA, EB, n_surrogate, rng, metric=metric, n_boot=n_boot)
    verdict = level0_verdict(float(np.nanmean(c_sub)), float(np.nanmean(c_coord)), gap_lo, joint)
    return {
        "subspace": res_sub.to_dict() | {"conservation_fraction": float(np.nanmean(c_sub))},
        "coordinate": res_coord.to_dict() | {"conservation_fraction": float(np.nanmean(c_coord))},
        "conservation_gap": {"estimate": gap, "ci_low": gap_lo, "ci_high": gap_hi,
                             "note": "mean C_sub - mean C_coord; >0 means the plane survives the "
                                     "context change better than the coordinates do"},
        "excess_joint_conservation": joint,
        "verdict": verdict,
    }


def detectability(kind, spec, distortion, n_grid, n_rep, rng, n_null, n_surrogate, n_boot):
    rows = []
    for n_pairs in n_grid:
        hits = 0
        for _ in range(n_rep):
            EA, SA, EB, SB = build_world(kind, n_pairs, spec, distortion, rng)
            out = analyse(EA, SA, EB, SB, rng, n_null, n_surrogate, n_boot)
            hits += int(out["verdict"]["verdict"] == EXPECTED[kind])
        rows.append(
            {
                "world": kind,
                "n_pairs": n_pairs,
                "replicates": n_rep,
                "correct_verdict_rate": hits / n_rep,
                "expected_verdict": EXPECTED[kind],
            }
        )
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seed", type=int, default=20260921)
    ap.add_argument("--pairs", type=int, default=300)
    ap.add_argument("--r", type=int, default=8)
    ap.add_argument("--sigma1", type=float, default=8.0)
    ap.add_argument("--sigma2", type=float, default=4.0)
    ap.add_argument("--noise", type=float, default=1.0)
    ap.add_argument("--distortion", default="gl_8", help="GL(2) level defining worlds S and D")
    ap.add_argument("--null-draws", type=int, default=60)
    ap.add_argument("--surrogate-draws", type=int, default=60)
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--detectability-grid", type=int, nargs="*", default=[50, 100, 200, 400])
    ap.add_argument("--detectability-reps", type=int, default=10)
    ap.add_argument("--skip-detectability", action="store_true")
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parents[1] / "results")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    spec = TruthSpec(r=args.r, sigma1=args.sigma1, sigma2=args.sigma2, noise=args.noise)

    worlds = {}
    for kind in ("S", "C", "N", "D"):
        EA, SA, EB, SB = build_world(kind, args.pairs, spec, args.distortion, rng)
        worlds[kind] = analyse(
            EA, SA, EB, SB, rng, args.null_draws, args.surrogate_draws, args.boot
        )
        v = worlds[kind]["verdict"]
        print(
            f"[L0-cal] world {kind}: C_sub={v['conservation_sub']:+.2f} "
            f"C_coord={v['conservation_coord']:+.2f} "
            f"joint={worlds[kind]['excess_joint_conservation']['verdict']} -> {v['verdict']}"
        )

    checks = {
        f"world_{k}_recovered": worlds[k]["verdict"]["verdict"] == EXPECTED[k]
        for k in EXPECTED
    }
    checks["S_and_D_separated"] = (
        worlds["S"]["excess_joint_conservation"]["verdict"] == "EXCESS_JOINT_CONSERVATION"
        and worlds["D"]["excess_joint_conservation"]["verdict"] == "NO_EXCESS_JOINT_CONSERVATION"
    )

    detect = []
    if not args.skip_detectability:
        for kind in ("S", "D"):
            detect += detectability(
                kind, spec, args.distortion, args.detectability_grid,
                args.detectability_reps, rng, args.null_draws, args.surrogate_draws,
                min(args.boot, 1000),
            )
            for row in detect[-len(args.detectability_grid):]:
                print(
                    f"[L0-cal] detectability world {kind} n_pairs={row['n_pairs']:>4} "
                    f"correct={row['correct_verdict_rate']:.2f}"
                )

    report = {
        "status": STATUS,
        "level": 0,
        "role": "estimator calibration on synthetic truth; fixes n_pairs for the real run",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "config": vars(args) | {"out": str(args.out)},
        "environment": {"python": platform.python_version(), "numpy": np.__version__},
        "worlds": worlds,
        "checks": checks,
        "gate": "PASS" if all(checks.values()) else "FAIL",
        "detectability": detect,
        "not_authorised": [
            "reading these numbers as evidence about genetics",
            "choosing the real run's margin from these results",
            "skipping the reality gate because the estimator works here",
        ],
    }

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "L0_CALIBRATION.json").write_text(json.dumps(report, indent=2))
    if detect:
        with (args.out / "l0_detectability.csv").open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(detect[0].keys()))
            w.writeheader()
            w.writerows(detect)

    print(f"[L0-cal] {report['gate']}  ->  {args.out / 'L0_CALIBRATION.json'}")
    for k, v in checks.items():
        print(f"  {'PASS' if v else 'FAIL'}  {k}")
    return 0 if report["gate"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
