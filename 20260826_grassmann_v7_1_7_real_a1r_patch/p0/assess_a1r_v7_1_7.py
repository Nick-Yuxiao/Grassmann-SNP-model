from __future__ import annotations

import argparse
import json
from pathlib import Path


MASKS = {"ld_block_0p90", "within_chrom_longrange_0p90"}
FAIRNESS = {"matched_parameter", "matched_compute"}


def main() -> None:
    parser = argparse.ArgumentParser(description="Apply the frozen v7.1.7 three-state A1-R decision")
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    data = json.loads(args.summary.read_text(encoding="utf-8"))
    if args.output.exists():
        raise SystemExit(f"refusing to overwrite output: {args.output}")
    if data.get("protocol_version") != "v7.1.7" or data.get("hgdp_accessed") is not False:
        verdict = "INVALID_PROTOCOL_VIOLATION"
        reasons = ["protocol version mismatch or HGDP was accessed"]
    else:
        cells = data.get("primary_inference_cells", [])
        keyed = {(c.get("mask"), c.get("fairness")): c for c in cells}
        required = {(m, f) for m in MASKS for f in FAIRNESS}
        complete = set(keyed) == required and data.get("primary_completed_runs") == 120
        converged = complete and data.get("primary_converged_runs") == 120
        delta_min = 0.010
        go = complete and all(
            keyed[k]["mask_seed_ci_lcb"] > delta_min and keyed[k]["population_ci_lcb"] > delta_min
            for k in required
        )
        excludes_meaningful_advantage = complete and all(
            keyed[k]["mask_seed_ci_ucb"] <= delta_min and keyed[k]["population_ci_ucb"] <= delta_min
            for k in required
        )
        size_rows = data.get("sample_size_diagnostic", [])
        size_keyed = {r.get("mask"): r for r in size_rows}
        size_complete = set(size_keyed) == MASKS and data.get("diagnostic_completed_runs") == 24
        trend_masks = []
        if size_complete:
            for mask in MASKS:
                r = size_keyed[mask]
                d25, d50, d100 = float(r["delta_25"]), float(r["delta_50"]), float(r["delta_100"])
                if d50 > d25 and d100 > d50 and d100 - d25 >= 0.002:
                    trend_masks.append(mask)
        size_trend = bool(trend_masks)
        reasons = []
        if go:
            verdict = "A1R_PRELIMINARY_GO"
            reasons.append("all simultaneous lower bounds exceed +0.010")
        elif converged and excludes_meaningful_advantage and size_complete and not size_trend:
            verdict = "NO_GO_EQUIVALENT_OR_WORSE"
            reasons.append("all primary curves converged and all simultaneous upper bounds are <= +0.010")
            reasons.append("complete nested-size diagnostic did not trigger sample limitation")
        else:
            verdict = "INCONCLUSIVE_SAMPLE_LIMITED_HAPNEST"
            if not complete:
                reasons.append("primary grid incomplete")
            if not converged:
                reasons.append("all-primary convergence not established")
            if not excludes_meaningful_advantage:
                reasons.append("precision insufficient to exclude a meaningful advantage")
            if not size_complete:
                reasons.append("nested-size diagnostic incomplete")
            if size_trend:
                reasons.append("positive donor-size trend: " + ",".join(sorted(trend_masks)))
    out = {
        "schema_version": "1.0",
        "protocol_version": "v7.1.7",
        "verdict": verdict,
        "project_termination": False if verdict == "INCONCLUSIVE_SAMPLE_LIMITED_HAPNEST" else None,
        "next_branch": "HAPNEST" if verdict == "INCONCLUSIVE_SAMPLE_LIMITED_HAPNEST" else None,
        "reasons": reasons,
    }
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(out, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

