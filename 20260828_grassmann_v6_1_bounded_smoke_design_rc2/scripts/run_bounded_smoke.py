from __future__ import annotations

import argparse
import hashlib
import json
import os
import resource
import sys
import time
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json_atomic(path: Path, value: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n"
    )
    temporary.replace(path)


def load_approval(path: Path) -> dict[str, object]:
    approval = json.loads(path.read_text(encoding="utf-8"))
    expected_manifest_hash = sha256_file(ROOT / "MANIFEST.sha256")
    checks = {
        "authorization": approval.get("authorization") == "APPROVE_BOUNDED_SMOKE_RC2_RUN",
        "approved_by_role": approval.get("approved_by_role") == "project_owner",
        "package_manifest_sha256": approval.get("package_manifest_sha256") == expected_manifest_hash,
        "approved_at_present": bool(approval.get("approved_at")),
    }
    if not all(checks.values()):
        raise PermissionError(f"invalid or mismatched approval record: {checks}")
    return approval


def add_parent_sources(parent: Path) -> None:
    r1 = parent.parent / json.loads((parent / "PARENT_EVIDENCE.json").read_text(encoding="utf-8"))["code_parent"]
    gc_screen = r1.parent / json.loads((r1 / "PARENT_EVIDENCE.json").read_text(encoding="utf-8"))["code_parent"]
    sys.path[:0] = [str(parent / "src"), str(r1 / "src"), str(gc_screen / "src")]


def rotate_family(family, strength, SharedFamily):
    gamma = np.zeros_like(family.Gamma[0])
    gamma[4, 0] = strength
    gamma[5, 1] = strength
    y = family.y + (family.g[:, None] * family.regions[0]) @ gamma
    gammas = list(family.Gamma)
    gammas[0] = gamma
    return SharedFamily(
        family.subject_ids,
        y,
        family.g,
        family.covariates,
        family.regions,
        family.B,
        tuple(gammas),
        family.seed,
    )


def scale_family_signal(family, scale, SharedFamily):
    if scale < 1:
        raise ValueError("smoke signal scale must be at least one")
    extra = np.zeros_like(family.y)
    for x, b, gamma in zip(family.regions, family.B, family.Gamma):
        extra += (scale - 1.0) * (x @ b + (family.g[:, None] * x) @ gamma)
    return SharedFamily(
        family.subject_ids,
        family.y + extra,
        family.g,
        family.covariates,
        family.regions,
        tuple(scale * b for b in family.B),
        tuple(scale * gamma for gamma in family.Gamma),
        family.seed,
    )


def p_on_grid(value: float, resamples: int) -> bool:
    return 1 <= round(value * (resamples + 1)) <= resamples + 1 and np.isclose(
        value * (resamples + 1), round(value * (resamples + 1)), atol=1e-12
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--approval-record", required=True)
    parser.add_argument("--parent-package")
    parser.add_argument("--r1-5-run")
    args = parser.parse_args()

    if os.environ.get("CUDA_VISIBLE_DEVICES") not in (None, "", "-1"):
        raise RuntimeError("bounded smoke is CPU-only")
    for variable in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        if os.environ.get(variable, "1") not in ("", "1"):
            raise RuntimeError(f"{variable} must be 1")

    approval = load_approval(Path(args.approval_record).resolve())
    config = json.loads((ROOT / "config" / "BOUNDED_SMOKE_CONFIG.json").read_text(encoding="utf-8"))
    metadata = json.loads((ROOT / "PARENT_EVIDENCE.json").read_text(encoding="utf-8"))
    parent = (
        Path(args.parent_package).resolve()
        if args.parent_package
        else (ROOT.parent / metadata["code_parent"]).resolve()
    )
    if sha256_file(parent / "MANIFEST.sha256") != metadata["code_parent_manifest_sha256"]:
        raise RuntimeError("R1.5 code-parent manifest mismatch")
    if args.r1_5_run:
        result_manifest = Path(args.r1_5_run).resolve() / "RESULT_MANIFEST.sha256"
        if sha256_file(result_manifest) != metadata["r1_5_result_manifest_sha256"]:
            raise RuntimeError("R1.5 result-manifest mismatch")

    add_parent_sources(parent)
    from grassmann_v6_1.core import conditional_matrices, geometry_scores
    from rank_gate_r1_5 import run_rank_gated_maxT
    from shared_family_r1 import (
        SharedFamily,
        generate_shared_family,
        independent_target_selection,
    )

    output = Path(args.output_dir).resolve()
    if output.exists():
        raise FileExistsError("duplicate output directory refused")
    output.mkdir(parents=True)

    rows: list[dict[str, object]] = []
    failures = 0
    started = time.perf_counter()
    for cell in config["cells"]:
        for replicate in range(config["replicates_per_cell"]):
            data_seed = config["data_seed_base"] + replicate
            bootstrap_seed = config["bootstrap_seed_base"] + replicate
            selection_applied = None
            selection_disjoint = None
            try:
                if cell["kind"] == "selected_null":
                    selection = independent_target_selection(
                        selection_seed=config["selection_seed_base"] + replicate,
                        inference_seed=config["inference_seed_base"] + replicate,
                        n_inference=config["n"],
                    )
                    family = generate_shared_family(
                        seed=config["inference_seed_base"] + 10_000 + replicate,
                        family_size=config["family_size"],
                        conditional_ld=True,
                        heteroskedastic=True,
                        g_override=selection.inference_selected_g,
                        subject_ids_override=selection.inference_subject_ids,
                    )
                    selection_applied = bool(
                        np.array_equal(
                            family.g,
                            selection.inference_target_panel[:, selection.selected_index],
                        )
                    )
                    selection_disjoint = not bool(
                        set(selection.selection_subject_ids) & set(selection.inference_subject_ids)
                    )
                else:
                    amplification = (
                        config["amplification"] if cell["kind"] == "pure_amplification" else 0.0
                    )
                    family = generate_shared_family(
                        seed=data_seed,
                        n=config["n"],
                        family_size=config["family_size"],
                        conditional_ld=cell["conditional_ld"],
                        heteroskedastic=cell["heteroskedastic"],
                        amplification=amplification,
                    )
                family = scale_family_signal(family, config["signal_scale"], SharedFamily)
                if cell["kind"] == "high_snr_rotation":
                    family = rotate_family(family, config["rotation_strength"], SharedFamily)

                result = run_rank_gated_maxT(
                    family,
                    resamples=config["resamples"],
                    seed=bootstrap_seed,
                    rank=config["rank"],
                    ridge_lambda=config["ridge_lambda"],
                    minimum_gap=config["minimum_fitted_rank_gap"],
                )
                truth = max(
                    geometry_scores(conditional_matrices(b, gamma), config["rank"])["direction_score"]
                    for b, gamma in zip(family.B, family.Gamma)
                )
                observed_d29_ok = bool(
                    np.all(result.observed[~result.observed_eligible] == 0)
                    and np.all(result.candidate_p_values[~result.observed_eligible] == 1)
                )
                resampled_d29_ok = bool(
                    np.all(result.resampled[~result.resampled_eligible] == 0)
                )
                row = {
                    "run_id": f"smoke_rc2:{cell['name']}:{replicate}",
                    "cell": cell["name"],
                    "kind": cell["kind"],
                    "replicate": replicate,
                    "data_seed": data_seed,
                    "bootstrap_seed": bootstrap_seed,
                    "family_size": len(result.observed),
                    "subject_rows_aligned": all(len(x) == len(family.subject_ids) for x in family.regions),
                    "unique_subject_ids": len(np.unique(family.subject_ids)) == len(family.subject_ids),
                    "multiplier_fingerprint_count": len(result.multiplier_fingerprints),
                    "p_value": result.family_p_value,
                    "reject_0_05": result.family_p_value <= 0.05,
                    "candidate_p_values": result.candidate_p_values.tolist(),
                    "observed_statistics": result.observed.tolist(),
                    "observed_raw_statistics": result.observed_raw.tolist(),
                    "observed_gaps": result.observed_gaps.tolist(),
                    "observed_eligible": result.observed_eligible.tolist(),
                    "bootstrap_eligible_counts": result.resampled_eligible.sum(axis=1).tolist(),
                    "observed_d29_ok": observed_d29_ok,
                    "resampled_d29_ok": resampled_d29_ok,
                    "truth_direction_max": truth,
                    "selection_applied": selection_applied,
                    "selection_disjoint": selection_disjoint,
                    "success": True,
                    "error": None,
                }
            except Exception as error:  # Preserve failures in the planned denominator.
                failures += 1
                row = {
                    "run_id": f"smoke_rc2:{cell['name']}:{replicate}",
                    "cell": cell["name"],
                    "kind": cell["kind"],
                    "replicate": replicate,
                    "data_seed": data_seed,
                    "bootstrap_seed": bootstrap_seed,
                    "success": False,
                    "error": f"{type(error).__name__}: {error}",
                }
            rows.append(row)

    with (output / "per_run.jsonl").open("x", encoding="utf-8", newline="\n") as stream:
        for row in rows:
            stream.write(json.dumps(row, sort_keys=True) + "\n")

    successful = [row for row in rows if row["success"]]
    nulls = [row for row in successful if row["kind"] in ("null", "selected_null")]
    amplification = [row for row in successful if row["kind"] == "pure_amplification"]
    rotation = [row for row in successful if row["kind"] == "high_snr_rotation"]
    run_ids = [row["run_id"] for row in rows]
    p_values = [row["p_value"] for row in successful]
    candidate_p_values = [value for row in successful for value in row["candidate_p_values"]]

    checks = {
        "all_21_families_present": len(rows) == config["planned_counts"]["independent_families"],
        "unique_run_ids": len(set(run_ids)) == len(run_ids),
        "no_failed_runs": failures == 0,
        "family_dimension_retained": all(row["family_size"] == config["family_size"] for row in successful),
        "subject_rows_aligned": all(row["subject_rows_aligned"] for row in successful),
        "unique_subject_ids": all(row["unique_subject_ids"] for row in successful),
        "one_shared_fingerprint_per_resample": all(
            row["multiplier_fingerprint_count"] == config["resamples"] for row in successful
        ),
        "all_p_values_finite_and_on_grid": all(
            np.isfinite(value) and p_on_grid(value, config["resamples"])
            for value in p_values + candidate_p_values
        ),
        "d29_observed_mapping_exact": all(row["observed_d29_ok"] for row in successful),
        "d29_resampled_mapping_exact": all(row["resampled_d29_ok"] for row in successful),
        "selected_target_applied": all(
            row["selection_applied"] is True for row in successful if row["kind"] == "selected_null"
        ),
        "selection_inference_subject_disjoint": all(
            row["selection_disjoint"] is True for row in successful if row["kind"] == "selected_null"
        ),
        "null_controls_not_minimum_grid_degenerate": bool(nulls)
        and sum(row["p_value"] == 1 / (config["resamples"] + 1) for row in nulls) < len(nulls),
        "pure_amplification_truth_direction_zero": bool(amplification)
        and all(row["truth_direction_max"] < 1e-12 for row in amplification),
        "pure_amplification_rejections_at_most_one": bool(amplification)
        and sum(row["reject_0_05"] for row in amplification) <= 1,
        "rotation_truth_direction_positive": bool(rotation)
        and all(row["truth_direction_max"] > 0.10 for row in rotation),
        "high_snr_rotation_detected_at_least_two": bool(rotation)
        and sum(row["reject_0_05"] for row in rotation) >= 2,
    }

    elapsed_seconds = time.perf_counter() - started
    peak_rss_kib = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    resource_checks = {
        "cpu_time_within_planning_ceiling": elapsed_seconds <= config["resource_ceiling"]["cpu_minutes"] * 60,
        "peak_rss_within_planning_ceiling": peak_rss_kib <= config["resource_ceiling"]["peak_rss_mib"] * 1024,
    }
    status = "BOUNDED_SMOKE_PASS" if all(checks.values()) and all(resource_checks.values()) else "BOUNDED_SMOKE_FAIL"
    gate = {
        "gate": "bounded-shared-family-smoke-rc2",
        "status": status,
        "checks": checks,
        "resource_checks": resource_checks,
        "elapsed_seconds": elapsed_seconds,
        "peak_rss_kib": peak_rss_kib,
        "failed_runs": failures,
        "approval_record_sha256": sha256_file(Path(args.approval_record).resolve()),
        "approval": approval,
        "authorizes": ["prospective_gc_screen_rc2_protocol_draft"] if status.endswith("PASS") else [],
        "does_not_authorize": config["does_not_authorize"],
    }
    write_json_atomic(output / "GATE_BOUNDED_SHARED_FAMILY_SMOKE_RC2.json", gate)
    write_json_atomic(
        output / "ENVIRONMENT.json",
        {"python": sys.version, "numpy": np.__version__, "accelerator_used": False},
    )
    print(json.dumps(gate, indent=2, sort_keys=True))
    return 0 if status.endswith("PASS") else 1


if __name__ == "__main__":
    raise SystemExit(main())
