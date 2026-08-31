from __future__ import annotations

import hashlib
import json
from pathlib import Path, PurePosixPath


ROOT = Path(__file__).resolve().parents[1]


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    config = json.loads((ROOT / "config" / "BOUNDED_SMOKE_CONFIG.json").read_text(encoding="utf-8"))
    parent = json.loads((ROOT / "PARENT_EVIDENCE.json").read_text(encoding="utf-8"))
    decision_lines = (ROOT / "DECISIONS_SMOKE_R2.tsv").read_text(encoding="utf-8").splitlines()
    manifest_lines = (ROOT / "MANIFEST.sha256").read_text(encoding="utf-8").splitlines()
    entries = [line.split("  ", 1) for line in manifest_lines if line]
    manifest_paths = [relative for _, relative in entries]
    actual_paths = sorted(
        path.relative_to(ROOT).as_posix()
        for path in ROOT.rglob("*")
        if path.is_file()
        and path.name != "MANIFEST.sha256"
        and "__pycache__" not in path.relative_to(ROOT).parts
    )
    parent_package = ROOT.parent / parent["code_parent"]

    checks = {
        "design_only_status": config["design_status"] == "FROZEN_PENDING_RUN_APPROVAL",
        "seven_cells": len(config["cells"]) == 7,
        "three_family_replicates": config["replicates_per_cell"] == 3,
        "thirty_nine_resamples": config["resamples"] == 39,
        "four_candidates": config["family_size"] == 4,
        "d29_gap_frozen": config["minimum_fitted_rank_gap"] == 0.10,
        "planned_counts_exact": config["planned_counts"]
        == {
            "cells": 7,
            "independent_families": 21,
            "family_level_resamples": 819,
            "candidate_observed_evaluations": 84,
            "candidate_bootstrap_evaluations": 3276,
        },
        "new_seed_namespaces": len(
            {
                config["data_seed_base"],
                config["selection_seed_base"],
                config["inference_seed_base"],
                config["bootstrap_seed_base"],
            }
        )
        == 4
        and min(
            config["data_seed_base"],
            config["selection_seed_base"],
            config["inference_seed_base"],
            config["bootstrap_seed_base"],
        )
        >= 628000000,
        "decisions_d22_through_d33_frozen": len(decision_lines) == 13
        and [line.split("\t", 1)[0] for line in decision_lines[1:]]
        == [f"D{number}" for number in range(22, 34)]
        and all("\tFROZEN\t" in line for line in decision_lines[1:]),
        "parent_package_present": parent_package.is_dir(),
        "parent_manifest_bound": parent_package.is_dir()
        and sha256_file(parent_package / "MANIFEST.sha256") == parent["code_parent_manifest_sha256"],
        "formal_r1_5_result_hash_bound": parent["r1_5_result_manifest_sha256"]
        == "6081484beb28de8641b606ff89aef4edf12c635e8762e1e94b9fb3f9d160de2b",
        "run_requires_approval_record": "--approval-record" in (ROOT / "scripts" / "run_bounded_smoke.py").read_text(encoding="utf-8"),
        "manifest_has_no_placeholder": all("PLACEHOLDER" not in line for line in manifest_lines),
        "manifest_uses_posix_paths": all(
            "\\" not in relative and PurePosixPath(relative).as_posix() == relative
            for relative in manifest_paths
        ),
        "manifest_paths_complete": sorted(manifest_paths) == actual_paths,
        "manifest_hashes_match": all(
            (ROOT / relative).is_file() and sha256_file(ROOT / relative) == digest
            for digest, relative in entries
        ),
        "text_files_use_lf": all(
            b"\r\n" not in path.read_bytes()
            for path in ROOT.rglob("*")
            if path.is_file() and path.suffix in {".md", ".json", ".py", ".tsv", ".sha256"}
        ),
        "v7_not_authorized": "gpu_or_v7_work" in config["does_not_authorize"],
        "formal_power_not_authorized": "t14_or_t16_power" in config["does_not_authorize"],
    }
    status = "PASS" if all(checks.values()) else "FAIL"
    print(
        json.dumps(
            {
                "package": ROOT.name,
                "status": status,
                "checks": checks,
                "manifest_entries": len(entries),
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0 if status == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
