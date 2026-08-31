from __future__ import annotations

import hashlib
import json
from pathlib import Path, PurePosixPath


ROOT = Path(__file__).resolve().parents[1]


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    config = json.loads((ROOT / "config" / "FEASIBILITY_CONFIG.json").read_text(encoding="utf-8"))
    firewall = json.loads((ROOT / "EVIDENCE_FIREWALL.json").read_text(encoding="utf-8"))
    parent = json.loads((ROOT / "PARENT_EVIDENCE.json").read_text(encoding="utf-8"))
    decisions = (ROOT / "DECISIONS_FEASIBILITY.tsv").read_text(encoding="utf-8").splitlines()
    manifest_lines = (ROOT / "MANIFEST.sha256").read_text(encoding="utf-8").splitlines()
    entries = [line.split("  ", 1) for line in manifest_lines if line]
    manifest_paths = [relative for _, relative in entries]
    static_paths = sorted(
        path.relative_to(ROOT).as_posix()
        for path in ROOT.rglob("*")
        if path.is_file()
        and path.name != "MANIFEST.sha256"
        and "__pycache__" not in path.relative_to(ROOT).parts
        and not ("results" in path.relative_to(ROOT).parts and path.name != "NON_EVIDENCE_NOTICE.md")
    )
    reference = ROOT.parent / parent["implementation_reference"]
    planned_cells = (
        len(config["sample_sizes"])
        * len(config["mafs"])
        * len(config["population_relative_gaps"])
        * len(config["true_max_principal_angles_deg"])
    )
    checks = {
        "classification_non_evidence": config["classification"] == "EXPLORATORY_NON_EVIDENCE",
        "firewall_blocks_promotion": firewall["formal_evidence_eligible"] is False
        and firewall["promotion_permitted"] is False,
        "no_formal_p_values": firewall["formal_p_values_generated"] is False,
        "no_power_ranking": firewall["power_ranking_permitted"] is False,
        "eight_decisions_frozen": len(decisions) == 9
        and [line.split("\t", 1)[0] for line in decisions[1:]] == [f"D{number}" for number in range(34, 42)]
        and all("\tFROZEN\t" in line for line in decisions[1:]),
        "grid_cells_exact": planned_cells == 240,
        "independent_replicates_frozen": config["replicates_per_cell"] == 12,
        "sample_and_gap_thresholds_frozen": config["minimum_group_count"] == 50
        and config["minimum_fitted_rank_gap"] == 0.10,
        "seed_namespace_separate": config["seed_base"] == 629000000,
        "implementation_reference_present": reference.is_dir(),
        "implementation_reference_bound": reference.is_dir()
        and sha256_file(reference / "MANIFEST.sha256")
        == parent["implementation_reference_manifest_sha256"],
        "output_forced_under_non_evidence_results": "results_root not in output.parents"
        in (ROOT / "scripts" / "run_feasibility.py").read_text(encoding="utf-8"),
        "manifest_has_no_placeholder": all("PLACEHOLDER" not in line for line in manifest_lines),
        "manifest_uses_posix_paths": all(
            "\\" not in relative and PurePosixPath(relative).as_posix() == relative
            for relative in manifest_paths
        ),
        "manifest_paths_complete": sorted(manifest_paths) == static_paths,
        "manifest_hashes_match": all(
            (ROOT / relative).is_file() and sha256_file(ROOT / relative) == digest
            for digest, relative in entries
        ),
        "static_text_files_use_lf": all(
            b"\r\n" not in path.read_bytes()
            for path in ROOT.rglob("*")
            if path.is_file()
            and path.suffix in {".md", ".json", ".py", ".tsv", ".sha256"}
            and not ("results" in path.relative_to(ROOT).parts and path.name != "NON_EVIDENCE_NOTICE.md")
        ),
        "v7_not_authorized": "gpu_or_v7_work" in config["does_not_authorize"],
    }
    status = "PASS" if all(checks.values()) else "FAIL"
    print(
        json.dumps(
            {"package": ROOT.name, "status": status, "checks": checks, "manifest_entries": len(entries)},
            indent=2,
            sort_keys=True,
        )
    )
    return 0 if status == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
