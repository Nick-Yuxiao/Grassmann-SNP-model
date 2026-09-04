#!/usr/bin/env python3
"""Design-time INTEGRITY validator for the Gate 0A package.

This checks the design package's own internal consistency and manifest ONLY.
It executes no part of the experiment, trains nothing, and touches no data.
A pass means the document set is self-consistent -- NOT that the scientific
design is valid.

Usage:
    python scripts/validate_design_package.py
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

PKG = Path(__file__).resolve().parent.parent

REQUIRED_FILES = [
    "README.md",
    "GATE_0A_PROTOCOL.md",
    "GATE_LADDER_ROADMAP.md",
    "DATA_CONTRACT.json",
    "DECISIONS_COMPRESSION.tsv",
    "EVIDENCE_FIREWALL.json",
    "PARENT_EVIDENCE.json",
    "config/GATE_0A_CONFIG.json",
    "config/DGP_REGIMES.json",
    "results/RUN_NOT_AUTHORIZED.md",
    "scripts/validate_design_package.py",
]

PRIMARY_REGIMES = {
    "major_LD_aligned",
    "spectral_tail",
    "low_MAF",
    "within_block_interaction",
    "between_block_additive",
    "between_block_interaction",
}


def _fail(msg: str) -> None:
    print(f"FAIL: {msg}")
    sys.exit(1)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_files_present() -> None:
    for rel in REQUIRED_FILES:
        if not (PKG / rel).is_file():
            _fail(f"missing required file: {rel}")
    print("OK: all required design files present")


def check_config_consistency() -> None:
    cfg = json.loads((PKG / "config/GATE_0A_CONFIG.json").read_text())
    if cfg.get("classification") != "DESIGN_ONLY_RUN_NOT_AUTHORIZED":
        _fail("config classification must be DESIGN_ONLY_RUN_NOT_AUTHORIZED")
    if cfg.get("gate") != "GATE_0A_nonlinear_compression_headroom":
        _fail("config gate must be GATE_0A_nonlinear_compression_headroom")
    if cfg.get("compression_objective") != "unsupervised_reconstruction":
        _fail("compression_objective must be unsupervised_reconstruction")
    if cfg.get("primary_linear_null") != "block_pca_on_maf_standardized_input":
        _fail("primary_linear_null must be block_pca_on_maf_standardized_input")

    arms = cfg.get("arms", {})
    if arms.get("B_pca_z", {}).get("role") != "PRIMARY_LINEAR_NULL":
        _fail("B_pca_z must be the PRIMARY_LINEAR_NULL")
    if arms.get("B_rand_bilinear", {}).get("role") != "EXPLORATORY_ONLY_no_kill":
        _fail("B_rand_bilinear must be EXPLORATORY_ONLY_no_kill at Gate 0A")

    metrics = cfg.get("metrics", {})
    if "R2_genetic" not in metrics.get("primary_decision_metric", ""):
        _fail("primary decision metric must be R2_genetic")

    kill = cfg.get("kill_criterion", {})
    if "geometry" in json.dumps(kill).lower() and "DELETE" not in json.dumps(kill):
        _fail("the Grassmann geometry-specificity kill must be deleted at Gate 0A")
    if "no_pooling_rule" not in kill:
        _fail("the no-pooling-across-regimes rule must be present")

    summary = cfg.get("primary_summary_across_k", {})
    if summary.get("rule") != "FIXED_BUDGET_POINTS":
        _fail("primary cross-k summary rule must be FIXED_BUDGET_POINTS")

    budget = cfg.get("representation_budget", {})
    if "sum_b min(k, m_b)" not in budget.get("total_budget_definition", ""):
        _fail("total budget must be defined as sum_b min(k, m_b)")

    if not cfg.get("cost_ledger", {}).get("record"):
        _fail("cost ledger (retention x cost axis) must be present")
    print("OK: config consistent (Gate 0A scope, null, no Grassmann kill, budget, cost)")


def check_dgp_regimes() -> None:
    regimes = json.loads((PKG / "config/DGP_REGIMES.json").read_text())
    present = set(regimes.get("primary_regimes", {}).keys())
    if present != PRIMARY_REGIMES:
        _fail(f"primary DGP regimes mismatch: {sorted(present)}")
    if regimes.get("trait_type") != "quantitative_only":
        _fail("round 1 must be quantitative_only")
    print("OK: six pre-registered DGP regimes present (quantitative only)")


def check_firewall_closed() -> None:
    fw = json.loads((PKG / "EVIDENCE_FIREWALL.json").read_text())
    for key in (
        "run_execution_permitted",
        "transformer_training_permitted",
        "gate_0b_1_2_permitted",
        "phenotype_gating_permitted",
        "grassmann_specificity_claim_permitted",
    ):
        if fw.get(key) is not False:
            _fail(f"firewall must set {key}=false")
    print("OK: evidence firewall is closed (design-only, Gate 0A scope)")


def check_no_run_outputs() -> None:
    stray = [
        p.name
        for p in (PKG / "results").iterdir()
        if p.is_file() and p.name != "RUN_NOT_AUTHORIZED.md"
    ]
    if stray:
        _fail(f"unexpected run outputs in results/: {stray}")
    print("OK: no unauthorized run outputs present")


def check_manifest() -> None:
    manifest = PKG / "MANIFEST.sha256"
    if not manifest.is_file():
        print("WARN: MANIFEST.sha256 not yet generated (expected before commit)")
        return
    mismatches = []
    for line in manifest.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        digest, rel = line.split(None, 1)
        target = PKG / rel
        if not target.is_file():
            mismatches.append(f"missing {rel}")
        elif _sha256(target) != digest:
            mismatches.append(f"hash mismatch {rel}")
    if mismatches:
        _fail("manifest verification failed: " + "; ".join(mismatches))
    print("OK: MANIFEST.sha256 verified")


def main() -> None:
    check_files_present()
    check_config_consistency()
    check_dgp_regimes()
    check_firewall_closed()
    check_no_run_outputs()
    check_manifest()
    print("DESIGN_PACKAGE_INTEGRITY_VALID")
    print("NOTE: integrity only; this does NOT assert scientific-design validity")


if __name__ == "__main__":
    main()
