#!/usr/bin/env python3
"""Design-time validator for the phenotype-signal compression decision package.

This checks the *design* package's own integrity only. It executes no part of
the decision experiment, trains nothing, and touches no data. It verifies that
the frozen design files are present, internally consistent, and match the
committed MANIFEST.sha256.

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
    "COMPRESSION_DECISION_PROTOCOL.md",
    "DATA_CONTRACT.json",
    "DECISIONS_COMPRESSION.tsv",
    "EVIDENCE_FIREWALL.json",
    "PARENT_EVIDENCE.json",
    "config/COMPRESSION_DECISION_CONFIG.json",
    "results/RUN_NOT_AUTHORIZED.md",
    "scripts/validate_design_package.py",
]


def _fail(msg: str) -> None:
    print(f"FAIL: {msg}")
    sys.exit(1)


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def check_files_present() -> None:
    for rel in REQUIRED_FILES:
        if not (PKG / rel).is_file():
            _fail(f"missing required file: {rel}")
    print("OK: all required design files present")


def check_config_consistency() -> None:
    cfg = json.loads((PKG / "config/COMPRESSION_DECISION_CONFIG.json").read_text())
    if cfg.get("classification") != "DESIGN_ONLY_RUN_NOT_AUTHORIZED":
        _fail("config classification must be DESIGN_ONLY_RUN_NOT_AUTHORIZED")
    if cfg.get("compression_objective") != "unsupervised_reconstruction":
        _fail("config compression_objective must be unsupervised_reconstruction")
    if cfg.get("correct_linear_null") != "block_pca":
        _fail("reconstruction objective requires block_pca as the linear null")
    arms = cfg.get("arms", {})
    if arms.get("B_pca", {}).get("role") != "PRIMARY_LINEAR_NULL":
        _fail("B_pca must be the PRIMARY_LINEAR_NULL")
    if "B_rand_bilinear" not in arms:
        _fail("the reinstated random-bilinear ablation arm is required")
    kc = cfg.get("kill_criteria", {})
    for key in ("primary_stage_1_2_kill", "geometry_specificity_kill"):
        if not kc.get(key):
            _fail(f"missing pre-registered kill criterion: {key}")
    print("OK: config internally consistent (objective, null, ablation, kills)")


def check_firewall_closed() -> None:
    fw = json.loads((PKG / "EVIDENCE_FIREWALL.json").read_text())
    for key in (
        "run_execution_permitted",
        "transformer_training_permitted",
        "gpu_or_v7_work_permitted",
        "phenotype_gating_permitted",
    ):
        if fw.get(key) is not False:
            _fail(f"firewall must set {key}=false in a design-only package")
    print("OK: evidence firewall is closed (design-only)")


def check_no_run_outputs() -> None:
    results = PKG / "results"
    stray = [
        p.name
        for p in results.iterdir()
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
    check_firewall_closed()
    check_no_run_outputs()
    check_manifest()
    print("DESIGN_PACKAGE_VALID")


if __name__ == "__main__":
    main()
