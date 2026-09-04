#!/usr/bin/env python3
"""Integrity validator for the Gate 0A sigma-pilot package. Does not run the sim."""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

PKG = Path(__file__).resolve().parent.parent

REQUIRED = [
    "README.md",
    "SIGMA_PILOT_PROTOCOL.md",
    "DECISIONS_SIGMA_PILOT.tsv",
    "EVIDENCE_FIREWALL.json",
    "PARENT_EVIDENCE.json",
    "config/SIGMA_PILOT_CONFIG.json",
    "src/gate0a/__init__.py",
    "src/gate0a/estimator.py",
    "src/gate0a/loader.py",
    "src/gate0a/sigma_map.py",
    "scripts/run_sigma_pilot.py",
    "scripts/write_export_manifest.py",
    "scripts/validate_package.py",
    "tests/test_gate0a.py",
    "server_ops/SERVER_STEPS.sigma_pilot.md",
    "results/RUN_BLOCKED_UNTIL_BOUND.md",
]


def _fail(m):
    print(f"FAIL: {m}")
    sys.exit(1)


def _sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def check_files():
    for rel in REQUIRED:
        if not (PKG / rel).is_file():
            _fail(f"missing {rel}")
    print("OK: required files present")


def check_config():
    cfg = json.loads((PKG / "config/SIGMA_PILOT_CONFIG.json").read_text())
    if cfg.get("classification") != "VARIANCE_CALIBRATION_ONLY":
        _fail("classification must be VARIANCE_CALIBRATION_ONLY")
    if cfg["requires_before_run"]["binding_status"] != "BOUND":
        _fail("must require binding status BOUND before run")
    if "spectral_tail_adversarial" not in cfg["representative_regimes"]:
        _fail("spectral_tail_adversarial must be included")
    if "max" not in cfg["replicate_count_rule"]["R_formal"].lower():
        _fail("R_formal must be the max over cells, not the average")
    if cfg.get("inner_cv_folds") != 5:
        _fail("inner_cv_folds must be 5 (real nested CV)")
    if cfg.get("polygenic_background_fraction") is None:
        _fail("polygenic_background_fraction must be set")
    if "train" not in str(cfg.get("trained_on", "")).lower():
        _fail("config must state training-only population (val untouched)")
    mustnot = set(cfg["firewall_results_may_not_decide"])
    for k in ("margin", "k", "arm", "success_threshold", "which_regime_is_primary"):
        if k not in mustnot:
            _fail(f"firewall must forbid deciding {k}")
    print("OK: config consistent (variance-only, BOUND-gated, R_formal=max)")


def check_firewall():
    fw = json.loads((PKG / "EVIDENCE_FIREWALL.json").read_text())
    if fw.get("run_permitted_only_when_bound") is not True:
        _fail("firewall must gate run on BOUND")
    if set(fw.get("results_may_only_decide", [])) != {"replicate_count_R", "compute_feasibility"}:
        _fail("results_may_only_decide must be exactly replicate_count_R + compute_feasibility")
    print("OK: firewall consistent (variance calibration only)")


def check_no_result_before_bound():
    if (PKG / "results" / "sigma_pilot_result.json").is_file():
        # allowed only if a real bound run produced it; flag for human confirmation
        print("NOTE: sigma_pilot_result.json present -- confirm it came from a BOUND run")
    else:
        print("OK: no pilot result yet (expected until a BOUND server run)")


def check_manifest():
    man = PKG / "MANIFEST.sha256"
    if not man.is_file():
        print("WARN: MANIFEST.sha256 not yet generated")
        return
    bad = []
    for line in man.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        digest, rel = line.split(None, 1)
        t = PKG / rel
        if not t.is_file():
            bad.append(f"missing {rel}")
        elif _sha(t) != digest:
            bad.append(f"hash mismatch {rel}")
    if bad:
        _fail("manifest: " + "; ".join(bad))
    print("OK: MANIFEST.sha256 verified")


def main():
    check_files()
    check_config()
    check_firewall()
    check_no_result_before_bound()
    check_manifest()
    print("SIGMA_PILOT_PACKAGE_INTEGRITY_VALID")


if __name__ == "__main__":
    main()
