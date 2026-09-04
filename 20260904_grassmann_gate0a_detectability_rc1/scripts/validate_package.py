#!/usr/bin/env python3
"""Integrity validator for the Gate 0A detectability package.

Checks file presence, config consistency, that the run outputs exist and are
consistent with the frozen grids, and the manifest. Does not re-run the sim.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

PKG = Path(__file__).resolve().parent.parent

REQUIRED = [
    "README.md",
    "DETECTABILITY_PROTOCOL.md",
    "DECISIONS_DETECTABILITY.tsv",
    "EVIDENCE_FIREWALL.json",
    "PARENT_EVIDENCE.json",
    "config/DETECTABILITY_CONFIG.json",
    "src/detectability/__init__.py",
    "src/detectability/simulator.py",
    "scripts/run_detectability.py",
    "scripts/validate_package.py",
    "tests/test_detectability.py",
    "results/PLANNING_PROXY_NOTICE.md",
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
    cfg = json.loads((PKG / "config/DETECTABILITY_CONFIG.json").read_text())
    if cfg.get("classification") != "PLANNING_PROXY":
        _fail("classification must be PLANNING_PROXY")
    m = cfg["margins"]
    if m.get("primary_practical_margin") != 0.005 or not m.get("primary_is_fixed_a_priori"):
        _fail("primary margin must be fixed a priori at 0.005")
    if 0.005 not in cfg["grids"]["delta_grid"]:
        _fail("delta grid must contain the primary margin 0.005")
    if 0.0 not in cfg["grids"]["delta_grid"]:
        _fail("delta grid must contain 0.0 for the FPR check")
    print("OK: config consistent (fixed 0.005 margin, delta grid includes 0 and 0.005)")


def check_results():
    j = PKG / "results/detectability_table.json"
    if not j.is_file():
        print("WARN: results/detectability_table.json not present (run scripts/run_detectability.py)")
        return
    payload = json.loads(j.read_text())
    cfg = json.loads((PKG / "config/DETECTABILITY_CONFIG.json").read_text())
    sig = set(str(s) for s in cfg["grids"]["sigma_grid"])
    if set(payload["min_R_for_0.90_power_at_0.005_per_sigma"].keys()) != sig:
        _fail("results sigma keys do not match config sigma grid")
    # FPR sanity: one-sided decision on a 95% CI should be well under 0.10
    for rec in payload["records"]:
        if rec["FPR@0"] > 0.10:
            _fail(f"FPR@0 too high ({rec['FPR@0']}) at sigma={rec['sigma']} R={rec['R']}")
    print("OK: results present and consistent with frozen grids (FPR@0 within bound)")


def check_firewall():
    fw = json.loads((PKG / "EVIDENCE_FIREWALL.json").read_text())
    if fw.get("formal_gate_0A_evidence_permitted") is not False:
        _fail("firewall must forbid formal Gate 0A evidence use")
    if fw.get("real_data_used") is not False:
        _fail("firewall must assert no real data used")
    print("OK: firewall consistent (planning proxy, no real data)")


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
    check_results()
    check_firewall()
    check_manifest()
    print("DETECTABILITY_PACKAGE_INTEGRITY_VALID")


if __name__ == "__main__":
    main()
