#!/usr/bin/env python3
"""Bind the frozen data contract: compute the real SHA-256 of the panel manifest
and the LD-block version artifact and write them into BINDING.json.

This is the ONLY sanctioned way to populate BINDING.json. It computes hashes from
real files; it never accepts hand-typed hashes and never fabricates. Binding is
all-or-nothing: both artifacts must be provided and must exist.

Run where the frozen artifacts actually live (the server). Then commit the
resulting BINDING.json.

Usage:
    python scripts/bind_data_contract.py \
        --panel-manifest /path/to/frozen_panel_MANIFEST.sha256 \
        --block-version  /path/to/ld_block_version.<ext> \
        [--block-version-id <string>] \
        [--verify-preprocessing] \
        [--bound-by <name>]
"""
from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import sys
from pathlib import Path

PKG = Path(__file__).resolve().parent.parent
BINDING = PKG / "BINDING.json"
EXPECTED_LENGTH = 154850


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _sha256_tree(path: Path) -> str:
    """Deterministic hash of a directory: hash of sorted (relpath, filehash) lines."""
    if path.is_file():
        return _sha256(path)
    lines = []
    for p in sorted(path.rglob("*")):
        if p.is_file():
            lines.append(f"{p.relative_to(path).as_posix()}  {_sha256(p)}")
    return hashlib.sha256("\n".join(lines).encode()).hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel-manifest", required=True)
    ap.add_argument("--block-version", required=True)
    ap.add_argument("--block-version-id", default=None)
    ap.add_argument("--verify-preprocessing", action="store_true")
    ap.add_argument("--bound-by", default=None)
    args = ap.parse_args()

    panel = Path(args.panel_manifest)
    block = Path(args.block_version)
    if not panel.exists():
        sys.exit(f"FAIL: panel manifest not found: {panel}")
    if not block.exists():
        sys.exit(f"FAIL: block version artifact not found: {block}")

    data = json.loads(BINDING.read_text())
    data["panel_manifest"]["path"] = str(panel)
    data["panel_manifest"]["sha256"] = _sha256(panel)
    data["block_version"]["path"] = str(block)
    data["block_version"]["sha256"] = _sha256_tree(block)
    data["block_version"]["id"] = args.block_version_id or block.name
    data["preprocessing_frozen"]["verified"] = bool(args.verify_preprocessing)
    data["bound_by"] = args.bound_by
    data["bound_utc"] = _dt.datetime.now(_dt.timezone.utc).isoformat()

    complete = (
        data["panel_manifest"]["sha256"]
        and data["block_version"]["sha256"]
        and data["preprocessing_frozen"]["verified"]
    )
    data["status"] = "BOUND" if complete else "PARTIAL"

    BINDING.write_text(json.dumps(data, indent=2))
    print(f"panel_manifest.sha256 = {data['panel_manifest']['sha256']}")
    print(f"block_version.sha256  = {data['block_version']['sha256']}")
    print(f"status = {data['status']}")
    if data["status"] != "BOUND":
        print("NOTE: status is PARTIAL until preprocessing is verified "
              "(rerun with --verify-preprocessing once centering/MAF-z is confirmed).")
        if EXPECTED_LENGTH:
            print(f"REMINDER: confirm the panel manifest describes L={EXPECTED_LENGTH} sites.")


if __name__ == "__main__":
    main()
