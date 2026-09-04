#!/usr/bin/env python3
"""Write EXPORT_MANIFEST.json tying an exported .npy panel to the BOUND contract.

Run after exporting the frozen panel to <data_dir>/{G,blocks,train_idx,val_idx}.npy.
This records the source panel/block SHA-256 (copied from BINDING.json), the SHA-256
of each exported array, and the counts, and checks them against the contract
(L=154850, train=2247, val=249). run_sigma_pilot.py refuses to run unless this
manifest is present and consistent, so the pilot cannot silently run on an export
that was not derived from the bound artifacts.

Usage:
    python scripts/write_export_manifest.py \
        --binding /path/to/rc2/BINDING.json \
        --data-dir /path/to/exported_panel_npy
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

EXPECTED = {"L": 154850, "train_n": 2247, "val_n": 249}
FILES = ["G.npy", "blocks.npy", "train_idx.npy", "val_idx.npy"]


def _sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--binding", required=True)
    ap.add_argument("--data-dir", required=True)
    args = ap.parse_args()

    b = json.loads(Path(args.binding).read_text())
    if b.get("status") != "BOUND":
        sys.exit(f"FAIL: BINDING status={b.get('status')} (need BOUND)")

    d = Path(args.data_dir)
    for f in FILES:
        if not (d / f).is_file():
            sys.exit(f"FAIL: missing exported array {f}")

    G = np.load(d / "G.npy", mmap_mode="r")
    blocks = np.load(d / "blocks.npy")
    train_idx = np.load(d / "train_idx.npy")
    val_idx = np.load(d / "val_idx.npy")

    counts = {"N": int(G.shape[0]), "L": int(G.shape[1]),
              "train_n": int(train_idx.shape[0]), "val_n": int(val_idx.shape[0])}
    problems = []
    if counts["L"] != EXPECTED["L"]:
        problems.append(f"L={counts['L']} != {EXPECTED['L']}")
    if counts["train_n"] != EXPECTED["train_n"]:
        problems.append(f"train_n={counts['train_n']} != {EXPECTED['train_n']}")
    if counts["val_n"] != EXPECTED["val_n"]:
        problems.append(f"val_n={counts['val_n']} != {EXPECTED['val_n']}")
    if blocks.shape[0] != counts["L"]:
        problems.append(f"blocks length {blocks.shape[0]} != L {counts['L']}")
    if set(np.unique(train_idx)).intersection(set(np.unique(val_idx))):
        problems.append("train_idx and val_idx overlap")
    if problems:
        sys.exit("FAIL: export does not match contract: " + "; ".join(problems))

    manifest = {
        "classification": "EXPORT_PROVENANCE",
        "source_panel_manifest_sha256": b["panel_manifest"]["sha256"],
        "source_block_version_sha256": b["block_version"]["sha256"],
        "source_block_version_id": b["block_version"]["id"],
        "arrays": {f: _sha(d / f) for f in FILES},
        "counts": counts,
        "note": "run_sigma_pilot.py verifies this against BINDING.json and the arrays before running",
    }
    (d / "EXPORT_MANIFEST.json").write_text(json.dumps(manifest, indent=2))
    print("wrote", d / "EXPORT_MANIFEST.json")
    print(json.dumps(manifest["counts"], indent=2))


if __name__ == "__main__":
    main()
