"""Load the frozen panel for the sigma pilot from a documented .npy interface.

The frozen v7 chr22 panel lives on the server, not in this repo. To keep the
estimator independent of the server's on-disk materialization format, the pilot
reads a small, explicit set of arrays that the server operator exports ONCE from
the bound frozen panel:

    <data_dir>/G.npy            int8/int16, shape (N_samples, L_sites), dosage 0/1/2
    <data_dir>/blocks.npy       int,        shape (L_sites,), LD-block id per site
    <data_dir>/train_idx.npy    int,        donor-train row indices (2247)
    <data_dir>/val_idx.npy      int,        donor-validation row indices (249)

The export MUST come from the artifact whose SHA-256 is recorded in BINDING.json
(status BOUND). The loader checks shapes and the expected site count; it does not
and cannot re-verify the upstream hash -- that is the binding step's job.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

EXPECTED_L = 154850


def load_frozen_panel(data_dir, expect_L=EXPECTED_L, strict=True):
    d = Path(data_dir)
    G = np.load(d / "G.npy")
    blocks = np.load(d / "blocks.npy")
    train_idx = np.load(d / "train_idx.npy")
    val_idx = np.load(d / "val_idx.npy")

    if G.ndim != 2:
        raise ValueError(f"G must be 2-D (N, L); got {G.shape}")
    if blocks.shape[0] != G.shape[1]:
        raise ValueError(f"blocks length {blocks.shape[0]} != L {G.shape[1]}")
    if strict and expect_L is not None and G.shape[1] != expect_L:
        raise ValueError(
            f"panel has L={G.shape[1]} sites but the contract expects {expect_L}; "
            "refusing to run against a mismatched panel"
        )
    if set(np.unique(G)).difference({0, 1, 2}):
        raise ValueError("G must be dosage-coded 0/1/2")
    return {
        "G": G,
        "blocks": blocks,
        "train_idx": train_idx,
        "val_idx": val_idx,
        "N": int(G.shape[0]),
        "L": int(G.shape[1]),
    }
