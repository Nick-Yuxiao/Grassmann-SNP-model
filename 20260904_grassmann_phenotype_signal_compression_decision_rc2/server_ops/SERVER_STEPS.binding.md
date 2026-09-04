# Server steps — bind the frozen data contract (pre-run steps 1–3)

_Run on the server where the frozen v7 chr22 panel and LD-block definition live. These steps compute real SHA-256 hashes and produce the values that go into `BINDING.json`. They read only; they do not modify any frozen artifact._

---

## 0. Locate the frozen artifacts (do not guess paths)

```bash
GRASS_ROOT=/data1/home/tanyuxiao/Grassmann_model

# panel manifest for the frozen chr22 donor-only panel (L=154850)
find "$GRASS_ROOT" -type f \( -iname '*MANIFEST*' -o -iname '*panel*manifest*' \) 2>/dev/null | sort

# LD-block version definition over the 154850-site panel
find "$GRASS_ROOT" -type f \( -iname '*block*' -o -iname '*ld_block*' \) 2>/dev/null | sort
```

Confirm the panel manifest describes **L=154850** sites, donor train 2247 / validation 249, before binding. Do not bind a manifest whose length differs.

## 1. Compute the hashes and write BINDING.json

From the checkout of this repo on the server (or after copying this package there):

```bash
V7_PY=<path to the frozen python with no extra deps needed; stdlib only is fine>

python 20260904_grassmann_phenotype_signal_compression_decision_rc2/scripts/bind_data_contract.py \
  --panel-manifest <PANEL_MANIFEST_PATH_FROM_STEP_0> \
  --block-version  <BLOCK_VERSION_PATH_FROM_STEP_0> \
  --block-version-id <e.g. ld_blocks_v1> \
  --bound-by "<your id>"
```

This prints `panel_manifest.sha256`, `block_version.sha256`, and `status`. It will be `PARTIAL` until preprocessing is verified.

## 2. Verify preprocessing, then finalize

Confirm that per-block centering (`G - 2p`) and MAF-z standardization (`(G-2p)/sqrt(2p(1-p))`) use **outer-training-fold statistics only** (no leakage), and that the previously flagged centering defect is fixed. Then rerun with the flag:

```bash
python 20260904_grassmann_phenotype_signal_compression_decision_rc2/scripts/bind_data_contract.py \
  --panel-manifest <...> --block-version <...> --block-version-id <...> \
  --bound-by "<your id>" --verify-preprocessing
```

`status` should now read `BOUND`.

## 3. Return the result

Commit the resulting `BINDING.json` (or paste its `panel_manifest.sha256`, `block_version.sha256`, and `status` back here). Only once `status == BOUND` may the σ-pilot package be built against this contract and run.

> ⚠️ Never hand-type or fabricate a hash. If the frozen artifacts are not reachable, binding does not proceed — the chain stops here rather than running against an unverified contract.
