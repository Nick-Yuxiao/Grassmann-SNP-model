# RUN NOT AUTHORIZED

_Design-only package (Gate 0A). No run has been executed and none is authorized here._

---

## 🚫 Why this directory is empty

Gate 0A is frozen (`GATE_0A_PROTOCOL.md`) but its run code is intentionally not written yet. Execution is blocked behind five pre-run gates:

1. centering AND standardization fixed and verified (training-fold stats only);
2. data-contract hashes bound (frozen v7 chr22 panel + block version);
3. the **detectability gate** fixes the margin and replicate count (20 seeds is likely underpowered at h²=0.05, N≈2247);
4. a compute smoke (3 blocks × 2 k × 1 fold) passes, or the eigensolver fallback is pre-registered;
5. explicit project-lead run authorization.

## 📋 What running it will produce (once authorized)

Per-regime, per-budget `R²_genetic` (primary) and `R²_pheno` (secondary) with replicate-bootstrap CIs, plus the cost ledger. The output is a **per-regime headroom verdict** feeding the gate ladder — never a real-phenotype claim, never a Grassmann-specificity claim, and never a pooled cross-regime average.
