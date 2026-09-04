# Server steps — run the Gate 0A sigma pilot (pre-run chain steps 4–8)

_Run on the server after `BINDING.json` reads `status == BOUND` (see `SERVER_STEPS.binding.md`). Estimates variance only._

---

## 4. Export the frozen panel to the .npy interface

From the bound frozen panel (the artifact whose SHA-256 is recorded in `BINDING.json`), write four arrays into one directory `$EXPORT`:

```text
$EXPORT/G.npy          int8/int16 (N_samples, 154850) dosage 0/1/2
$EXPORT/blocks.npy     int (154850,) LD-block id per site (the bound block version)
$EXPORT/train_idx.npy  int donor-train row indices (2247)
$EXPORT/val_idx.npy    int donor-validation row indices (249)
```

Use the same materialization that produced the bound panel; do not re-derive MAF, coordinates, or blocks. The loader refuses any panel whose length is not 154850.

## 5. Run the pilot

```bash
RC2=20260904_grassmann_phenotype_signal_compression_decision_rc2
python 20260904_grassmann_gate0a_sigma_pilot_rc1/scripts/run_sigma_pilot.py \
  --binding "$RC2/BINDING.json" \
  --data-dir "$EXPORT" \
  --out 20260904_grassmann_gate0a_sigma_pilot_rc1/results/sigma_pilot_result.json
```

Needs numpy only. It refuses with `RUN_BLOCKED` unless the binding is `BOUND`. It runs three representative regimes (`spectral_tail_adversarial`, `major_LD_aligned`, `between_block_interaction`) at per-block `k ∈ {8,16}`, `R_pilot=30`, outer 5-fold CV.

## 6–8. Read sigma, map to R, freeze R_formal

The script prints and writes, per `(regime, k)`: `sigma_hat` and `R_needed`, plus `R_formal = max` over cells (never the average — size to the hardest regime). Paste the result JSON back, or commit it. `R_formal` then becomes the frozen Gate 0A replicate count.

If `R_formal` is `None`, at least one primary-relevant cell is underpowered on the candidate grid — increase the candidate `R` grid; never lower the 0.005 margin.

> ⚠️ This pilot estimates variance only. Its output may set `R` and expose compute feasibility. It may NOT change the margin, drop a regime, change `k`/arms/threshold, or reassign which regime is primary — those are frozen.
