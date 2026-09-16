# GEUVADIS chr19 development Pilot rc1

This directory executes the frozen 0526 Stage 2 bridge protocol in
`../../protocols/20260915_0526_stage2_bridge_rc2/`.

Primary estimand: held-out predictive-R2 contrast `C - B`, where `B` is the
same-panel additive dosage ridge and `C` adds the complete coordinate-aligned
hidden state from a genotype-only pretrained and then frozen encoder.

Run exactly once from the repository root:

```powershell
python snp_spectral_foundation/pilots/20260915_geuvadis_chr19_development_rc1/run_pilot.py
```

Generated artifacts go to `results/`. A completed `RESULTS.json` makes the runner
refuse a second execution. This is a real molecular-phenotype development Pilot,
not an external replication or a claim that Grassmann is the preferred geometry.
