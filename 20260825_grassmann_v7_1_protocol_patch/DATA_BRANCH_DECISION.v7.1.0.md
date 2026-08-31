# T01 Branch B decision v7.1.0

- Decision: **CONDITIONAL_GO**
- Decision date: 2026-08-25
- Scope: public 1KGP donor haplotypes, HAPNEST synthetic augmentation, real HGDP
  external evaluation, chr22 first round
- Primary generator: HAPNEST
- Generator sensitivity: msprime population-history simulation
- Real holdout: HGDP only; forbidden as generator donor or preprocessing-fit data
- Activation condition: `PANEL_MANIFEST.v7.1.0.json` passes all disjointness,
  provenance, site-rule and hash checks
- Scientific activation condition: T02 re-sign, T03 GPU1 profile and T04 capacity
  contract all pass

This signs the branch and design constraints, not the eventual panel contents or a
scientific GO. The exact panel remains pending until its immutable files are hashed.
