# Gate ladder roadmap (0A → 0B → 1 → 2)

_On record so the attribution chain is explicit. Only Gate 0A is designed in this package. 0B/1/2 are NOT authorized and NOT designed here._

---

The core error corrected in rc2: three different propositions were collapsed into one experiment. They are separated below. Each gate can only decide its own row, and each is a precondition for the next.

```text
                    GATE 0A  (this package, design only; no Transformer)
      Does nonlinear compression have headroom over PCA?
                         │
         MAF-standardized / raw-centered genotype
                         │
        PCA_z  vs  KPCA / AE   (+ random-proj, random-bilinear, LD-prune, block-mean)
                         │
                 matched total budget, fixed ridge / matched interaction head
                         │
             R²_genetic per pre-registered DGP regime (decided separately)
                         │
             FAIL all regimes ── STOP unsupervised nonlinear path
                         │
             PASS some regimes ─→ carry that headroom forward
                         ▼
                    GATE 0B  (requires training an encoder; NOT no-Transformer)
      Does phenotype-FREE pretraining add information?
                         │
        raw-genotype PCA/KPCA   vs   masked-SSL pretrained embeddings + same compression
                         │
        fairness: encoder trained with NO label leakage; MATCHED total budget
                         │
             FAIL ── no foundation-model value on this axis
                         │
                        PASS
                         ▼
                    GATE 1
      Does GRASSMANN geometry add anything over a generic quadratic feature?
                         │
        actual Grassmann invariant feature   vs   MATCHED random bilinear feature
                         │
             FAIL ── keep generic spectral, drop Grassmann-specific claim
                         │
                        PASS
                         ▼
                    GATE 2
      phenotype gating / global mixer / fusion  (most expensive, least attributable; last)
```

## Why the order is forced

- **0A ≠ 0B.** KPCA on raw genotype cannot stand in for a pretrained encoder's `H`; a negative at 0A does not condemn pretraining, and a positive does not endorse it.
- **0B ≠ 1.** Pretraining value is not Grassmann value; Gate 1 needs an actual Grassmann arm against a matched random bilinear, which Gate 0A does not contain.
- **Gate 2 last.** Phenotype gating is the most expensive and least ablatable module and is easiest to mistake for learned feature selection; it waits until 0A/0B/1 have verdicts.

## Firewall

This roadmap authorizes nothing beyond documenting the ladder. Designing or running Gate 0B, 1, or 2 requires its own prospective, hash-bound protocol and explicit authorization.
