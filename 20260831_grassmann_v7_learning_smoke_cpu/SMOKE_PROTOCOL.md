# Learning-smoke protocol (frozen before running)

_Classification: `ENGINEERING_NON_EVIDENCE`. This document freezes the task, success rule and
boundaries **before** any numbers are produced, so the run cannot be tuned to a desired result._

---

## 1. Object under test

The **frozen** architecture in `profile_models_v7_1.py`, imported (not copied):
`Architecture()`, `MaskedGenotypeModel(kind, cfg)` for each `kind` in

```
local_attn_8m_w256
local_attn_gpc_8m_w256
grassmann_full_8m_w256
```

No architecture parameter is changed. `cfg.validate()` must pass (`attention_window == 256`, etc.).

## 2. Data-generating process (frozen)

Per sample `b`, site `j`, with `period = cfg.attention_window = 256`:

```
pos_pattern ~ fixed U{0,1,2}^period      (one draw, seed = --seed)
shift[b]    ~ U{0,1,2}                    (per sample)
target[b,j] = (pos_pattern[j mod period] + shift[b]) mod 3
tokens      = target, with a mask_rate fraction set to the mask token (== genotype_states = 3)
pcs[b]      = 0.10 * N(0, I_16);  pcs[b,0] := shift[b] - 1     (clean global signal on ch. 0)
```

Loss: `cross_entropy(logits[mask], target[mask])`. Held-out eval uses a fixed batch
(`seed + 777777`) never trained on.

## 3. Frozen run settings (defaults)

| Setting | Value | Note |
| --- | ---: | --- |
| `length` | 512 | ≥ 2 blocks ⇒ Grassmann wedge channel non-degenerate |
| `steps` | 300 | per arm |
| `batch_size` | 8 | |
| `lr` | 1e-3 (AdamW) | smoke LR, **not** the protocol `4e-4` |
| `mask_rate` | 0.15 | smoke masking, **not** the protocol `0.90` |
| `seed` | 70101 | fixes `pos_pattern` and batch stream |
| `success_threshold` | 0.30 | per-arm final held-out eval loss |
| `device` | cpu | cuda allowed but unnecessary |

## 4. Success rule (frozen)

- Random baseline (no learning) = `ln(genotype_states) = ln(3) ≈ 1.0986`.
- An arm **learned** iff its final held-out eval loss `< success_threshold` (0.30).
- Terminal `status`:
  - `LEARNING_SMOKE_PASS` — all three arms learned.
  - `LEARNING_SMOKE_INCOMPLETE` — at least one arm did not reach the threshold (exit code 4).

`INCOMPLETE` is not a model failure verdict; it is a smoke that did not converge under these
settings and should be diagnosed (steps/lr/mask_rate), never re-labeled as an architecture result.

## 5. Boundaries — `does_not_authorize`

- architecture comparison or ranking (the arms are not matched for a fair contest here)
- any A1-EFFICIENCY / A1-CAPACITY / architecture GO–NO-GO claim
- any GPU, server evidence directory, or 1KGP / holdout / HGDP read
- promotion of any number here into the v7 evidence chain

## 6. What a reader may conclude

Only: "the frozen model runs forward/backward on CPU and reduces masked-genotype loss well
below the random baseline on learnable structure." Nothing about real genomes, LD, imputation
headroom, or which architecture is better.
