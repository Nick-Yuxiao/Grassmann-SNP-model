# V7 CPU learning-smoke (ENGINEERING_NON_EVIDENCE)

_Purpose: prove the frozen Grassmann architecture **runs and learns** — forward/backward
work and masked-genotype loss falls well below the `ln(3)` random baseline — on a small,
structured, learnable synthetic task. CPU-only, minutes, bind-only to the frozen model._

---

## What this is (and is not)

| It **is** | It is **not** |
| --- | --- |
| A "does the model train and learn?" engineering check | An architecture comparison / ranking |
| Bound to the frozen `profile_models_v7_1.MaskedGenotypeModel` (no re-implementation) | A read of real 1KGP data, PCs, or any holdout |
| A learnable synthetic DGP so loss *must* drop if the model works | An A1-EFFICIENCY / A1-CAPACITY / GO-NO-GO claim |
| Runnable on CPU (no GPU, no server evidence dir) | Authorization for any GPU / server / evidence-chain run |

Classification is stamped in every output as `ENGINEERING_NON_EVIDENCE` with an explicit
`does_not_authorize` list. Nothing here may enter the v7 evidence chain.

## Why a new smoke instead of the T03 profiler

The T03 profiler (`profile_models_v7_1.py::make_batch`) draws **i.i.d. uniform** genotypes.
Masked cross-entropy on independent uniform targets is pinned at `ln(3) ≈ 1.0986` and cannot
fall — the profiler is a throughput/memory check, not a learning check. This smoke supplies a
**recoverable signal**, so a working model is forced to drive the loss down. That is the
"有反馈 / it actually learns" evidence.

## The learnable DGP

For sample `b`, site `j`:

```
target[b, j] = ( pos_pattern[j mod 256] + shift[b] ) mod 3
```

- `pos_pattern`: fixed random vector over one attention window (period = `attention_window` = 256),
  representable by the model's periodic position embedding.
- `shift[b] ∈ {0,1,2}`: per-sample global offset. Exposed to the PC arm via `pcs[:,0]`, and
  inferable from any unmasked token by the local / grassmann arms — so **all three arms have a
  learnable path**.

## Two deliberate settings (documented traps)

1. **`length = 512` (≥ 2 blocks of 256).** The `GrassmannBlockMixer` global wedge channel is
   **identically zero at L ≤ 256** (a single 256-block makes `context == reduced`, so the
   exterior product cancels). The unit test `test_grassmann_wedge_degenerate_at_one_block`
   pins this. The profiler's CPU dry-run caps L at 256, which would silently disable the
   Grassmann channel — this smoke uses 512 so the channel is live.
2. **`mask_rate = 0.15` (not the protocol's 0.90).** Low masking leaves ample context for a
   fast, clean learning signal. This is a smoke, not the protocol training horizon.

Neither setting is a protocol change; both are local to this non-evidence engineering check.

## Files

```
20260831_grassmann_v7_learning_smoke_cpu/
├── README.md                      # this file
├── SMOKE_PROTOCOL.md              # frozen DGP, success rule, boundaries
├── MANIFEST.sha256                # LF, forward-slash
├── p0/
│   ├── run_learning_smoke.py      # binds to frozen model; trains 3 arms; writes JSON report
│   └── tests/
│       └── test_learning_smoke.py # DGP determinism, shapes, wedge trap, learns-a-little
└── server_ops/
    └── SERVER_STEPS.md            # exact commands to run + verify on the server
```

## Expected feedback

`run_learning_smoke.py` writes `LEARNING_SMOKE_REPORT.json` with, per arm: parameter count,
initial vs final held-out eval loss, min loss, and `learned` (final < success threshold).
A healthy run prints three arms with `learned=True` and terminal status
`LEARNING_SMOKE_PASS` (each final eval loss well under `ln(3) ≈ 1.0986`). See `server_ops`.
