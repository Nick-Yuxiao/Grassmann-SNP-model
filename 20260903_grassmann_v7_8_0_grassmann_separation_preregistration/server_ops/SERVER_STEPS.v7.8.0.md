# v7.8.0 server steps

Deploy only after the v7.7.7 re-plan contract readiness verifies. Copy validated R41
to R42, overlay this patch with `--strip-components=1`, then run the validator,
manifest verification, and unit tests. This is a CPU-first pre-registration: it
inspects, reserves, and allocates no GPU, trains nothing, and fits no Grassmann.

```
V7_PY="$GRASS_ROOT/v7/code/releases/v7_20260825_p0_r1/.conda-v7/bin/python"
cd "$R42"
"$V7_PY" validate_v7_8_0.py
( cd "$R42" && sha256sum -c MANIFEST.v7.8.0.sha256 )
"$V7_PY" -m unittest discover -s p0/tests -p 'test_v7_8_0.py' -v

REPLAN="<v7.7.7 output dir with LONG_RANGE_DIFFICULTY_REPLAN_READINESS.v7.7.7.json>"
OUT="$GRASS_ROOT/v7/results/grassmann_separation_preregistration/v7.8.0/$(date -u +%Y%m%dT%H%M%SZ)_grassmann_separation_preregistration_v7_8_0"
"$V7_PY" p0/build_separation_prereg_readiness_v7_8_0.py \
  --replan-readiness "$REPLAN/LONG_RANGE_DIFFICULTY_REPLAN_READINESS.v7.7.7.json" \
  --output-dir "$OUT"
( cd "$OUT" && sha256sum -c GRASSMANN_SEPARATION_PREREG_READINESS_MANIFEST.v7.8.0.sha256 )
```

The builder accepts the source only when it is the v7.7.7 re-plan contract readiness
(`gpu_authorized = false`, `grassmann_fitted = false`). It emits the frozen separation
pre-registration and the arm map (fair conventional suite + Grassmann + oracle / local
/ target-shuffled controls); it fits nothing.

Expected terminal status:
`LONG_RANGE_GRASSMANN_SEPARATION_PREREGISTRATION_SIGNED_NO_LAUNCH`. The only next stage
is `IMPLEMENT_V7_8_1_CPU_GRASSMANN_SEPARATION_HARNESS_NO_LAUNCH`, which requires a
CPU-runnable Grassmann module and the fair conventional suite, and still does not
authorize GPU. GPU-scale training is authorized only by a later, separately signed
stage and only if the CPU-first primary separation GO passes with all fairness and
control gates. If no separation appears under a fair, demonstrably converged
conventional suite, close the v7 Grassmann-primary route rather than expand the task
or escalate to GPU.
