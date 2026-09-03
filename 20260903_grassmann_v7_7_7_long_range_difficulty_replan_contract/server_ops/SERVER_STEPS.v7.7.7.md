# v7.7.7 server steps

Deploy only after the v7.7.6 execution and its manifest both verify and the result
is `LONG_RANGE_TASK_DIFFICULTY_UNRESOLVED` (`selected_k = null`). Copy validated R40
to R41, overlay this patch with `--strip-components=1`, then run the validator,
manifest verification, and unit tests. This stage is CPU planning-only: it inspects,
reserves, and allocates no GPU, trains nothing, and fits no Grassmann.

```
V7_PY="$GRASS_ROOT/v7/code/releases/v7_20260825_p0_r1/.conda-v7/bin/python"
cd "$R41"
"$V7_PY" validate_v7_7_7.py
( cd "$R41" && sha256sum -c MANIFEST.v7.7.7.sha256 )
"$V7_PY" -m unittest discover -s p0/tests -p 'test_v7_7_7.py' -v

TASK_RUN="<v7.7.6 output dir with TASK_DIFFICULTY_EXECUTION.v7.7.6.json>"
OUT="$GRASS_ROOT/v7/results/long_range_difficulty_replan_contract/v7.7.7/$(date -u +%Y%m%dT%H%M%SZ)_long_range_difficulty_replan_contract_v7_7_7"
"$V7_PY" p0/build_difficulty_replan_readiness_v7_7_7.py \
  --task-difficulty-execution "$TASK_RUN/TASK_DIFFICULTY_EXECUTION.v7.7.6.json" \
  --output-dir "$OUT"
( cd "$OUT" && sha256sum -c LONG_RANGE_DIFFICULTY_REPLAN_READINESS_MANIFEST.v7.7.7.sha256 )
```

The builder accepts the source only when it is a v7.7.6 record with status
`LONG_RANGE_TASK_DIFFICULTY_UNRESOLVED`, `selected_k = null`, `gpu_used = false`, and
`grassmann_fitted = false`. It emits a frozen pre-declared difficulty search grid and
the contract readiness; it runs no fits.

Expected terminal status:
`LONG_RANGE_DIFFICULTY_REPLAN_CONTRACT_SIGNED_NO_LAUNCH`. The only next stage is
`IMPLEMENT_V7_7_8_LONG_RANGE_DIFFICULTY_REPLAN_HARNESS_NO_LAUNCH`, which will run the
pre-declared grid search under the fairness and control rules and select at most one
seed-stable corridor-eligible configuration. If no such configuration exists, close
the v7 Grassmann-primary route rather than expand the task search. No GPU, Grassmann,
A3, HAPNEST, HGDP, phenotype, or biological work is authorized here.
