# Server steps for v7.2.4

## Preconditions

- Copy immutable release r21 to a new release r22, then overlay this patch with one
  stripped top-level component. Never modify r21 or the completed v7.2.3 run.
- Verify the v7.2.4 manifest, validator, unit tests and shell syntax from r22.
- Set `V7_BUDGET_EXTENSION_30K_RUN_ROOT` to the immutable completed v7.2.3 run:
  `/data1/home/tanyuxiao/Grassmann_model/v7/results/budget_extension/v7.2.3/20260828T023841Z_budget_extension_v7_2_3_3758608`.
- Confirm the v7.2.3 result manifest still verifies and no duplicate v7.2.4 task exists.
- Perform a fresh read-only GPU audit. GPUs 1,3,4,5,6,7 must all be idle. GPU0 and
  GPU2 are forbidden.

## Launch

```bash
nohup env \
  V7_SERVER_ROOT="$GRASS_ROOT" \
  V7_PY="$V7_PY" \
  V7_BUDGET_EXTENSION_30K_RUN_ROOT="$RUN_30K" \
  bash "$R22/p0/run_final_budget_v7_2_4_nonintrusive.sh" \
  > "$FINAL_LOG" 2>&1 < /dev/null &
```

The runner resumes only the six selected-LR cells, one per allowed GPU. It verifies
all immutable input manifests, restores model/optimizer/RNG state, rechecks every GPU
after acquiring its advisory lock, and never signals or terminates another process.

## Expected exits

- 0: all six cells are operationally adequate with no prospective shape flag. A
  separate formal-schedule contract may be prepared; it is not launched here.
- 4: run, lineage, sequence, degradation or other instability failure; diagnose.
- 6: one or more cells remain operationally inadequate at 40k. Stop; no automatic
  extension is authorized.
- 7: all primary criteria pass, but at least one cell is stable-but-accelerating.
  Stop for shape review; no schedule is frozen automatically.
- 8 or 9: missing precondition or advisory-lock refusal. Nothing should be killed;
  inspect the audit and retry only after the resource conflict clears.

No exit permits an architecture decision, HGDP access, holdout access or an automatic
continuation beyond step 40,000.
