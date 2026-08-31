# Server steps for v7.2.2

## Preconditions

- Deploy this patch by copying the immutable r19 release to a new release and
  overlaying this archive with one stripped top-level component.
- Verify `MANIFEST.v7.2.2.sha256`, validator, unit tests and shell syntax.
- Set `V7_BUDGET_BRIDGE_RUN_ROOT` to the immutable v7.2.1 run whose decision is
  `BUDGET_BRIDGE_EXTEND_ALL_TO_30K`.
- Confirm no duplicate v7.2.2 task exists.
- Perform a fresh read-only GPU audit. GPUs 1,3,4,5,6,7 must all be idle. GPU0
  and GPU2 are never used.

## Launch

```bash
nohup env \
  V7_SERVER_ROOT="$GRASS_ROOT" \
  V7_PY="$V7_PY" \
  V7_BUDGET_BRIDGE_RUN_ROOT="$BB_RUN_DIR" \
  bash "$RELEASE/p0/run_budget_extension_v7_2_2_nonintrusive.sh" \
  > "$LOG" 2>&1 < /dev/null &
```

The runner independently verifies the source manifest and decision, rechecks
each GPU after acquiring its advisory lock, and refuses to signal or terminate
any process.

## Expected exits

- 0: `BUDGET_EXTENSION_30K_ADEQUATE`.
- 4: lineage, run, sequence or instability failure; replan.
- 6: 30k still fails the terminal-change criterion; replan without automatic extension.
