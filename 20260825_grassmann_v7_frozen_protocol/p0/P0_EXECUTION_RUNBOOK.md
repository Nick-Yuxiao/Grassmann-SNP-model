# V7 rapid P0 execution runbook

Goal: obtain a defensible **GO_TO_A0 / NOT_READY** decision in roughly three
days without interrupting any existing server task. This is not the A1
scientific GO/NO-GO.

## Non-interference rule

Every GPU action starts with `audit_server.py`. A GPU is eligible only when it
has no listed compute process, used memory is at most 1024 MiB and utilization
is at most 5%. The scripts never call `kill`, `pkill`, scheduler cancellation,
preemption, or destructive cleanup. Existing output directories are not
reused. If no GPU is eligible, the correct action is to wait and re-audit.

## Upload

Upload this entire frozen-protocol directory to a new server directory. Do not
overlay an existing experiment directory. After transfer, compare the local
and remote hashes for the frozen files with `sha256sum -c MANIFEST.sha256`.

## Day 1: T00, T01 and T02

1. Run `bash p0/run_t00_nonintrusive.sh`. It performs the read-only audit,
   creates the isolated virtual environment and runs the CUDA forward/backward
   smoke test on an eligible GPU.
2. Freeze the actual data branch with `build_panel_manifest.py`. Branch B must
   include an explicitly named synthetic-generator/config input. Place outputs
   in `p0/runtime_t01/`.
3. Independently compare `DECISIONS.v7.0.1.tsv` with
   `METRIC_DEFINITIONS.md`. Do not tune constants from pilot outcomes.

## Day 2: T03

Run `bash p0/run_t03_nonintrusive.sh`. It profiles 100 synchronized steps for
all three model families at L=8192, 131072 and 262144, performs an OOM-boundary
search, and includes data generation/transfer plus checkpoint timing. The
default output directory is new and immutable after completion.

The model implementations are profiler candidates, not yet the T05 production
harness. Their parameter counts are kept near the frozen 8M budget; the report
must retain exact counts so matched-parameter fairness can be audited.

## Day 3: T04 and decision

The frozen protocol does not specify training steps/epochs per confirmatory
run. Do not invent this after seeing performance. First create a signed
training-budget record, then run `build_compute_contract.py` with that record,
the actual signed GPU-hours/storage/window and T03 report. It requires exact
L=131072 and L=262144 measurements and applies both the 5x2 repeats and 2x
engineering margin.

Finally run:

```bash
python p0/assess_p0.py --output-dir p0/runtime_verdict
```

`GO_TO_A0` means environment, data, metrics, measurements and capacity are
ready. Missing evidence yields `NOT_READY`, not `NO_GO`. A scientific Go needs
the preregistered A1 paired results and confidence intervals.
