# Server steps: v7.2.1 budget bridge

Run each block separately. The launcher is non-interrupting, forbids GPU0,
excludes GPU2, and requires physical GPUs 1,3,4,5,6,7 to be idle.

## Validate R18

```bash
set +e
set +u
set +o pipefail

GRASS_ROOT=/data1/home/tanyuxiao/Grassmann_model
R18="$GRASS_ROOT/v7/code/releases/v7_20260825_p0_r18"
V7_PY="$GRASS_ROOT/v7/code/releases/v7_20260825_p0_r1/.conda-v7/bin/python"

cd "$R18" || exit 1
python3 validate_v7_2_1.py
echo "validator_exit=$?"
sha256sum -c MANIFEST.v7.2.1.sha256
echo "manifest_exit=$?"
"$V7_PY" -m unittest discover -s p0/tests -p 'test_v7_2_1.py' -v
echo "test_exit=$?"
bash -n p0/run_budget_bridge_v7_2_1_nonintrusive.sh
echo "shell_syntax_exit=$?"
```

## Preflight

```bash
LR_PILOT_RUN_ROOT="$GRASS_ROOT/v7/results/lr_pilot/v7.2.0/20260827T082643Z_lr_pilot_v7_2_0_3424359"

cd "$LR_PILOT_RUN_ROOT"
sha256sum -c LR_PILOT_MANIFEST.v7.2.0.sha256
echo "lr_pilot_manifest_exit=$?"
```

```bash
pgrep -a -u "$(id -u)" \
  -f '[r]un_budget_bridge_v7_2_1|[t]rain_budget_bridge_v7_2_1' \
  || echo "none detected"

nvidia-smi --query-gpu=index,memory.used,memory.free,utilization.gpu,pstate --format=csv,noheader,nounits
```

## Launch

Launch only when no existing bridge is detected and GPUs 1,3,4,5,6,7 are idle.

```bash
mkdir -p "$GRASS_ROOT/v7/logs"
STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
BB_LOG="$GRASS_ROOT/v7/logs/budget_bridge_v7_2_1_${STAMP}.log"

nohup env \
  V7_SERVER_ROOT="$GRASS_ROOT" \
  V7_PY="$V7_PY" \
  V7_LR_PILOT_RUN_ROOT="$LR_PILOT_RUN_ROOT" \
  bash "$R18/p0/run_budget_bridge_v7_2_1_nonintrusive.sh" \
  > "$BB_LOG" 2>&1 < /dev/null &

BB_PID=$!
echo "budget_bridge_pid=$BB_PID"
echo "budget_bridge_log=$BB_LOG"
```

## Progress

```bash
BB_RUN_DIR="$(find "$GRASS_ROOT/v7/results/budget_bridge/v7.2.1" \
  -mindepth 1 -maxdepth 1 -type d -name '*_budget_bridge_v7_2_1_*' \
  | sort | tail -n 1)"
echo "budget_bridge_run_dir=$BB_RUN_DIR"
```

```bash
"$V7_PY" -c '
import json,sys
from pathlib import Path
root=Path(sys.argv[1]); steps=[]
for p in root.glob("*/BUDGET_CURVE.jsonl"):
    lines=[line for line in p.read_text().splitlines() if line]
    if lines: steps.append(json.loads(lines[-1])["step"])
print("started_runs:",len(steps))
print("aggregate_steps:",sum(steps),"/",240000)
print("aggregate_percent:",round(100*sum(steps)/240000,2))
print("completed_runs:",len(list(root.glob("*/RESULT.json"))))
print("failed_runs:",len(list(root.glob("*/FAILURE.json"))))
' "$BB_RUN_DIR"

tail -n 40 "$BB_LOG"
```

Assessment exit 0 means 20k adequate; exit 5 means all 12 must later extend to
30k; exit 4 means instability/reproducibility replan. All outputs are retained.

