# Server steps: v7.2.0 shared LR pilot

Run each block separately. GPU0 is forbidden and GPU2 is excluded. The launcher
requires physical GPUs 1,3,4,5,6,7 to be idle and never signals another process.

## Validate R17

```bash
set +e
set +u
set +o pipefail

GRASS_ROOT=/data1/home/tanyuxiao/Grassmann_model
R17="$GRASS_ROOT/v7/code/releases/v7_20260825_p0_r17"
V7_PY="$GRASS_ROOT/v7/code/releases/v7_20260825_p0_r1/.conda-v7/bin/python"

cd "$R17" || exit 1
python3 validate_v7_2_0.py
echo "validator_exit=$?"
sha256sum -c MANIFEST.v7.2.0.sha256
echo "manifest_exit=$?"
"$V7_PY" -m unittest discover -s p0/tests -p 'test_v7_2_0.py' -v
echo "test_exit=$?"
bash -n p0/run_lr_pilot_v7_2_0_nonintrusive.sh
echo "shell_syntax_exit=$?"
```

## Freeze the validation firewall (CPU only)

```bash
DATA_DIR="$GRASS_ROOT/v7/resources/panels/v7.1.0/a1r/preprocessed/v7.1.10"
VALIDATION_TSV="$GRASS_ROOT/v7/resources/panels/v7.1.0/frozen/v7.1.1/1KGP_DONOR_VALIDATION.tsv"
FIREWALL_DIR="$GRASS_ROOT/v7/resources/panels/v7.1.0/a1r/validation_firewall/v7.2.0"

if [ -e "$FIREWALL_DIR" ]; then
  echo "STOP: validation firewall already exists: $FIREWALL_DIR"
else
  "$V7_PY" "$R17/p0/freeze_lr_validation_firewall_v7_2_0.py" \
    --data-dir "$DATA_DIR" \
    --validation-tsv "$VALIDATION_TSV" \
    --output-dir "$FIREWALL_DIR"
fi
```

```bash
cd "$FIREWALL_DIR"
sha256sum -c VALIDATION_FIREWALL.v7.2.0.sha256
echo "firewall_manifest_exit=$?"
python3 -c 'import json; d=json.load(open("VALIDATION_FIREWALL.v7.2.0.json")); print("status:",d["status"]); print("tuning_count:",d["tuning_count"]); print("decision_holdout_count:",d["decision_holdout_count"]); print("decision_populations_absent:",d["decision_populations_absent"])'
```

## Preflight and launch

```bash
RUNNING="$(pgrep -a -u "$(id -u)" -f '[r]un_lr_pilot_v7_2_0|[t]rain_lr_pilot_v7_2_0' 2>/dev/null)"
if [ -n "$RUNNING" ]; then
  echo "$RUNNING"
  echo "STOP: existing LR pilot detected"
else
  echo "PASS: no existing LR pilot"
fi

nvidia-smi --query-gpu=index,memory.used,memory.free,utilization.gpu,pstate --format=csv,noheader,nounits
```

Launch only if `RUNNING` is empty and physical GPUs 1,3,4,5,6,7 are idle:

```bash
mkdir -p "$GRASS_ROOT/v7/logs"
STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
LR_LOG="$GRASS_ROOT/v7/logs/lr_pilot_v7_2_0_${STAMP}.log"

nohup env \
  V7_SERVER_ROOT="$GRASS_ROOT" \
  V7_PY="$V7_PY" \
  bash "$R17/p0/run_lr_pilot_v7_2_0_nonintrusive.sh" \
  > "$LR_LOG" 2>&1 < /dev/null &

LR_PID=$!
echo "lr_pilot_pid=$LR_PID"
echo "lr_pilot_log=$LR_LOG"
```

## Progress

```bash
RUN_DIR="$(find "$GRASS_ROOT/v7/results/lr_pilot/v7.2.0" -mindepth 1 -maxdepth 1 -type d -name '*_lr_pilot_v7_2_0_*' | sort | tail -n 1)"
echo "run_dir=$RUN_DIR"

"$V7_PY" -c '
import json,sys
from pathlib import Path
root=Path(sys.argv[1]); steps=[]
for p in root.glob("*/PILOT_CURVE.jsonl"):
    rows=[line for line in p.read_text().splitlines() if line]
    if rows: steps.append(json.loads(rows[-1])["step"])
print("started_runs:",len(steps))
print("aggregate_steps:",sum(steps),"/",96000)
print("aggregate_percent:",round(100*sum(steps)/96000,2))
print("completed_runs:",len(list(root.glob("*/RESULT.json"))))
print("failed_runs:",len(list(root.glob("*/FAILURE.json"))))
' "$RUN_DIR"

tail -n 40 "$LR_LOG"
```

The master may exit 4 when the prespecified selector returns REPLAN; that is a
scientific control status, not automatically a software failure. Inspect
`LR_PILOT_DECISION.v7.2.0.json` and the run manifest.

