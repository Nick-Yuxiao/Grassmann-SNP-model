# Server steps — V7 CPU learning-smoke

_CPU-only. Needs **no GPU** and must not inspect, signal, or touch any GPU process. Safe to run
while GPUs 1/3/4/5/6/7 are busy. Uses the frozen `V7_PY` (has torch). Non-evidence._

Run each block from any server directory. Re-define absolute paths in every block; `echo` first.
Use `set +e; set +u; set +o pipefail` in an interactive shell so an expected non-zero exit
(status 4 = INCOMPLETE) does not close the terminal.

---

## 0. Locate the frozen model and pick python

```bash
set +e; set +u; set +o pipefail
GRASS_ROOT=/data1/home/tanyuxiao/Grassmann_model
V7_PY="$GRASS_ROOT/v7/code/releases/v7_20260825_p0_r1/.conda-v7/bin/python"

# find the profile_models_v7_1.py inside the latest release (r22/r23)
MODEL_FILE="$(find "$GRASS_ROOT/v7/code/releases" -name profile_models_v7_1.py 2>/dev/null | sort | tail -n 1)"
echo "model_file=$MODEL_FILE"
export V7_MODEL_DIR="$(dirname "$MODEL_FILE")"
echo "V7_MODEL_DIR=$V7_MODEL_DIR"

"$V7_PY" -c "import torch,sys; print('torch',torch.__version__); sys.path.insert(0,'$V7_MODEL_DIR'); import profile_models_v7_1 as m; print('kinds',m.MODEL_KINDS)"
echo "import_check_exit=$?"
```

## 1. Place the downloaded smoke package

Download this folder from GitHub (`20260831_grassmann_v7_learning_smoke_cpu/`) and put it under
an ingest path — do **not** overwrite any release file.

```bash
GRASS_ROOT=/data1/home/tanyuxiao/Grassmann_model
SMOKE="$GRASS_ROOT/incoming/20260831_grassmann_v7_learning_smoke_cpu"
echo "smoke_dir=$SMOKE"
ls -1 "$SMOKE/p0" "$SMOKE/p0/tests"
```

## 2. Run the unit tests

```bash
set +e; set +u; set +o pipefail
GRASS_ROOT=/data1/home/tanyuxiao/Grassmann_model
V7_PY="$GRASS_ROOT/v7/code/releases/v7_20260825_p0_r1/.conda-v7/bin/python"
SMOKE="$GRASS_ROOT/incoming/20260831_grassmann_v7_learning_smoke_cpu"
# V7_MODEL_DIR must still be exported from step 0
echo "V7_MODEL_DIR=$V7_MODEL_DIR"

cd "$SMOKE" || exit
"$V7_PY" -m unittest discover -s p0/tests -p 'test_learning_smoke.py' -v
echo "test_exit=$?"
```

Expect all four tests to pass (they skip only if torch or the frozen model is unreachable).

## 3. Run the learning-smoke

```bash
set +e; set +u; set +o pipefail
GRASS_ROOT=/data1/home/tanyuxiao/Grassmann_model
V7_PY="$GRASS_ROOT/v7/code/releases/v7_20260825_p0_r1/.conda-v7/bin/python"
SMOKE="$GRASS_ROOT/incoming/20260831_grassmann_v7_learning_smoke_cpu"
STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
OUT="$GRASS_ROOT/incoming/learning_smoke_out/${STAMP}"

"$V7_PY" "$SMOKE/p0/run_learning_smoke.py" \
  --model-dir "$V7_MODEL_DIR" \
  --output-dir "$OUT" \
  --length 512 --steps 300 --batch-size 8 --lr 1e-3 --mask-rate 0.15 --seed 70101
echo "smoke_exit=$?"    # 0 = LEARNING_SMOKE_PASS ; 4 = INCOMPLETE (diagnose, not a verdict)
echo "report=$OUT/LEARNING_SMOKE_REPORT.json"
```

To force the CPU path even on a GPU box, leave `--device cpu` (the default). It never selects a
GPU and never calls into CUDA.

## 4. Read the result

```bash
GRASS_ROOT=/data1/home/tanyuxiao/Grassmann_model
V7_PY="$GRASS_ROOT/v7/code/releases/v7_20260825_p0_r1/.conda-v7/bin/python"
OUT="$(find "$GRASS_ROOT/incoming/learning_smoke_out" -mindepth 1 -maxdepth 1 -type d | sort | tail -n 1)"
echo "out=$OUT"

"$V7_PY" -c 'import json,sys; d=json.load(open(sys.argv[1])); print("status:",d["status"]); print("baseline ln3:",round(d["random_baseline_cross_entropy"],4)); [print(f"  {r[\"model\"]:<24} init={r[\"initial_eval_loss\"]:.4f} final={r[\"final_eval_loss\"]:.4f} learned={r[\"learned\"]}") for r in d["results"]]' \
  "$OUT/LEARNING_SMOKE_REPORT.json"
```

Paste this readout back. A `LEARNING_SMOKE_PASS` with three `learned=True` arms whose final loss
sits well under `1.0986` is the "model runs and learns" feedback. Anything else, paste the JSON
and we diagnose (it stays non-evidence either way).
