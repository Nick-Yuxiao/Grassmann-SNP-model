#!/usr/bin/env bash
set -euo pipefail
umask 077

ROOT="/data1/home/tanyuxiao/Grassmann_model"
PYTHON_BIN=""
MATRIX_SIZE=2048
ITERATIONS=10

while [[ $# -gt 0 ]]; do
  case "$1" in
    --root) ROOT="$2"; shift 2 ;;
    --python) PYTHON_BIN="$2"; shift 2 ;;
    --matrix-size) MATRIX_SIZE="$2"; shift 2 ;;
    --iterations) ITERATIONS="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done

if [[ "${ROOT}" != "/data1/home/tanyuxiao/Grassmann_model" ]]; then
  echo "Refusing unexpected root: ${ROOT}" >&2
  exit 2
fi
if [[ -z "${PYTHON_BIN}" ]]; then
  PYTHON_BIN="$(command -v python3 || true)"
fi
if [[ -z "${PYTHON_BIN}" || ! -x "${PYTHON_BIN}" ]]; then
  echo "A Python executable with Torch >=2.7/cu128 is required via --python." >&2
  exit 2
fi
if ! command -v nvidia-smi >/dev/null 2>&1; then
  echo "nvidia-smi is unavailable; refusing GPU test." >&2
  exit 3
fi
if ! command -v flock >/dev/null 2>&1; then
  echo "flock is unavailable; refusing an unlocked GPU test." >&2
  exit 3
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_ID="$(date -u +%Y%m%dT%H%M%SZ)_gpu_test_$$"
RESULT_DIR="${ROOT}/v7/results/gpu_test/${RUN_ID}"
RESOURCE_DIR="${ROOT}/v7/resources/gpu_test/${RUN_ID}"
mkdir -p "${RESULT_DIR}" "${RESOURCE_DIR}" "${ROOT}/v7/locks"

nvidia-smi --query-gpu=index,uuid,name,memory.total,memory.used,memory.free,utilization.gpu,temperature.gpu,pstate,driver_version \
  --format=csv,noheader,nounits > "${RESOURCE_DIR}/gpus_before.csv"
nvidia-smi --query-compute-apps=gpu_uuid,pid,process_name,used_memory \
  --format=csv,noheader,nounits > "${RESOURCE_DIR}/compute_processes_before.csv" || true
ps -eo pid,ppid,user,stat,etimes,pcpu,pmem,comm --sort=-pcpu > "${RESOURCE_DIR}/processes_before.txt"

GPU_INDEX="$("${PYTHON_BIN}" - "${RESOURCE_DIR}/gpus_before.csv" "${RESOURCE_DIR}/compute_processes_before.csv" <<'PY'
import csv, sys
gpu_rows = list(csv.reader(open(sys.argv[1], encoding="utf-8"), skipinitialspace=True))
proc_rows = list(csv.reader(open(sys.argv[2], encoding="utf-8"), skipinitialspace=True))
busy = {row[0] for row in proc_rows if row}
for row in gpu_rows:
    index, uuid = row[0], row[1]
    used_mib, utilization = float(row[4]), float(row[6])
    if uuid not in busy and used_mib <= 1024 and utilization <= 5:
        print(index)
        raise SystemExit(0)
raise SystemExit(3)
PY
)" || {
  printf '%s\n' 'NO_IDLE_GPU: no test was started and no process was signalled.' | tee "${RESULT_DIR}/STATUS.txt"
  exit 3
}

LOCK_FILE="${ROOT}/v7/locks/gpu_${GPU_INDEX}.lock"
exec 9>"${LOCK_FILE}"
if ! flock -n 9; then
  printf '%s\n' "PROJECT_LOCK_BUSY for GPU ${GPU_INDEX}: no test was started." | tee "${RESULT_DIR}/STATUS.txt"
  exit 4
fi

# Re-query after taking the project lock. Any process on the selected GPU aborts the test.
nvidia-smi --query-gpu=index,uuid,name,memory.total,memory.used,memory.free,utilization.gpu,temperature.gpu,pstate,driver_version \
  --format=csv,noheader,nounits > "${RESOURCE_DIR}/gpus_recheck.csv"
nvidia-smi --query-compute-apps=gpu_uuid,pid,process_name,used_memory \
  --format=csv,noheader,nounits > "${RESOURCE_DIR}/compute_processes_recheck.csv" || true

"${PYTHON_BIN}" - "${GPU_INDEX}" "${RESOURCE_DIR}/gpus_recheck.csv" "${RESOURCE_DIR}/compute_processes_recheck.csv" <<'PY'
import csv, sys
selected = sys.argv[1]
gpus = {row[0]: row for row in csv.reader(open(sys.argv[2], encoding="utf-8"), skipinitialspace=True)}
processes = list(csv.reader(open(sys.argv[3], encoding="utf-8"), skipinitialspace=True))
row = gpus.get(selected)
if row is None:
    raise SystemExit("selected GPU disappeared")
busy = {proc[0] for proc in processes if proc}
if row[1] in busy or float(row[4]) > 1024 or float(row[6]) > 5:
    raise SystemExit("selected GPU became busy; refusing test")
PY

export CUDA_VISIBLE_DEVICES="${GPU_INDEX}"
TEST_COMMAND=(
  "${PYTHON_BIN}" "${SCRIPT_DIR}/gpu_test_5090.py"
  --output "${RESULT_DIR}/GPU_TEST.json"
  --physical-gpu-index "${GPU_INDEX}"
  --matrix-size "${MATRIX_SIZE}"
  --iterations "${ITERATIONS}"
)
if command -v ionice >/dev/null 2>&1; then
  nice -n 10 ionice -c 2 -n 7 "${TEST_COMMAND[@]}" 2>&1 | tee "${RESULT_DIR}/GPU_TEST.log"
else
  nice -n 10 "${TEST_COMMAND[@]}" 2>&1 | tee "${RESULT_DIR}/GPU_TEST.log"
fi

nvidia-smi --query-gpu=index,uuid,name,memory.total,memory.used,memory.free,utilization.gpu,temperature.gpu,pstate,driver_version \
  --format=csv,noheader,nounits > "${RESOURCE_DIR}/gpus_after.csv"
nvidia-smi --query-compute-apps=gpu_uuid,pid,process_name,used_memory \
  --format=csv,noheader,nounits > "${RESOURCE_DIR}/compute_processes_after.csv" || true
printf '%s\n' "PASS: GPU ${GPU_INDEX}; results=${RESULT_DIR}; resources=${RESOURCE_DIR}" | tee "${RESULT_DIR}/STATUS.txt"
find "${RESULT_DIR}" "${RESOURCE_DIR}" -maxdepth 1 -type f ! -name 'RUN_MANIFEST.sha256' -print0 \
  | sort -z | xargs -0 sha256sum > "${RESULT_DIR}/RUN_MANIFEST.sha256"
