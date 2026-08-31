#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SERVER_ROOT="${V7_SERVER_ROOT:-/data1/home/tanyuxiao/Grassmann_model}"
mkdir -p "${SERVER_ROOT}/v7/resources/audits" "${SERVER_ROOT}/v7/locks"
command -v flock >/dev/null 2>&1 || { echo "flock is required for non-interference" >&2; exit 4; }
RUN_ID="$(date -u +%Y%m%dT%H%M%SZ)_t03_$$"
AUDIT="${SERVER_ROOT}/v7/resources/audits/${RUN_ID}_before.json"
OUT="${1:-${ROOT}/p0/profile_runtime}"

python3 "${ROOT}/p0/audit_server.py" --output "${AUDIT}" --require-idle-gpu
GPU_INDEX="$(python3 - "${AUDIT}" <<'PY'
import json, sys
payload = json.load(open(sys.argv[1], encoding="utf-8"))
print(payload["idle_gpu_indices"][0])
PY
)"

# No kill, pkill, preemption, or overwrite of another run is performed.
if [[ -e "${OUT}" ]]; then
  echo "Refusing to reuse existing output directory: ${OUT}" >&2
  exit 5
fi
exec 9>"${SERVER_ROOT}/v7/locks/gpu_${GPU_INDEX}.lock"
if ! flock -n 9; then
  echo "Project GPU lock is busy; no profiler was started." >&2
  exit 4
fi
RECHECK="${SERVER_ROOT}/v7/resources/audits/${RUN_ID}_recheck.json"
python3 "${ROOT}/p0/audit_server.py" --output "${RECHECK}" --require-idle-gpu
python3 - "${RECHECK}" "${GPU_INDEX}" <<'PY'
import json, sys
payload = json.load(open(sys.argv[1], encoding="utf-8"))
if int(sys.argv[2]) not in payload["idle_gpu_indices"]:
    raise SystemExit("selected GPU became busy; refusing T03")
PY
mkdir -p "${OUT}"
export CUDA_VISIBLE_DEVICES="${GPU_INDEX}"
source "${V7_ENV_DIR:-${ROOT}/.venv}/bin/activate"
exec nice -n 10 python "${ROOT}/p0/profile_models.py" \
  --lengths 8192,131072,262144 \
  --steps 100 \
  --output-dir "${OUT}"
