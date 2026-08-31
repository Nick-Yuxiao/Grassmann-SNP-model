#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SERVER_ROOT="${V7_SERVER_ROOT:-/data1/home/tanyuxiao/Grassmann_model}"
GPU_INDEX="${V7_GPU_INDEX:-1}"
V7_PY="${V7_PY:-${ROOT}/.conda-v7/bin/python}"
: "${V7_CHR22_L:?Set V7_CHR22_L to the data-derived chr22 variant count from PANEL_MANIFEST.v7.1.0.json}"

if [[ "${GPU_INDEX}" == "0" ]]; then
  echo "REFUSE: physical GPU 0 is forbidden by the frozen v7.1.0 policy." >&2
  exit 8
fi
if ! [[ "${GPU_INDEX}" =~ ^[0-9]+$ && "${V7_CHR22_L}" =~ ^[1-9][0-9]*$ ]]; then
  echo "GPU index and V7_CHR22_L must be positive integers." >&2
  exit 8
fi
if [[ ! -x "${V7_PY}" ]]; then
  echo "Python is not executable: ${V7_PY}" >&2
  exit 8
fi
command -v flock >/dev/null 2>&1 || { echo "flock is required" >&2; exit 8; }

mkdir -p "${SERVER_ROOT}/v7/resources/audits" "${SERVER_ROOT}/v7/locks"
RUN_ID="$(date -u +%Y%m%dT%H%M%SZ)_t03_v7_1_gpu${GPU_INDEX}_$$"
AUDIT_BEFORE="${SERVER_ROOT}/v7/resources/audits/${RUN_ID}_before.json"
AUDIT_RECHECK="${SERVER_ROOT}/v7/resources/audits/${RUN_ID}_recheck.json"
OUT_PARENT="${1:-${ROOT}/p0/profile_runtime_v7_1}"
OUT="${OUT_PARENT}/${RUN_ID}"

python3 "${ROOT}/p0/audit_server.py" --output "${AUDIT_BEFORE}"
python3 - "${AUDIT_BEFORE}" "${GPU_INDEX}" <<'PY'
import json, sys
p = json.load(open(sys.argv[1], encoding="utf-8"))
gpu = int(sys.argv[2])
if gpu not in p.get("idle_gpu_indices", []):
    raise SystemExit(f"REFUSE: requested physical GPU {gpu} is not idle by policy")
PY

exec 9>"${SERVER_ROOT}/v7/locks/gpu_${GPU_INDEX}.lock"
if ! flock -n 9; then
  echo "REFUSE: project lock for GPU ${GPU_INDEX} is busy." >&2
  exit 9
fi

python3 "${ROOT}/p0/audit_server.py" --output "${AUDIT_RECHECK}"
python3 - "${AUDIT_RECHECK}" "${GPU_INDEX}" <<'PY'
import json, sys
p = json.load(open(sys.argv[1], encoding="utf-8"))
gpu = int(sys.argv[2])
if gpu not in p.get("idle_gpu_indices", []):
    raise SystemExit(f"REFUSE: physical GPU {gpu} became busy before launch")
PY

if [[ -e "${OUT}" ]]; then
  echo "REFUSE: output already exists: ${OUT}" >&2
  exit 9
fi
export CUDA_VISIBLE_DEVICES="${GPU_INDEX}"
export PYTHONUNBUFFERED=1
echo "Starting non-interrupting T03 on physical GPU ${GPU_INDEX}; process-local device will be cuda:0."
exec nice -n 10 "${V7_PY}" "${ROOT}/p0/profile_models_v7_1.py" \
  --lengths "8192,${V7_CHR22_L}" \
  --steps 100 \
  --profile-mask-rate 0.90 \
  --output-dir "${OUT}"
