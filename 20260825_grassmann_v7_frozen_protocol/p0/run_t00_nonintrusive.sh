#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
AUDIT="${ROOT}/p0/REMOTE_SAFETY_AUDIT.json"
SERVER_ROOT="${V7_SERVER_ROOT:-/data1/home/tanyuxiao/Grassmann_model}"

# This wrapper never sends signals and never removes another run's files.
python3 "${ROOT}/p0/audit_server.py" --output "${AUDIT}" --require-idle-gpu
GPU_INDEX="$(python3 - "${AUDIT}" <<'PY'
import json, sys
payload = json.load(open(sys.argv[1], encoding="utf-8"))
print(payload["idle_gpu_indices"][0])
PY
)"

python3 "${ROOT}/environment/probe_server.py"
bash "${ROOT}/environment/create_env_5090.sh"

# Re-audit immediately before the CUDA smoke test to reduce the observation/use race.
mkdir -p "${SERVER_ROOT}/v7/locks"
command -v flock >/dev/null 2>&1 || { echo "flock is required for non-interference" >&2; exit 4; }
exec 9>"${SERVER_ROOT}/v7/locks/gpu_${GPU_INDEX}.lock"
flock -n 9 || { echo "Project GPU lock is busy; no smoke test was started." >&2; exit 4; }
python3 "${ROOT}/p0/audit_server.py" --output "${AUDIT}" --require-idle-gpu
python3 - "${AUDIT}" "${GPU_INDEX}" <<'PY'
import json, sys
payload = json.load(open(sys.argv[1], encoding="utf-8"))
if int(sys.argv[2]) not in payload["idle_gpu_indices"]:
    raise SystemExit("originally selected GPU became busy; refusing smoke test")
PY
export CUDA_VISIBLE_DEVICES="${GPU_INDEX}"
source "${V7_ENV_DIR:-${ROOT}/.venv}/bin/activate"
python "${ROOT}/environment/smoke_cuda.py" --output "${ROOT}/environment/ENV_SMOKE.json"

sha256sum \
  "${ROOT}/environment/requirements-cu128.lock" \
  "${ROOT}/environment/ENV_SMOKE.json" \
  "${ROOT}/environment/SERVER_RESOURCE.json" \
  "${AUDIT}" \
  > "${ROOT}/p0/T00_RUNTIME_MANIFEST.sha256"
echo "T00 completed on an idle GPU selected by policy; no existing process was signalled."
