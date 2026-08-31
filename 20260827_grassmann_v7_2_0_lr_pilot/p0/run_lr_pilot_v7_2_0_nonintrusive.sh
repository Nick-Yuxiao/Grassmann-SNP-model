#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
GRASS_ROOT="${V7_SERVER_ROOT:-/data1/home/tanyuxiao/Grassmann_model}"
V7_PY="${V7_PY:-${GRASS_ROOT}/v7/code/releases/v7_20260825_p0_r1/.conda-v7/bin/python}"
DATA_DIR="${V7_DATA_DIR:-${GRASS_ROOT}/v7/resources/panels/v7.1.0/a1r/preprocessed/v7.1.10}"
FIREWALL_DIR="${V7_LR_FIREWALL_DIR:-${GRASS_ROOT}/v7/resources/panels/v7.1.0/a1r/validation_firewall/v7.2.0}"
SCHEDULE="${V7_LR_PILOT_SCHEDULE:-${ROOT}/LR_PILOT_SCHEDULE.v7.2.0.tsv}"
OUT_PARENT="${V7_LR_PILOT_OUT_PARENT:-${GRASS_ROOT}/v7/results/lr_pilot/v7.2.0}"
GPUS=(1 3 4 5 6 7)

for path in "${V7_PY}" "${SCHEDULE}" "${DATA_DIR}/PREPROCESS_MANIFEST.v7.1.10.sha256" "${FIREWALL_DIR}/VALIDATION_FIREWALL.v7.2.0.sha256"; do
  [[ -e "${path}" ]] || { echo "missing: ${path}" >&2; exit 8; }
done
mkdir -p "${OUT_PARENT}" "${GRASS_ROOT}/v7/locks" "${GRASS_ROOT}/v7/resources/audits"
(
  cd "${DATA_DIR}"
  sha256sum -c PREPROCESS_MANIFEST.v7.1.10.sha256
)
(
  cd "${FIREWALL_DIR}"
  sha256sum -c VALIDATION_FIREWALL.v7.2.0.sha256
)

RUN_ID="$(date -u +%Y%m%dT%H%M%SZ)_lr_pilot_v7_2_0_$$"
OUT="${OUT_PARENT}/${RUN_ID}"
mkdir "${OUT}"
AUDIT="${OUT}/GPU_AUDIT_BEFORE.json"
python3 "${ROOT}/p0/audit_server.py" --output "${AUDIT}"
python3 - "${AUDIT}" <<'PY'
import json,sys
payload=json.load(open(sys.argv[1]))
required={1,3,4,5,6,7}
idle=set(payload.get('idle_gpu_indices',[]))
if 0 in required or 2 in required: raise SystemExit('forbidden GPU in required set')
if not required <= idle: raise SystemExit(f'REFUSE: required GPUs not idle: {sorted(required-idle)}')
PY

pids=()
worker_failed=0
for gpu in "${GPUS[@]}"; do
  (
    exec 9>"${GRASS_ROOT}/v7/locks/gpu_${gpu}.lock"
    flock -n 9 || exit 9
    RECHECK="${OUT}/GPU_${gpu}_RECHECK.json"
    python3 "${ROOT}/p0/audit_server.py" --output "${RECHECK}"
    python3 - "${RECHECK}" "${gpu}" <<'PY'
import json,sys
payload=json.load(open(sys.argv[1])); gpu=int(sys.argv[2])
if gpu not in payload.get('idle_gpu_indices',[]): raise SystemExit(f'REFUSE: GPU {gpu} became busy')
PY
    export CUDA_VISIBLE_DEVICES="${gpu}" PYTHONUNBUFFERED=1
    nice -n 10 "${V7_PY}" "${ROOT}/p0/run_lr_pilot_gpu_worker_v7_2_0.py" \
      --schedule "${SCHEDULE}" --gpu "${gpu}" --data-dir "${DATA_DIR}" \
      --firewall-dir "${FIREWALL_DIR}" --output-root "${OUT}" \
      --trainer "${ROOT}/p0/train_lr_pilot_v7_2_0.py"
  ) >"${OUT}/GPU_${gpu}_WORKER.log" 2>&1 &
  pids+=("$!")
done
for pid in "${pids[@]}"; do
  wait "${pid}" || worker_failed=1
done

selector_status=0
if [[ "${worker_failed}" -eq 0 ]]; then
  "${V7_PY}" "${ROOT}/p0/select_shared_lr_v7_2_0.py" --run-root "${OUT}" || selector_status=$?
else
  echo "one or more LR pilot workers failed; outputs retained" >&2
  selector_status=4
fi
(
  cd "${OUT}"
  find . -type f ! -name 'LR_PILOT_MANIFEST.v7.2.0.sha256' -print0 | sort -z | xargs -0 sha256sum > LR_PILOT_MANIFEST.v7.2.0.sha256
)
echo "output_dir=${OUT}"
echo "selector_exit=${selector_status}"
exit "${selector_status}"

