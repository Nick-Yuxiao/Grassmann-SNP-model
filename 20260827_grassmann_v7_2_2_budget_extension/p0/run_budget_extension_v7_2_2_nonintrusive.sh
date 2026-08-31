#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
GRASS_ROOT="${V7_SERVER_ROOT:-/data1/home/tanyuxiao/Grassmann_model}"
V7_PY="${V7_PY:-${GRASS_ROOT}/v7/code/releases/v7_20260825_p0_r1/.conda-v7/bin/python}"
DATA_DIR="${V7_DATA_DIR:-${GRASS_ROOT}/v7/resources/panels/v7.1.0/a1r/preprocessed/v7.1.10}"
FIREWALL_DIR="${V7_LR_FIREWALL_DIR:-${GRASS_ROOT}/v7/resources/panels/v7.1.0/a1r/validation_firewall/v7.2.0}"
SOURCE_ROOT="${V7_BUDGET_BRIDGE_RUN_ROOT:?V7_BUDGET_BRIDGE_RUN_ROOT must name the immutable completed v7.2.1 bridge run}"
SCHEDULE="${V7_BUDGET_EXTENSION_SCHEDULE:-${ROOT}/BUDGET_EXTENSION_SCHEDULE.v7.2.2.tsv}"
OUT_PARENT="${V7_BUDGET_EXTENSION_OUT_PARENT:-${GRASS_ROOT}/v7/results/budget_extension/v7.2.2}"
GPUS=(1 3 4 5 6 7)

for path in "${V7_PY}" "${SCHEDULE}" \
  "${DATA_DIR}/PREPROCESS_MANIFEST.v7.1.10.sha256" \
  "${FIREWALL_DIR}/VALIDATION_FIREWALL.v7.2.0.sha256" \
  "${SOURCE_ROOT}/BUDGET_BRIDGE_MANIFEST.v7.2.1.sha256" \
  "${SOURCE_ROOT}/BUDGET_BRIDGE_DECISION.v7.2.1.json"; do
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
(
  cd "${SOURCE_ROOT}"
  sha256sum -c BUDGET_BRIDGE_MANIFEST.v7.2.1.sha256
)
python3 - "${SOURCE_ROOT}/BUDGET_BRIDGE_DECISION.v7.2.1.json" <<'PY'
import json,sys
d=json.load(open(sys.argv[1]))
if d.get('status')!='BUDGET_BRIDGE_EXTEND_ALL_TO_30K':
    raise SystemExit('source bridge does not authorize extension')
if d.get('result_count')!=12 or d.get('failure_count')!=0:
    raise SystemExit('source bridge result-count mismatch')
if d.get('decision_holdout_read') is not False:
    raise SystemExit('source bridge firewall violation')
if d.get('architecture_decision_permitted') is not False:
    raise SystemExit('source bridge architecture firewall violation')
PY

RUN_ID="$(date -u +%Y%m%dT%H%M%SZ)_budget_extension_v7_2_2_$$"
OUT="${OUT_PARENT}/${RUN_ID}"
mkdir "${OUT}"
AUDIT="${OUT}/GPU_AUDIT_BEFORE.json"
python3 "${ROOT}/p0/audit_server.py" --output "${AUDIT}"
python3 - "${AUDIT}" <<'PY'
import json,sys
p=json.load(open(sys.argv[1])); required={1,3,4,5,6,7}; idle=set(p.get('idle_gpu_indices',[]))
if not required<=idle:
    raise SystemExit(f'REFUSE: required GPUs not idle: {sorted(required-idle)}')
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
p=json.load(open(sys.argv[1])); gpu=int(sys.argv[2])
if gpu not in p.get('idle_gpu_indices',[]):
    raise SystemExit(f'REFUSE: GPU {gpu} became busy')
PY
    export CUDA_VISIBLE_DEVICES="${gpu}" PYTHONUNBUFFERED=1
    nice -n 10 "${V7_PY}" "${ROOT}/p0/run_budget_extension_gpu_worker_v7_2_2.py" \
      --schedule "${SCHEDULE}" --gpu "${gpu}" --data-dir "${DATA_DIR}" \
      --firewall-dir "${FIREWALL_DIR}" --source-root "${SOURCE_ROOT}" \
      --output-root "${OUT}" --trainer "${ROOT}/p0/train_budget_extension_v7_2_2.py"
  ) >"${OUT}/GPU_${gpu}_WORKER.log" 2>&1 &
  pids+=("$!")
done
for pid in "${pids[@]}"; do
  wait "${pid}" || worker_failed=1
done

assessment_status=0
"${V7_PY}" "${ROOT}/p0/assess_budget_extension_v7_2_2.py" \
  --run-root "${OUT}" --source-root "${SOURCE_ROOT}" || assessment_status=$?
if [[ "${worker_failed}" -ne 0 && "${assessment_status}" -eq 0 ]]; then
  echo "worker failure was not reflected by assessment" >&2
  assessment_status=4
fi
(
  cd "${OUT}"
  find . -type f ! -name 'BUDGET_EXTENSION_MANIFEST.v7.2.2.sha256' -print0 \
    | sort -z | xargs -0 sha256sum > BUDGET_EXTENSION_MANIFEST.v7.2.2.sha256
)
echo "output_dir=${OUT}"
echo "assessment_exit=${assessment_status}"
exit "${assessment_status}"
