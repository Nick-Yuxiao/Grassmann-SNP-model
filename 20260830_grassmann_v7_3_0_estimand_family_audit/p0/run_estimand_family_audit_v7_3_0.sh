#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
GRASS_ROOT="${V7_SERVER_ROOT:-/data1/home/tanyuxiao/Grassmann_model}"
V7_PY="${V7_PY:-${GRASS_ROOT}/v7/code/releases/v7_20260825_p0_r1/.conda-v7/bin/python}"
BRIDGE_20K_ROOT="${V7_BUDGET_BRIDGE_20K_RUN_ROOT:-${GRASS_ROOT}/v7/results/budget_bridge/v7.2.1/20260827T095917Z_budget_bridge_v7_2_1_3451087}"
EXTENSION_30K_ROOT="${V7_BUDGET_EXTENSION_30K_RUN_ROOT:-${GRASS_ROOT}/v7/results/budget_extension/v7.2.3/20260828T023841Z_budget_extension_v7_2_3_3758608}"
FINAL_40K_ROOT="${V7_FINAL_BUDGET_40K_RUN_ROOT:?V7_FINAL_BUDGET_40K_RUN_ROOT must name the immutable completed v7.2.4 run}"
C0_20K_ROOT="${V7_C0_20K_RUN_ROOT:-${GRASS_ROOT}/v7/results/c0/v7.1.13/20260827T060221Z_c0_extension_v7_1_13_3381139}"
OUT_PARENT="${V7_ESTIMAND_FAMILY_AUDIT_OUT_PARENT:-${GRASS_ROOT}/v7/results/estimand_family_audit/v7.3.0}"

for path in "${V7_PY}" "${ROOT}/p0/audit_estimand_family_v7_3_0.py"; do
  [[ -e "${path}" ]] || { echo "missing required path: ${path}" >&2; exit 8; }
done
for directory in "${BRIDGE_20K_ROOT}" "${EXTENSION_30K_ROOT}" "${FINAL_40K_ROOT}" "${C0_20K_ROOT}"; do
  [[ -d "${directory}" ]] || { echo "missing required input directory: ${directory}" >&2; exit 8; }
done

mkdir -p "${GRASS_ROOT}/v7/locks" "${OUT_PARENT}"
exec 9>"${GRASS_ROOT}/v7/locks/v7_3_0_estimand_family_cpu_audit.lock"
flock -n 9 || { echo "STOP: existing v7.3.0 CPU audit detected" >&2; exit 9; }

verify_manifest() {
  local directory="$1"
  local manifest="$2"
  [[ -f "${directory}/${manifest}" ]] || {
    echo "missing source manifest: ${directory}/${manifest}" >&2
    exit 8
  }
  (
    cd "${directory}"
    sha256sum --quiet -c "${manifest}"
  )
}

verify_manifest "${BRIDGE_20K_ROOT}" "BUDGET_BRIDGE_MANIFEST.v7.2.1.sha256"
verify_manifest "${EXTENSION_30K_ROOT}" "BUDGET_EXTENSION_MANIFEST.v7.2.3.sha256"
verify_manifest "${FINAL_40K_ROOT}" "FINAL_BUDGET_MANIFEST.v7.2.4.sha256"
verify_manifest "${C0_20K_ROOT}" "C0_EXTENSION_MANIFEST.v7.1.13.sha256"

RUN_ID="$(date -u +%Y%m%dT%H%M%SZ)_estimand_family_audit_v7_3_0_$$"
OUT="${OUT_PARENT}/${RUN_ID}"
mkdir "${OUT}"

"${V7_PY}" "${ROOT}/p0/audit_estimand_family_v7_3_0.py" \
  --bridge-20k-root "${BRIDGE_20K_ROOT}" \
  --extension-30k-root "${EXTENSION_30K_ROOT}" \
  --final-40k-root "${FINAL_40K_ROOT}" \
  --c0-20k-root "${C0_20K_ROOT}" \
  --output-dir "${OUT}"

(
  cd "${OUT}"
  find . -type f ! -name 'ESTIMAND_FAMILY_AUDIT_MANIFEST.v7.3.0.sha256' -print0 \
    | sort -z | xargs -0 sha256sum > ESTIMAND_FAMILY_AUDIT_MANIFEST.v7.3.0.sha256
)

echo "status=ESTIMAND_FAMILY_AUDIT_PASS_NO_GPU_AUTHORIZED"
echo "output_dir=${OUT}"
echo "gpu_used=false"
echo "next_authorized_stage=DRAFT_NEW_EXPERIMENTAL_CONTRACT_ONLY"
