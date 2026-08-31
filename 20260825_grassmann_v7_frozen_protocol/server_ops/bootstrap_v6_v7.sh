#!/usr/bin/env bash
set -euo pipefail
umask 077

ROOT="${1:-/data1/home/tanyuxiao/Grassmann_model}"
EXPECTED_ROOT="/data1/home/tanyuxiao/Grassmann_model"
if [[ "${ROOT}" != "${EXPECTED_ROOT}" ]]; then
  echo "Refusing unexpected root: ${ROOT}" >&2
  echo "Expected exactly: ${EXPECTED_ROOT}" >&2
  exit 2
fi

for version in v6 v7; do
  mkdir -p \
    "${ROOT}/${version}/code/releases" \
    "${ROOT}/${version}/results" \
    "${ROOT}/${version}/resources/audits" \
    "${ROOT}/${version}/resources/gpu_test" \
    "${ROOT}/${version}/logs" \
    "${ROOT}/${version}/locks" \
    "${ROOT}/${version}/runbooks"
done
mkdir -p "${ROOT}/incoming"

RUN_ID="$(date -u +%Y%m%dT%H%M%SZ)_bootstrap_$$"
SNAPSHOT="${ROOT}/v7/resources/audits/${RUN_ID}"
mkdir "${SNAPSHOT}"

# Read-only inventory. Failures are retained but do not erase successful probes.
{
  date -u +%FT%TZ
  uname -a
  id
  uptime
} > "${SNAPSHOT}/host.txt" 2>&1 || true

{
  nvidia-smi --query-gpu=index,uuid,name,memory.total,memory.used,memory.free,utilization.gpu,temperature.gpu,pstate,driver_version \
    --format=csv,noheader,nounits
} > "${SNAPSHOT}/gpus.csv" 2>&1 || true

{
  nvidia-smi --query-compute-apps=gpu_uuid,pid,process_name,used_memory \
    --format=csv,noheader,nounits
} > "${SNAPSHOT}/gpu_compute_processes.csv" 2>&1 || true

ps -eo pid,ppid,user,stat,etimes,pcpu,pmem,comm --sort=-pcpu \
  > "${SNAPSHOT}/processes.txt" 2>&1 || true
df -B1 "${ROOT}" > "${SNAPSHOT}/disk_bytes.txt" 2>&1 || true
df -h "${ROOT}" > "${SNAPSHOT}/disk_human.txt" 2>&1 || true
free -b > "${SNAPSHOT}/memory_bytes.txt" 2>&1 || true
lscpu > "${SNAPSHOT}/cpu.txt" 2>&1 || true

if command -v squeue >/dev/null 2>&1; then
  squeue -u "$(id -un)" -o '%.18i %.9P %.32j %.8T %.10M %.6D %.30R' \
    > "${SNAPSHOT}/slurm_queue.txt" 2>&1 || true
else
  printf '%s\n' 'SLURM_NOT_AVAILABLE' > "${SNAPSHOT}/slurm_queue.txt"
fi

find "${ROOT}" -maxdepth 3 -type d -printf '%M\t%u\t%g\t%p\n' \
  > "${SNAPSHOT}/directory_tree.txt"
find "${SNAPSHOT}" -maxdepth 1 -type f ! -name 'MANIFEST.sha256' -print0 \
  | sort -z | xargs -0 sha256sum > "${SNAPSHOT}/MANIFEST.sha256"

printf '%s\n' \
  "Bootstrap complete: ${ROOT}" \
  "Resource snapshot: ${SNAPSHOT}" \
  "No process was signalled, stopped, reprioritized, or removed."
