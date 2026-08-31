#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
python3 "${ROOT}/environment/probe_server.py"
bash "${ROOT}/environment/create_env_5090.sh"
sha256sum \
  "${ROOT}/environment/requirements-cu128.lock" \
  "${ROOT}/environment/ENV_SMOKE.json" \
  "${ROOT}/environment/SERVER_RESOURCE.json" \
  > "${ROOT}/environment/T00_MANIFEST.sha256"
