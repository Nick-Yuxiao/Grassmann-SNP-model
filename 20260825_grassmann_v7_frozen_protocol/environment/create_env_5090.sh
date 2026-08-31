#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_DIR="${V7_ENV_DIR:-${ROOT}/.venv}"

python3 - <<'PY'
import sys
if sys.version_info < (3, 10):
    raise SystemExit(f"Python >=3.10 required, found {sys.version}")
PY

python3 -m venv "${ENV_DIR}"
source "${ENV_DIR}/bin/activate"
python -m pip install --upgrade pip setuptools wheel
python -m pip install "torch==2.7.1" --index-url https://download.pytorch.org/whl/cu128
python -m pip install "numpy==1.26.4" "scikit-learn==1.5.2" "pytest==8.3.5" "psutil==6.1.1"
python -m pip freeze | LC_ALL=C sort > "${ROOT}/environment/requirements-cu128.lock"
python "${ROOT}/environment/smoke_cuda.py" --output "${ROOT}/environment/ENV_SMOKE.json"
