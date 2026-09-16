#!/usr/bin/env bash
# Chain R1 inventory and validation in one read-only pass.
#
# Usage:
#   OUT=/mnt/bigdisk/e0_runs ./run_r1.sh --root /ukb/geno --root /ukb/kinship \
#       --cohort-id UKB --authorization-id 12345
#
# Everything is written under "$OUT"; nothing is written inside a scanned root.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="${OUT:?set OUT to a writable directory outside the cohort roots}"
RUN_ID="${RUN_ID:-r1_inventory_rc0}"
RUN_DIR="$OUT/$RUN_ID"

mkdir -p "$RUN_DIR"

echo "[1/2] read-only inventory -> $RUN_DIR"
python3 "$HERE/scripts/r1_inventory.py" --out-dir "$RUN_DIR" "$@" \
  > "$RUN_DIR/r1_inventory.stdout.log" 2> "$RUN_DIR/r1_inventory.stderr.log" || {
    echo "inventory failed; see $RUN_DIR/r1_inventory.stderr.log" >&2
    exit 1
  }
tail -n 20 "$RUN_DIR/r1_inventory.stdout.log"

echo "[2/2] independent validation"
if python3 "$HERE/scripts/r1_validate_inventory.py" \
      --inventory "$RUN_DIR/INVENTORY.json" \
      --output "$RUN_DIR/R1_VALIDATION.json"; then
  echo "R1 validation PASS; review $RUN_DIR/INVENTORY_SUMMARY.md before R2"
else
  echo "R1 validation FAIL; R2 stays blocked" >&2
  exit 1
fi
