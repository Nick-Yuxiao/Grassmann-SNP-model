# Server steps for v7.3.0

_CPU-only deployment and audit steps; no command in this document starts or queries a GPU task._

---

## 📋 Preconditions

- Treat release r22 and all completed v7.2.1, v7.2.3 and v7.2.4 result directories as immutable.
- Copy r22 to a new release r23, then overlay this patch with one stripped top-level component.
- Verify the patch SHA-256 supplied with the transfer, release manifest, validator, tests and shell syntax.
- Do not run this audit while another v7.3.0 CPU audit holds the advisory lock.
- This audit does not require an idle GPU and must not inspect, signal or terminate any GPU process.

## 🚀 Deploy release r23

Run this block from any server directory. Replace only `PATCH_EXPECTED_SHA256` with the checksum supplied beside the patch.

```bash
set +e
set +u
set +o pipefail

GRASS_ROOT=/data1/home/tanyuxiao/Grassmann_model
R22="$GRASS_ROOT/v7/code/releases/v7_20260825_p0_r22"
R23="$GRASS_ROOT/v7/code/releases/v7_20260825_p0_r23"
PATCH="$GRASS_ROOT/incoming/v7_3_0_estimand_family_audit_patch_20260830.tar.gz"
PATCH_EXPECTED_SHA256=REPLACE_WITH_SUPPLIED_SHA256

printf '%s  %s\n' "$PATCH_EXPECTED_SHA256" "$PATCH" | sha256sum -c -
echo "patch_hash_exit=$?"

if [ -e "$R23" ]; then
  echo "STOP: release already exists: $R23"
else
  cp -a "$R22" "$R23"
  tar -xzf "$PATCH" --strip-components=1 -C "$R23"
  echo "deploy_exit=$?"
fi
```

## ✅ Validate release r23

```bash
GRASS_ROOT=/data1/home/tanyuxiao/Grassmann_model
R23="$GRASS_ROOT/v7/code/releases/v7_20260825_p0_r23"
V7_PY="$GRASS_ROOT/v7/code/releases/v7_20260825_p0_r1/.conda-v7/bin/python"

cd "$R23" || exit
"$V7_PY" validate_v7_3_0.py
echo "validator_exit=$?"

sha256sum -c MANIFEST.v7.3.0.sha256
echo "manifest_exit=$?"

"$V7_PY" -m unittest discover -s p0/tests -p 'test_v7_3_0.py' -v
echo "test_exit=$?"

bash -n p0/run_estimand_family_audit_v7_3_0.sh
echo "shell_syntax_exit=$?"
```

## 🔎 Locate immutable 40k input

```bash
GRASS_ROOT=/data1/home/tanyuxiao/Grassmann_model

FINAL_40K_ROOT="$(find "$GRASS_ROOT/v7/results/final_budget/v7.2.4" \
  -mindepth 1 -maxdepth 1 -type d \
  -name '*_final_budget_v7_2_4_*' | sort | tail -n 1)"

echo "final_40k_root=$FINAL_40K_ROOT"
test -f "$FINAL_40K_ROOT/FINAL_BUDGET_DECISION.v7.2.4.json"
echo "final_decision_exists=$?"
```

## 🧪 Run the CPU-only audit

This command normally finishes quickly and may run in the foreground. It verifies the 20k/30k/40k lineage manifests and the historical C0 manifest before reading curves.

```bash
GRASS_ROOT=/data1/home/tanyuxiao/Grassmann_model
R23="$GRASS_ROOT/v7/code/releases/v7_20260825_p0_r23"
V7_PY="$GRASS_ROOT/v7/code/releases/v7_20260825_p0_r1/.conda-v7/bin/python"

FINAL_40K_ROOT="$(find "$GRASS_ROOT/v7/results/final_budget/v7.2.4" \
  -mindepth 1 -maxdepth 1 -type d \
  -name '*_final_budget_v7_2_4_*' | sort | tail -n 1)"

env \
  V7_SERVER_ROOT="$GRASS_ROOT" \
  V7_PY="$V7_PY" \
  V7_FINAL_BUDGET_40K_RUN_ROOT="$FINAL_40K_ROOT" \
  bash "$R23/p0/run_estimand_family_audit_v7_3_0.sh"

echo "audit_exit=$?"
```

## 📊 Verify and summarize the output

```bash
GRASS_ROOT=/data1/home/tanyuxiao/Grassmann_model
V7_PY="$GRASS_ROOT/v7/code/releases/v7_20260825_p0_r1/.conda-v7/bin/python"

AUDIT_ROOT="$(find "$GRASS_ROOT/v7/results/estimand_family_audit/v7.3.0" \
  -mindepth 1 -maxdepth 1 -type d \
  -name '*_estimand_family_audit_v7_3_0_*' | sort | tail -n 1)"

echo "audit_root=$AUDIT_ROOT"
(
  cd "$AUDIT_ROOT" || exit
  sha256sum -c ESTIMAND_FAMILY_AUDIT_MANIFEST.v7.3.0.sha256
)
echo "audit_manifest_exit=$?"

"$V7_PY" -c 'import json,sys; d=json.load(open(sys.argv[1])); print("status:",d["status"]); print("lineage:",d["lineage_20k_30k_40k"]); print("false_plateau_cells:",d["false_plateau_cells_30k_pass_40k_fail"]); print("ld_capacity:",d["families"]["ld_block_0p90"]["capacity_status"]); print("longrange_capacity:",d["families"]["within_chrom_longrange_0p90"]["capacity_status"]); print("variance_role:",d["historical_variance_proxy"]["role"]); print("gpu_used:",d["gpu_used"]); print("next:",d["next_authorized_stage"])' \
  "$AUDIT_ROOT/ESTIMAND_FAMILY_AUDIT.v7.3.0.json"
```

Expected terminal status: `ESTIMAND_FAMILY_AUDIT_PASS_NO_GPU_AUTHORIZED`.
The only authorized next stage is drafting a new experimental contract. Do not start an n=5 pilot, formal A1-R, HAPNEST or another horizon extension from this status.
