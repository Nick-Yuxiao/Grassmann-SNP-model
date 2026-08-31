# Server steps for v7.1.0

These commands do not connect to the server from the local machine automatically,
do not signal or kill any process, and refuse GPU 0.

If a JSON view is stuck at `(END)`, press `q` once. It is the `less` pager, not a
hung Python process.

## 1. Optional PowerShell upload

Run locally after replacing the SSH alias:

```powershell
$Bundle = Join-Path (Get-Location) 'deliverables\v7_1_0_protocol_patch_20260825.tar.gz'
Get-FileHash -Algorithm SHA256 -LiteralPath $Bundle
$SshTarget = '<YOUR_SSH_ALIAS>'
scp -- $Bundle "${SshTarget}:/data1/home/tanyuxiao/Grassmann_model/incoming/"
```

Uploading through the VS remote file explorer is also fine. Do not extract into the
existing `v7_20260825_p0_r1` release.

## 2. Create release r2 from the original base bundle plus this additive patch

Paste the complete block into one server shell:

```bash
set -euo pipefail
export GRASS_ROOT=/data1/home/tanyuxiao/Grassmann_model
export OLD_RELEASE_DIR="$GRASS_ROOT/v7/code/releases/v7_20260825_p0_r1"
export NEW_RELEASE_ID=v7_20260825_p0_r2
export NEW_RELEASE_DIR="$GRASS_ROOT/v7/code/releases/$NEW_RELEASE_ID"
export BASE_BUNDLE="$GRASS_ROOT/incoming/v7_p0_server_bundle_20260825.tar.gz"
export PATCH_BUNDLE="$GRASS_ROOT/incoming/v7_1_0_protocol_patch_20260825.tar.gz"

test -d "$OLD_RELEASE_DIR"
test -x "$OLD_RELEASE_DIR/.conda-v7/bin/python"
test -f "$BASE_BUNDLE"
test -f "$PATCH_BUNDLE"
if [ -e "$NEW_RELEASE_DIR" ]; then
  echo "STOP: r2 already exists: $NEW_RELEASE_DIR"
  return 1 2>/dev/null || exit 1
fi
mkdir "$NEW_RELEASE_DIR"
tar -xzf "$BASE_BUNDLE" --strip-components=1 -C "$NEW_RELEASE_DIR"
tar -xzf "$PATCH_BUNDLE" --strip-components=1 -C "$NEW_RELEASE_DIR"
cd "$NEW_RELEASE_DIR"
python3 validate_v7_1.py
sha256sum -c MANIFEST.v7.1.0.sha256
```

Use the already tested isolated environment from r1 without changing it:

```bash
export V7_PY="$OLD_RELEASE_DIR/.conda-v7/bin/python"
export V7_SERVER_ROOT="$GRASS_ROOT"
export V7_GPU_INDEX=1
export CUDA_VISIBLE_DEVICES=1
"$V7_PY" -c 'import torch; print(torch.__version__, torch.version.cuda, torch.cuda.is_available(), torch.cuda.get_device_name(0))'
```

The last device is process-local `0`; because `CUDA_VISIBLE_DEVICES=1`, it maps to
physical GPU 1.

## 3. Finish T01 before profiling

Prepare immutable files under a new data-manifest directory. Required inputs are:

- normalized chr22 joint panel;
- final ordered variant-ID list;
- disjoint donor-train and donor-validation sample lists;
- disjoint HGDP primary and exact 216-person calibration lists;
- `sample_id`, `population`, `cohort` TSV;
- frozen HAPNEST config;
- frozen donor-PC neighbour classification;
- public source-release/checksum record.

Then run, replacing only the paths in angle brackets:

```bash
cd "$NEW_RELEASE_DIR"
T01_OUT="$GRASS_ROOT/v7/results/T01_v7_1_0_$(date -u +%Y%m%dT%H%M%SZ)"
"$V7_PY" p0/build_panel_manifest_v7_1.py \
  --panel-spec PANEL_SPEC.v7.1.0.yaml \
  --joint-chr22-panel '<CHR22_PANEL.bcf>' \
  --variants '<CHR22_VARIANTS.tsv>' \
  --donor-train-samples '<1KGP_DONOR_TRAIN.tsv>' \
  --donor-validation-samples '<1KGP_DONOR_VALIDATION.tsv>' \
  --hgdp-primary-samples '<HGDP_PRIMARY.tsv>' \
  --hgdp-snpbag-calibration-samples '<HGDP_SNPBAG_216.tsv>' \
  --sample-populations '<SAMPLE_POPULATION_COHORT.tsv>' \
  --hapnest-config '<HAPNEST_FROZEN_CONFIG.yaml>' \
  --neighbor-classification '<HGDP_1KGP_NEIGHBOR_CLASSIFICATION.tsv>' \
  --source-release-record '<SOURCE_RELEASES.sha256>' \
  --compatibility-note 'exact release/filter comparison recorded here' \
  --output-dir "$T01_OUT"
"$V7_PY" -m json.tool "$T01_OUT/PANEL_MANIFEST.v7.1.0.json" | less
```

Press `q` to leave `less`. Do not start T03 if this builder reports overlap,
incorrect cohort membership, a non-216 calibration set or missing hashes.

## 4. Run T03 on physical GPU 1 without interrupting other work

```bash
export V7_CHR22_L="$("$V7_PY" -c 'import json,sys; print(json.load(open(sys.argv[1]))["chr22"]["variant_count"])' "$T01_OUT/PANEL_MANIFEST.v7.1.0.json")"
echo "data-derived chr22 L=$V7_CHR22_L"
bash p0/run_t03_gpu1_v7_1_nonintrusive.sh
```

The script performs two read-only audits, takes a nonblocking project lock, refuses
physical GPU 0 and refuses an existing output directory. It never invokes `kill`,
`pkill`, preemption or scheduler cancellation. If GPU 1 is busy, stop and wait.

From a second terminal, read-only monitoring is:

```bash
nvidia-smi -i 1
pgrep -a -u "$(id -u)" -f 'profile_models_v7_1|run_t03_gpu1_v7_1' || true
```

## 5. Build T04 only after signing the training step budget

Locate the new report and replace the signed capacity values:

```bash
PROFILE_REPORT="$(find "$NEW_RELEASE_DIR/p0/profile_runtime_v7_1" -name 'PROFILE_REPORT.v7.1.0.json' -type f | sort | tail -n 1)"
T04_OUT="$GRASS_ROOT/v7/results/T04_v7_1_0_$(date -u +%Y%m%dT%H%M%SZ)"
"$V7_PY" p0/build_compute_contract_v7_1.py \
  --grid GRID_SPEC.v7.1.0.yaml \
  --panel-manifest "$T01_OUT/PANEL_MANIFEST.v7.1.0.json" \
  --profile-report "$PROFILE_REPORT" \
  --train-steps-per-run '<SIGNED_INTEGER>' \
  --training-budget-reference '<SIGNED_RECORD_ID>' \
  --signed-gpu-hours '<SIGNED_GPU_HOURS>' \
  --signed-storage-gib '<SIGNED_STORAGE_GIB>' \
  --gpu-count '<AUTHORIZED_NONZERO_GPU_COUNT>' \
  --concurrency-limit '<NONINTERRUPTING_CONCURRENCY>' \
  --available-window-hours '<SIGNED_WINDOW_HOURS>' \
  --output-dir "$T04_OUT"
```

Do not invent the step budget from favorable profiler or pilot outcomes. A capacity
failure is `P0_CAPACITY_NO_GO`, not a scientific architecture NO-GO.

## 6. P0 readiness verdict

```bash
VERDICT_OUT="$GRASS_ROOT/v7/results/P0_v7_1_0_$(date -u +%Y%m%dT%H%M%SZ)"
"$V7_PY" p0/assess_p0_v7_1.py \
  --t00-smoke "$OLD_RELEASE_DIR/environment/ENV_SMOKE.json" \
  --t00-resource "$OLD_RELEASE_DIR/environment/SERVER_RESOURCE.json" \
  --t00-runtime-manifest "$OLD_RELEASE_DIR/p0/T00_RUNTIME_MANIFEST.sha256" \
  --panel-manifest "$T01_OUT/PANEL_MANIFEST.v7.1.0.json" \
  --profile-report "$PROFILE_REPORT" \
  --compute-contract "$T04_OUT/COMPUTE_CONTRACT.v7.1.0.json" \
  --output-dir "$VERDICT_OUT"
"$V7_PY" -m json.tool "$VERDICT_OUT/P0_VERDICT.v7.1.0.json"
```

Only `GO_TO_A0` permits implementation/integrity work. It is not the scientific
Grassmann GO, which still requires the preregistered A1 results.
