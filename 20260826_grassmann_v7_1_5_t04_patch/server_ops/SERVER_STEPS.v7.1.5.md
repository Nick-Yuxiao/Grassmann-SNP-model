# v7.1.5 server steps

Deploy additively from r6 to r7. This patch computes a T04 contract from the
already frozen T02 and T03 evidence. It does not run a GPU workload.

```bash
set +e
set +u
set +o pipefail

GRASS_ROOT=/data1/home/tanyuxiao/Grassmann_model
R6="$GRASS_ROOT/v7/code/releases/v7_20260825_p0_r6"
R7="$GRASS_ROOT/v7/code/releases/v7_20260825_p0_r7"
BUNDLE="$GRASS_ROOT/incoming/v7_1_5_t04_patch_20260826.tar.gz"

if [ -e "$R7" ]; then
  echo "STOP: r7 already exists"
else
  cp -a "$R6" "$R7"
  tar -xzf "$BUNDLE" --strip-components=1 -C "$R7"
fi

python3 "$R7/validate_v7_1_5.py"
(
  cd "$R7" || exit 1
  sha256sum -c MANIFEST.v7.1.5.sha256
)
```

Build the CPU-only contract:

```bash
PROFILE="$GRASS_ROOT/v7/results/t03/v7.1.4/20260826T093031Z_t03_v7_1_gpu1_2451537/PROFILE_REPORT.v7.1.0.json"
MATERIAL="$GRASS_ROOT/v7/resources/panels/v7.1.0/materialized/v7.1.4/MATERIALIZATION_AUDIT.v7.1.4.json"
OUT="$GRASS_ROOT/v7/results/t04/v7.1.5"
DATA_PY="$GRASS_ROOT/v7/resources/envs/v7.1.0/conda-v7-data/bin/python"

"$DATA_PY" "$R7/p0/build_t04_contract_v7_1_5.py" \
  --grid "$R7/PILOT_GRID.v7.1.5.yaml" \
  --profile-report "$PROFILE" \
  --materialization-audit "$MATERIAL" \
  --output-dir "$OUT"
```
