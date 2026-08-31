# v7.1.6 server steps

Deploy additively from r7 to r8. This patch builds the CPU-only primary-first
capacity contract and a balanced 120-run schedule; it does not launch GPUs.

```bash
set +e
set +u
set +o pipefail

GRASS_ROOT=/data1/home/tanyuxiao/Grassmann_model
R7="$GRASS_ROOT/v7/code/releases/v7_20260825_p0_r7"
R8="$GRASS_ROOT/v7/code/releases/v7_20260825_p0_r8"
BUNDLE="$GRASS_ROOT/incoming/v7_1_6_primary_first_patch_20260826.tar.gz"

if [ -e "$R8" ]; then
  echo "STOP: r8 already exists"
else
  cp -a "$R7" "$R8"
  tar -xzf "$BUNDLE" --strip-components=1 -C "$R8"
fi

python3 "$R8/validate_v7_1_6.py"
(
  cd "$R8" || exit 1
  sha256sum -c MANIFEST.v7.1.6.sha256
)
```

Build the contract and schedule:

```bash
PROFILE="$GRASS_ROOT/v7/results/t03/v7.1.4/20260826T093031Z_t03_v7_1_gpu1_2451537/PROFILE_REPORT.v7.1.0.json"
MATERIAL="$GRASS_ROOT/v7/resources/panels/v7.1.0/materialized/v7.1.4/MATERIALIZATION_AUDIT.v7.1.4.json"
OUT="$GRASS_ROOT/v7/results/t04/v7.1.6"
DATA_PY="$GRASS_ROOT/v7/resources/envs/v7.1.0/conda-v7-data/bin/python"

"$DATA_PY" "$R8/p0/build_t04_primary_first_v7_1_6.py" \
  --grid "$R8/PRIMARY_FIRST_GRID.v7.1.6.yaml" \
  --profile-report "$PROFILE" \
  --materialization-audit "$MATERIAL" \
  --output-dir "$OUT"
```
