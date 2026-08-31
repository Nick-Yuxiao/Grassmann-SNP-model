# Server steps for v7.1.1 panel correction

These steps create an additive r3 release. They do not modify r2, use a GPU,
kill a process or start a background task.

## Deploy

```bash
set +e
set +u
set +o pipefail

GRASS_ROOT=/data1/home/tanyuxiao/Grassmann_model
R2="$GRASS_ROOT/v7/code/releases/v7_20260825_p0_r2"
R3="$GRASS_ROOT/v7/code/releases/v7_20260825_p0_r3"
BUNDLE="$GRASS_ROOT/incoming/v7_1_1_panel_patch_20260826.tar.gz"

if [ -e "$R3" ]; then
  echo "STOP: r3 already exists: $R3"
else
  mkdir -p "$R3"
  cp -a "$R2/." "$R3/"
  tar -xzf "$BUNDLE" --strip-components=1 -C "$R3"
fi

python3 "$R3/validate_v7_1_1.py"
(
  cd "$R3" || return
  tr -d '\r' < MANIFEST.v7.1.1.sha256 | sha256sum -c -
)
```

## Freeze samples

```bash
set +e
set +u
set +o pipefail

GRASS_ROOT=/data1/home/tanyuxiao/Grassmann_model
R3="$GRASS_ROOT/v7/code/releases/v7_20260825_p0_r3"
PANEL_ROOT="$GRASS_ROOT/v7/resources/panels/v7.1.0"
DATA_PY="$GRASS_ROOT/v7/resources/envs/v7.1.0/conda-v7-data/bin/python"
META="$PANEL_ROOT/source_catalog/gnomad_meta_updated.tsv"
SAMPLES="$(find "$PANEL_ROOT/audits" -type f -name bcf_samples.txt | sort | tail -n 1)"
OUT="$PANEL_ROOT/frozen/v7.1.1"

"$DATA_PY" "$R3/p0/freeze_samples_v7_1_1.py" \
  --metadata "$META" \
  --bcf-samples "$SAMPLES" \
  --output-dir "$OUT"

freeze_exit=$?
echo "freeze_exit=$freeze_exit"

if [ "$freeze_exit" -eq 0 ]; then
  "$DATA_PY" -m json.tool "$OUT/SPLIT_SUMMARY.json"
  (
    cd "$OUT" || return
    sha256sum -c SAMPLE_FREEZE.sha256
  )
fi
```
