# v7.1.4 server steps

Deploy additively from r5 to r6. The runner is CPU/I/O-only, holds a file lock,
does not overwrite an existing final output, and never signals another process.

```bash
set +e
set +u
set +o pipefail

GRASS_ROOT=/data1/home/tanyuxiao/Grassmann_model
R5="$GRASS_ROOT/v7/code/releases/v7_20260825_p0_r5"
R6="$GRASS_ROOT/v7/code/releases/v7_20260825_p0_r6"
BUNDLE="$GRASS_ROOT/incoming/v7_1_4_materialization_patch_20260826.tar.gz"

if [ -e "$R6" ]; then
  echo "STOP: r6 already exists"
else
  cp -a "$R5" "$R6"
  tar -xzf "$BUNDLE" --strip-components=1 -C "$R6"
fi

python3 "$R6/validate_v7_1_4.py"
(
  cd "$R6" || exit 1
  sha256sum -c MANIFEST.v7.1.4.sha256
)
bash -n "$R6/p0/materialize_panel_v7_1_4.sh"
```

Run materialization in the foreground:

```bash
bash "$R6/p0/materialize_panel_v7_1_4.sh"
```

If the final status is PASS, verify the immutable result manifest:

```bash
OUT=/data1/home/tanyuxiao/Grassmann_model/v7/resources/panels/v7.1.0/materialized/v7.1.4
(
  cd "$OUT" || exit 1
  sha256sum -c MATERIALIZATION.v7.1.4.sha256
)
```
