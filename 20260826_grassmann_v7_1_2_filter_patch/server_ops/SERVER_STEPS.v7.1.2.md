# v7.1.2 server steps

Deploy additively from r3 to r4, then execute the stage-1 script with `bash`.
The script is a child process; its explicit exit codes cannot close the
interactive terminal. It uses no GPU and starts no background task.

```bash
set +e
set +u
set +o pipefail

GRASS_ROOT=/data1/home/tanyuxiao/Grassmann_model
R3="$GRASS_ROOT/v7/code/releases/v7_20260825_p0_r3"
R4="$GRASS_ROOT/v7/code/releases/v7_20260825_p0_r4"
BUNDLE="$GRASS_ROOT/incoming/v7_1_2_filter_patch_20260826.tar.gz"

if [ -e "$R4" ]; then
  echo "STOP: r4 already exists"
else
  cp -a "$R3" "$R4"
  tar -xzf "$BUNDLE" --strip-components=1 -C "$R4"
fi

python3 "$R4/validate_v7_1_2.py"
(
  cd "$R4" || return
  sha256sum -c MANIFEST.v7.1.2.sha256
)
bash -n "$R4/p0/run_site_stage1_v7_1_2.sh"
```

After all three checks pass:

```bash
bash "$R4/p0/run_site_stage1_v7_1_2.sh"
echo "site_stage1_exit=$?"
```
