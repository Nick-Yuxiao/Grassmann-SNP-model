# Server steps v7.1.7

Run each numbered block separately. These commands do not connect to or kill
existing jobs. Physical GPU 0 remains forbidden.

## 1. Deploy r9 from validated r8

```bash
set +e
set +u
set +o pipefail
GRASS_ROOT=/data1/home/tanyuxiao/Grassmann_model
R8="$GRASS_ROOT/v7/code/releases/v7_20260825_p0_r8"
R9="$GRASS_ROOT/v7/code/releases/v7_20260825_p0_r9"
BUNDLE="$GRASS_ROOT/incoming/v7_1_7_real_a1r_patch_20260826.tar.gz"
test -d "$R8" && test -f "$BUNDLE" && test ! -e "$R9"
echo "precheck_exit=$?"
```

```bash
cp -a "$R8" "$R9"
tar -xzf "$BUNDLE" --strip-components=1 -C "$R9"
find "$R9" -type f \( -name '*.py' -o -name '*.md' -o -name '*.tsv' -o -name '*.yaml' -o -name '*.sha256' \) -exec sed -i 's/\r$//' {} +
python3 "$R9/validate_v7_1_7.py"
```

```bash
cd "$R9"
sha256sum -c MANIFEST.v7.1.7.sha256
python3 -m unittest discover -s p0/tests -p 'test_v7_1_7.py' -v
```

## 2. Freeze nested real-donor subsets

```bash
GRASS_ROOT=/data1/home/tanyuxiao/Grassmann_model
R9="$GRASS_ROOT/v7/code/releases/v7_20260825_p0_r9"
FROZEN="$GRASS_ROOT/v7/resources/panels/v7.1.0/frozen/v7.1.1/1KGP_DONOR_TRAIN.tsv"
SUBSETS="$GRASS_ROOT/v7/resources/panels/v7.1.0/a1r/v7.1.7"
pgrep -a -u "$(id -u)" -f '[f]reeze_a1r_subsets_v7_1_7' || echo 'none detected'
```

```bash
python3 "$R9/p0/freeze_a1r_subsets_v7_1_7.py" --donor-train "$FROZEN" --output-dir "$SUBSETS"
subset_exit=$?
echo "subset_exit=$subset_exit"
```

```bash
cd "$SUBSETS"
sha256sum -c A1R_SUBSET_FREEZE.v7.1.7.sha256
wc -l DONOR_TRAIN_25P_562.samples.txt DONOR_TRAIN_50P_1124.samples.txt DONOR_TRAIN_100P_2247.samples.txt
```

## 3. Build the A1-R execution contract

```bash
PROFILE="$GRASS_ROOT/v7/results/t03/v7.1.4/20260826T093031Z_t03_v7_1_gpu1_2451537/PROFILE_REPORT.v7.1.0.json"
OUT="$GRASS_ROOT/v7/results/t04/v7.1.7"
python3 "$R9/p0/build_a1r_contract_v7_1_7.py" --grid "$R9/A1R_GRID.v7.1.7.yaml" --profile-report "$PROFILE" --subset-freeze "$SUBSETS/A1R_SUBSET_FREEZE.v7.1.7.json" --output-dir "$OUT"
echo "contract_exit=$?"
```

```bash
cd "$OUT"
sha256sum -c A1R_T04.v7.1.7.sha256
echo "scheduled_runs=$(($(wc -l < A1R_RUN_SCHEDULE.v7.1.7.tsv)-1))"
echo "pilot_runs=$(($(wc -l < A1R_C0_SCHEDULE.v7.1.7.tsv)-1))"
```

The training launcher must consume these schedules, perform a fresh read-only
GPU idle audit and obtain a non-blocking per-GPU lock before each block. It must
not access any file whose basename starts with `HGDP_`.
