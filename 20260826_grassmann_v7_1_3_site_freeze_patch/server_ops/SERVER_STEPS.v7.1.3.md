# v7.1.3 server steps

Deploy additively from r4 to r5. The finalizer reads existing Stage-1 and Stage-2
metrics, rejects duplicate positions, requires complete HGDP GT, and writes a new
immutable final-site directory. It does not use a GPU.

```bash
set +e
set +u
set +o pipefail

GRASS_ROOT=/data1/home/tanyuxiao/Grassmann_model
R4="$GRASS_ROOT/v7/code/releases/v7_20260825_p0_r4"
R5="$GRASS_ROOT/v7/code/releases/v7_20260825_p0_r5"
BUNDLE="$GRASS_ROOT/incoming/v7_1_3_site_freeze_patch_20260826.tar.gz"

if [ -e "$R5" ]; then
  echo "STOP: r5 already exists"
else
  cp -a "$R4" "$R5"
  tar -xzf "$BUNDLE" --strip-components=1 -C "$R5"
fi

python3 "$R5/validate_v7_1_3.py"
(
  cd "$R5" || return
  sha256sum -c MANIFEST.v7.1.3.sha256
)
```

After validation:

```bash
PANEL_ROOT="$GRASS_ROOT/v7/resources/panels/v7.1.0"
DATA_PY="$GRASS_ROOT/v7/resources/envs/v7.1.0/conda-v7-data/bin/python"
"$DATA_PY" "$R5/p0/finalize_sites_v7_1_3.py" \
  --donor-sites "$PANEL_ROOT/site_selection/v7.1.2/donor_FILTER_UNSET_biallelic_SNP_MAFgt0p01.tsv.gz" \
  --hgdp-metrics "$PANEL_ROOT/site_selection/v7.1.2/hgdp_FILTER_UNSET_biallelic_SNP.metrics.tsv.gz" \
  --output-dir "$PANEL_ROOT/site_selection/v7.1.3"
```
