#!/usr/bin/env bash

set +e
set +u
set +o pipefail

GRASS_ROOT=/data1/home/tanyuxiao/Grassmann_model
PANEL_ROOT="$GRASS_ROOT/v7/resources/panels/v7.1.0"
R4="$GRASS_ROOT/v7/code/releases/v7_20260825_p0_r4"
SITE_INPUT="$PANEL_ROOT/site_selection/v7.1.1"
SITE_OUT="$PANEL_ROOT/site_selection/v7.1.2"
BCF="$PANEL_ROOT/raw/hgdp1kgp_phased_v2/hgdp1kgp_chr22.filtered.SNV_INDEL.phased.shapeit5.bcf"
DATA_ENV="$GRASS_ROOT/v7/resources/envs/v7.1.0/conda-v7-data"
BCFTOOLS="$DATA_ENV/bin/bcftools"
BGZIP="$DATA_ENV/bin/bgzip"
TABIX="$DATA_ENV/bin/tabix"
DATA_PY="$DATA_ENV/bin/python"
DONORS="$SITE_INPUT/donor_all.samples.txt"
EXPECTED_SOURCE_SHA=09e50c698522d19bbc2c43a2ea2d29b290d63cf8c1d6983c5198dcb4289ff3fa

ALL_METRICS="$SITE_OUT/donor_FILTER_UNSET_biallelic_SNP.metrics.tsv.gz"
MAF_SITES="$SITE_OUT/donor_FILTER_UNSET_biallelic_SNP_MAFgt0p01.tsv.gz"
AUDIT="$SITE_OUT/SITE_STAGE1_AUDIT.v7.1.2.json"
LOG="$SITE_OUT/SITE_STAGE1.v7.1.2.log"

mkdir -p "$SITE_OUT"

actual_source_sha="$(sha256sum "$BCF" 2>/dev/null | awk '{print $1}')"
echo "expected_source_sha=$EXPECTED_SOURCE_SHA"
echo "actual_source_sha=$actual_source_sha"
if [ "$actual_source_sha" != "$EXPECTED_SOURCE_SHA" ]; then
    echo "STOP: source BCF hash mismatch"
    exit 10
fi
if [ "$(wc -l < "$DONORS")" -ne 2496 ]; then
    echo "STOP: donor sample count is not 2496"
    exit 11
fi
if [ -e "$ALL_METRICS" ] || [ -e "$MAF_SITES" ] || [ -e "$AUDIT" ]; then
    echo "STOP: v7.1.2 stage-1 output already exists"
    exit 12
fi

RUN_ID="$(date -u +%Y%m%dT%H%M%SZ)"
ALL_PART="$ALL_METRICS.part.$RUN_ID"
MAF_PART="$MAF_SITES.part.$RUN_ID"

"$BCFTOOLS" view \
  -S "$DONORS" \
  -f . \
  -m2 -M2 \
  -v snps \
  -Ou "$BCF" 2> "$LOG" \
| "$BCFTOOLS" +fill-tags -Ou -- \
  -t AC,AN,AF,MAF,F_MISSING 2>> "$LOG" \
| "$BCFTOOLS" query \
  -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/AF\t%INFO/MAF\t%INFO/F_MISSING\n' \
  2>> "$LOG" \
| "$BGZIP" -@ 2 -c > "$ALL_PART"

metric_codes=("${PIPESTATUS[@]}")
echo "metrics_pipeline_codes=${metric_codes[*]}"
metrics_ok=1
for code in "${metric_codes[@]}"; do
    if [ "$code" -ne 0 ]; then metrics_ok=0; fi
done
if [ "$metrics_ok" -ne 1 ] || [ "$("$BGZIP" -dc "$ALL_PART" 2>/dev/null | wc -l)" -eq 0 ]; then
    echo "STOP: metrics pipeline failed or was empty; part retained"
    exit 20
fi
mv "$ALL_PART" "$ALL_METRICS"

"$BGZIP" -dc "$ALL_METRICS" \
| awk -F '\t' '($8+0) > 0.01' \
| "$BGZIP" -@ 2 -c > "$MAF_PART"
maf_codes=("${PIPESTATUS[@]}")
echo "maf_pipeline_codes=${maf_codes[*]}"
maf_ok=1
for code in "${maf_codes[@]}"; do
    if [ "$code" -ne 0 ]; then maf_ok=0; fi
done
if [ "$maf_ok" -ne 1 ] || [ "$("$BGZIP" -dc "$MAF_PART" 2>/dev/null | wc -l)" -eq 0 ]; then
    echo "STOP: MAF pipeline failed or was empty; part retained"
    exit 21
fi
mv "$MAF_PART" "$MAF_SITES"

"$TABIX" -s 1 -b 2 -e 2 "$MAF_SITES"
tabix_exit=$?
echo "tabix_exit=$tabix_exit"
if [ "$tabix_exit" -ne 0 ]; then exit 22; fi

"$DATA_PY" "$R4/p0/audit_site_stage1_v7_1_2.py" \
  --all-metrics "$ALL_METRICS" \
  --maf-sites "$MAF_SITES" \
  --source-sha256 "$actual_source_sha" \
  --output "$AUDIT"
audit_exit=$?
echo "audit_exit=$audit_exit"
if [ "$audit_exit" -ne 0 ]; then exit 23; fi

sha256sum \
  "$DONORS" \
  "$ALL_METRICS" \
  "$MAF_SITES" \
  "$MAF_SITES.tbi" \
  "$AUDIT" \
  "$LOG" \
  > "$SITE_OUT/SITE_STAGE1.v7.1.2.sha256"

echo "stage1_status=PASS"
