#!/usr/bin/env bash

set +e
set +u
set +o pipefail

GRASS_ROOT=/data1/home/tanyuxiao/Grassmann_model
PANEL_ROOT="$GRASS_ROOT/v7/resources/panels/v7.1.0"
R6="$GRASS_ROOT/v7/code/releases/v7_20260825_p0_r6"
DATA_ENV="$GRASS_ROOT/v7/resources/envs/v7.1.0/conda-v7-data"
BCFTOOLS="$DATA_ENV/bin/bcftools"
BGZIP="$DATA_ENV/bin/bgzip"
TABIX="$DATA_ENV/bin/tabix"
DATA_PY="$DATA_ENV/bin/python"

SOURCE_BCF="$PANEL_ROOT/raw/hgdp1kgp_phased_v2/hgdp1kgp_chr22.filtered.SNV_INDEL.phased.shapeit5.bcf"
FINAL_DIR="$PANEL_ROOT/site_selection/v7.1.3"
FINAL_VARIANTS="$FINAL_DIR/FINAL_VARIANTS.v7.1.3.tsv"
FINAL_IDS="$FINAL_DIR/FINAL_VARIANT_IDS.v7.1.3.txt"
FROZEN="$PANEL_ROOT/frozen/v7.1.1"
SITE_LISTS="$PANEL_ROOT/site_selection/v7.1.1"
OUT_PARENT="$PANEL_ROOT/materialized"
OUT_FINAL="$OUT_PARENT/v7.1.4"
LOCK_DIR="$PANEL_ROOT/locks"
LOCK_FILE="$LOCK_DIR/materialize_panel_v7_1_4.lock"
EXPECTED_SOURCE_SHA=09e50c698522d19bbc2c43a2ea2d29b290d63cf8c1d6983c5198dcb4289ff3fa

stop_run() {
  echo "STOP: $1"
  if [ -n "${WORK_DIR:-}" ]; then
    echo "retained_work_dir=$WORK_DIR"
  fi
  exit 1
}

require_file() {
  if [ ! -f "$1" ]; then
    stop_run "required file missing: $1"
  fi
}

echo "=== V7.1.4 PRECHECK ==="
for required in "$BCFTOOLS" "$BGZIP" "$TABIX" "$DATA_PY" "$SOURCE_BCF" "$SOURCE_BCF.csi" \
  "$FINAL_VARIANTS" "$FINAL_IDS" "$FINAL_DIR/FINAL_SITE_FREEZE.v7.1.3.sha256" \
  "$FROZEN/SAMPLE_POPULATION_COHORT.tsv" \
  "$SITE_LISTS/donor_train.samples.txt" "$SITE_LISTS/donor_validation.samples.txt" \
  "$SITE_LISTS/hgdp_primary.samples.txt" \
  "$R6/p0/build_materialization_audit_v7_1_4.py"; do
  require_file "$required"
done

mkdir -p "$OUT_PARENT" "$LOCK_DIR"
mkdir_exit=$?
if [ "$mkdir_exit" -ne 0 ]; then
  stop_run "could not create output/lock parent"
fi

if ! command -v flock >/dev/null 2>&1; then
  stop_run "flock is unavailable"
fi
exec 9>"$LOCK_FILE"
flock -n 9
lock_exit=$?
if [ "$lock_exit" -ne 0 ]; then
  stop_run "another v7.1.4 materialization holds the lock"
fi

if [ -e "$OUT_FINAL" ]; then
  stop_run "final output already exists: $OUT_FINAL"
fi

RUNNING="$(pgrep -a -u "$(id -u)" -f '[b]cftools.*hgdp1kgp_chr22' 2>/dev/null)"
if [ -n "$RUNNING" ]; then
  echo "$RUNNING"
  stop_run "an existing source-panel bcftools task was detected; nothing was interrupted"
fi

(
  cd "$FINAL_DIR" || exit 1
  sha256sum -c FINAL_SITE_FREEZE.v7.1.3.sha256
)
freeze_manifest_exit=$?
if [ "$freeze_manifest_exit" -ne 0 ]; then
  stop_run "v7.1.3 site-freeze manifest failed"
fi

actual_source_sha="$(sha256sum "$SOURCE_BCF" | awk '{print $1}')"
echo "expected_source_sha=$EXPECTED_SOURCE_SHA"
echo "actual_source_sha=$actual_source_sha"
if [ "$actual_source_sha" != "$EXPECTED_SOURCE_SHA" ]; then
  stop_run "source BCF SHA-256 mismatch"
fi

final_rows="$(awk 'NR>1{n++} END{print n+0}' "$FINAL_VARIANTS")"
id_rows="$(awk 'NF{n++} END{print n+0}' "$FINAL_IDS")"
if [ "$final_rows" -ne 154850 ] || [ "$id_rows" -ne 154850 ]; then
  stop_run "frozen site count is not 154850"
fi

WORK_DIR="$(mktemp -d "$OUT_PARENT/.v7.1.4_build.XXXXXX")"
mktemp_exit=$?
if [ "$mktemp_exit" -ne 0 ] || [ -z "$WORK_DIR" ] || [ ! -d "$WORK_DIR" ]; then
  stop_run "could not create unique build directory"
fi
echo "work_dir=$WORK_DIR"
LOG="$WORK_DIR/MATERIALIZATION.v7.1.4.log"

echo "=== BUILD EXACT SITE VCF ==="
{
  printf '##fileformat=VCFv4.2\n'
  printf '##contig=<ID=chr22>\n'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
  awk -F '\t' 'NR>1 {print $2"\t"$3"\t"$1"\t"$4"\t"$5"\t.\t.\t."}' "$FINAL_VARIANTS"
} | "$BGZIP" -c > "$WORK_DIR/FINAL_SITES.v7.1.4.vcf.gz"
site_vcf_codes=("${PIPESTATUS[@]}")
echo "site_vcf_codes=${site_vcf_codes[*]}" | tee -a "$LOG"
for code in "${site_vcf_codes[@]}"; do
  if [ "$code" -ne 0 ]; then
    stop_run "exact-site VCF construction failed"
  fi
done
"$TABIX" -f -p vcf "$WORK_DIR/FINAL_SITES.v7.1.4.vcf.gz" 2>>"$LOG"
if [ "$?" -ne 0 ]; then
  stop_run "exact-site VCF indexing failed"
fi

echo "=== EXACT INTERSECTION AND GENOTYPE-ONLY CLEANUP ==="
"$BCFTOOLS" isec -c none -n=2 -w1 -Ou "$SOURCE_BCF" "$WORK_DIR/FINAL_SITES.v7.1.4.vcf.gz" 2>>"$LOG" | \
  "$BCFTOOLS" annotate -x INFO,FORMAT -Ob -o "$WORK_DIR/ALL_SOURCE_EXACT_4091.bcf" 2>>"$LOG"
isec_codes=("${PIPESTATUS[@]}")
echo "isec_pipeline_codes=${isec_codes[*]}" | tee -a "$LOG"
for code in "${isec_codes[@]}"; do
  if [ "$code" -ne 0 ]; then
    stop_run "exact intersection or PP cleanup failed"
  fi
done
"$BCFTOOLS" index -f "$WORK_DIR/ALL_SOURCE_EXACT_4091.bcf" 2>>"$LOG"
if [ "$?" -ne 0 ]; then
  stop_run "all-source BCF indexing failed"
fi

"$BCFTOOLS" query -l "$SOURCE_BCF" > "$WORK_DIR/EXPECTED_SAMPLES.all_source.txt" 2>>"$LOG"
if [ "$?" -ne 0 ]; then
  stop_run "source sample extraction failed"
fi
awk -F '\t' 'NR>1 {print $1}' "$FROZEN/SAMPLE_POPULATION_COHORT.tsv" > "$WORK_DIR/EXPECTED_SAMPLES.joint_release.txt"
cp "$SITE_LISTS/donor_train.samples.txt" "$WORK_DIR/EXPECTED_SAMPLES.donor_train.txt"
cp "$SITE_LISTS/donor_validation.samples.txt" "$WORK_DIR/EXPECTED_SAMPLES.donor_validation.txt"
cp "$SITE_LISTS/hgdp_primary.samples.txt" "$WORK_DIR/EXPECTED_SAMPLES.hgdp_primary.txt"
if [ "$?" -ne 0 ]; then
  stop_run "expected sample-list staging failed"
fi

echo "=== MATERIALIZE SAMPLE PARTITIONS ==="
"$BCFTOOLS" view -S "$WORK_DIR/EXPECTED_SAMPLES.joint_release.txt" -Ob \
  -o "$WORK_DIR/JOINT_RELEASE_3264.bcf" "$WORK_DIR/ALL_SOURCE_EXACT_4091.bcf" 2>>"$LOG"
if [ "$?" -ne 0 ]; then stop_run "joint release materialization failed"; fi
"$BCFTOOLS" view -S "$WORK_DIR/EXPECTED_SAMPLES.donor_train.txt" -Ob \
  -o "$WORK_DIR/DONOR_TRAIN_2247.bcf" "$WORK_DIR/ALL_SOURCE_EXACT_4091.bcf" 2>>"$LOG"
if [ "$?" -ne 0 ]; then stop_run "donor train materialization failed"; fi
"$BCFTOOLS" view -S "$WORK_DIR/EXPECTED_SAMPLES.donor_validation.txt" -Ob \
  -o "$WORK_DIR/DONOR_VALIDATION_249.bcf" "$WORK_DIR/ALL_SOURCE_EXACT_4091.bcf" 2>>"$LOG"
if [ "$?" -ne 0 ]; then stop_run "donor validation materialization failed"; fi
"$BCFTOOLS" view -S "$WORK_DIR/EXPECTED_SAMPLES.hgdp_primary.txt" -Ob \
  -o "$WORK_DIR/HGDP_PRIMARY_768.bcf" "$WORK_DIR/ALL_SOURCE_EXACT_4091.bcf" 2>>"$LOG"
if [ "$?" -ne 0 ]; then stop_run "HGDP primary materialization failed"; fi

for bcf in "$WORK_DIR/JOINT_RELEASE_3264.bcf" "$WORK_DIR/DONOR_TRAIN_2247.bcf" \
  "$WORK_DIR/DONOR_VALIDATION_249.bcf" "$WORK_DIR/HGDP_PRIMARY_768.bcf"; do
  "$BCFTOOLS" index -f "$bcf" 2>>"$LOG"
  if [ "$?" -ne 0 ]; then stop_run "indexing failed: $bcf"; fi
done

echo "=== AUDIT COUNTS AND SAMPLE SETS ==="
printf 'role\tartifact\texpected_variants\tactual_variants\texpected_samples\tactual_samples\tsample_set_match\n' \
  > "$WORK_DIR/ARTIFACT_INVENTORY.v7.1.4.tsv"

audit_artifact() {
  role="$1"
  artifact="$2"
  expected_samples="$3"
  expected_list="$4"
  actual_list="$WORK_DIR/SAMPLES.${role}.txt"
  "$BCFTOOLS" query -l "$artifact" > "$actual_list" 2>>"$LOG"
  if [ "$?" -ne 0 ]; then stop_run "$role sample audit failed"; fi
  actual_samples="$(awk 'NF{n++} END{print n+0}' "$actual_list")"
  actual_variants="$("$BCFTOOLS" index -n "$artifact" 2>>"$LOG")"
  if [ "$?" -ne 0 ]; then stop_run "$role variant-count audit failed"; fi
  sort -u "$expected_list" > "$WORK_DIR/.expected.${role}.sorted"
  sort -u "$actual_list" > "$WORK_DIR/.actual.${role}.sorted"
  cmp -s "$WORK_DIR/.expected.${role}.sorted" "$WORK_DIR/.actual.${role}.sorted"
  if [ "$?" -eq 0 ]; then sample_match=PASS; else sample_match=FAIL; fi
  artifact_name="$(basename "$artifact")"
  printf '%s\t%s\t154850\t%s\t%s\t%s\t%s\n' \
    "$role" "$artifact_name" "$actual_variants" "$expected_samples" "$actual_samples" "$sample_match" \
    >> "$WORK_DIR/ARTIFACT_INVENTORY.v7.1.4.tsv"
  if [ "$actual_variants" -ne 154850 ] || [ "$actual_samples" -ne "$expected_samples" ] || [ "$sample_match" != PASS ]; then
    stop_run "$role inventory contract failed"
  fi
}

audit_artifact all_source "$WORK_DIR/ALL_SOURCE_EXACT_4091.bcf" 4091 "$WORK_DIR/EXPECTED_SAMPLES.all_source.txt"
audit_artifact joint_release "$WORK_DIR/JOINT_RELEASE_3264.bcf" 3264 "$WORK_DIR/EXPECTED_SAMPLES.joint_release.txt"
audit_artifact donor_train "$WORK_DIR/DONOR_TRAIN_2247.bcf" 2247 "$WORK_DIR/EXPECTED_SAMPLES.donor_train.txt"
audit_artifact donor_validation "$WORK_DIR/DONOR_VALIDATION_249.bcf" 249 "$WORK_DIR/EXPECTED_SAMPLES.donor_validation.txt"
audit_artifact hgdp_primary "$WORK_DIR/HGDP_PRIMARY_768.bcf" 768 "$WORK_DIR/EXPECTED_SAMPLES.hgdp_primary.txt"

"$BCFTOOLS" query -f '%CHROM:%POS:%REF:%ALT\n' "$WORK_DIR/ALL_SOURCE_EXACT_4091.bcf" \
  > "$WORK_DIR/MATERIALIZED_VARIANT_IDS.v7.1.4.txt" 2>>"$LOG"
if [ "$?" -ne 0 ]; then stop_run "materialized variant-key extraction failed"; fi
cmp -s "$FINAL_IDS" "$WORK_DIR/MATERIALIZED_VARIANT_IDS.v7.1.4.txt"
if [ "$?" -ne 0 ]; then stop_run "materialized variant key/order differs from frozen IDs"; fi

"$BCFTOOLS" view -h "$WORK_DIR/ALL_SOURCE_EXACT_4091.bcf" \
  > "$WORK_DIR/ALL_SOURCE_EXACT_4091.header.vcf" 2>>"$LOG"
if [ "$?" -ne 0 ]; then stop_run "output header audit failed"; fi
FORMAT_IDS="$(sed -n 's/^##FORMAT=<ID=\([^,>]*\).*/\1/p' "$WORK_DIR/ALL_SOURCE_EXACT_4091.header.vcf" | sort -u | paste -sd, -)"
INFO_IDS="$(sed -n 's/^##INFO=<ID=\([^,>]*\).*/\1/p' "$WORK_DIR/ALL_SOURCE_EXACT_4091.header.vcf" | sort -u | paste -sd, -)"
echo "format_ids=$FORMAT_IDS" | tee -a "$LOG"
echo "info_ids=$INFO_IDS" | tee -a "$LOG"
if [ "$FORMAT_IDS" != GT ]; then stop_run "output FORMAT set is not exactly GT"; fi
if [ -n "$INFO_IDS" ]; then stop_run "output retained INFO annotations"; fi

"$BCFTOOLS" query -i 'GT="mis"' -f '%CHROM\t%POS\n' "$WORK_DIR/JOINT_RELEASE_3264.bcf" \
  > "$WORK_DIR/MISSING_GT_SITES.joint_release.txt" 2>>"$LOG"
if [ "$?" -ne 0 ]; then stop_run "missing-GT audit query failed"; fi
MISSING_SITES="$(awk 'NF{n++} END{print n+0}' "$WORK_DIR/MISSING_GT_SITES.joint_release.txt")"

"$BCFTOOLS" view -P -Ou "$WORK_DIR/JOINT_RELEASE_3264.bcf" 2>>"$LOG" | \
  "$BCFTOOLS" query -f '%CHROM\t%POS\n' > "$WORK_DIR/UNPHASED_GT_SITES.joint_release.txt" 2>>"$LOG"
unphased_codes=("${PIPESTATUS[@]}")
echo "unphased_pipeline_codes=${unphased_codes[*]}" | tee -a "$LOG"
for code in "${unphased_codes[@]}"; do
  if [ "$code" -ne 0 ]; then stop_run "unphased-GT audit query failed"; fi
done
UNPHASED_SITES="$(awk 'NF{n++} END{print n+0}' "$WORK_DIR/UNPHASED_GT_SITES.joint_release.txt")"
echo "missing_gt_sites=$MISSING_SITES" | tee -a "$LOG"
echo "unphased_gt_sites=$UNPHASED_SITES" | tee -a "$LOG"

"$DATA_PY" "$R6/p0/build_materialization_audit_v7_1_4.py" \
  --inventory "$WORK_DIR/ARTIFACT_INVENTORY.v7.1.4.tsv" \
  --source-bcf "$SOURCE_BCF" \
  --frozen-ids "$FINAL_IDS" \
  --actual-ids "$WORK_DIR/MATERIALIZED_VARIANT_IDS.v7.1.4.txt" \
  --format-ids "$FORMAT_IDS" \
  --info-ids "$INFO_IDS" \
  --missing-sites "$MISSING_SITES" \
  --unphased-sites "$UNPHASED_SITES" \
  --output "$WORK_DIR/MATERIALIZATION_AUDIT.v7.1.4.json"
audit_exit=$?
if [ "$audit_exit" -ne 0 ]; then stop_run "materialization audit failed"; fi

find "$WORK_DIR" -maxdepth 1 -type f \
  ! -name 'MATERIALIZATION.v7.1.4.sha256' \
  -printf '%f\0' | sort -z | (
    cd "$WORK_DIR" || exit 1
    xargs -0 sha256sum
  ) > "$WORK_DIR/MATERIALIZATION.v7.1.4.sha256"
manifest_write_exit=$?
if [ "$manifest_write_exit" -ne 0 ]; then stop_run "materialization manifest write failed"; fi
(
  cd "$WORK_DIR" || exit 1
  sha256sum -c MATERIALIZATION.v7.1.4.sha256
) >/dev/null
manifest_verify_exit=$?
if [ "$manifest_verify_exit" -ne 0 ]; then stop_run "materialization manifest verification failed"; fi

mv "$WORK_DIR" "$OUT_FINAL"
publish_exit=$?
if [ "$publish_exit" -ne 0 ]; then stop_run "atomic publication failed"; fi
WORK_DIR=""

echo "=== V7.1.4 MATERIALIZATION PASS ==="
echo "output_dir=$OUT_FINAL"
cat "$OUT_FINAL/ARTIFACT_INVENTORY.v7.1.4.tsv"
cat "$OUT_FINAL/MATERIALIZATION_AUDIT.v7.1.4.json"
echo "GPU_USED=false"
