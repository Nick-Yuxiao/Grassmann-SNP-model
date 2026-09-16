#!/usr/bin/env bash
set -euo pipefail

source /home/yuxiao_tan/miniforge3/etc/profile.d/conda.sh
conda activate epinformer_repro

project_root='/mnt/c/Users/Yuxiao Tan/OneDrive/桌面/Qinghua_bioinfo/Grassmann_model'
protocol_rel='snp_spectral_foundation/protocols/20260916_e0_foundation_v0_rc0'
output_root='/mnt/f/Yuxiao Tan/Document/Grassmann_model_external/e0_runs/m0_validation_2blocks_rc0'

cd "$project_root"
python "$protocol_rel/scripts/run_m0_development.py" \
  --vcf '/mnt/f/Yuxiao Tan/Document/Grassmann_model_external/1000g_phase3_chr22/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.parallel.vcf.gz' \
  --family-manifest "$protocol_rel/development/1000g_chr22_rc3_cm/family_manifest.tsv" \
  --block-manifest "$protocol_rel/development/1000g_chr22_rc3_cm/block_manifest.tsv" \
  --output-dir "$output_root" \
  --evaluation-family-split validation \
  --evaluation-block-split validation \
  --mask-seeds 20260916 \
  --max-blocks 2
