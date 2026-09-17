# Task B 服务器执行步骤

_路径按 `c20-5090` 写。除 `b2` 抽列外，全部只读原始数据。_

---

## 0. 上传与自检

```bash
cd /home/tyuxiao/Grassmann_model/SNP_model/incoming
tar -xzf taskb_tq.tar.gz && cd 20260917_taskb_tq_unrelated_rc0
sha256sum -c MANIFEST.sha256

export OUT=/home/tyuxiao/Grassmann_model/e0_runs
export SAS=/data3/ukb_all/phenotype/mnt_gsq_xinjianwenjianjia
export COV=/data3/ukb_all/phenotype/derived/ukb_covariates_chr1fam_aligned.tsv
export FAM=/data3/ukb_all/genotype/mnt_gsq/ukb_imp_chr1_v3.fam
mkdir -p "$OUT"
```

`b1`/`b3` 不需要 numpy，可以立刻跑。`b4`/`b5` 需要 numpy，见第 4 步前的环境准备。

## 1. 找 22021 在哪个 SAS 文件里

```bash
python3 scripts/b1_scan_sas_fields.py \
  --sas-dir "$SAS" --out-dir "$OUT/taskb_fields" --scan-bytes $((300*1024*1024))

cat "$OUT/taskb_fields/FIELD_SCAN.md"
```

小文件（`data.sas7bdat`、`ukb4/ukb5`）先跑完，大文件慢一些。**输出里 `22021 是否存在` 那一行决定下一步**：

- 存在 → 继续第 2 步
- 不存在 → 停，走 `KING_FALLBACK.zh-CN.md`，或向数据管理方要 `ukb44430_rel_*.dat`

## 2. 抽列（需要 pyreadstat 或 pandas）

先用 `--row-limit` 做 200 行冒烟，确认列名对：

```bash
python3 scripts/b2_extract_sas_columns.py \
  --sas "$SAS/<b1 指出的文件>.sas7bdat" \
  --column n_eid --column n_22021_0_0 --column n_22027_0_0 --column n_22019_0_0 \
  --out "$OUT/taskb_extract/kinship_flag_smoke.tsv" --row-limit 200
head -3 "$OUT/taskb_extract/kinship_flag_smoke.tsv"
```

确认无误后跑全量（后台，大文件慢）：

```bash
nohup python3 scripts/b2_extract_sas_columns.py \
  --sas "$SAS/<文件>.sas7bdat" \
  --column n_eid --column n_22021_0_0 --column n_22027_0_0 --column n_22019_0_0 \
  --out "$OUT/taskb_extract/kinship_flag.tsv" \
  > "$OUT/taskb_extract/kinship.log" 2>&1 &
```

Trait 同理，单独抽一份（**这份文件只属于 Task B**）：

```bash
nohup python3 scripts/b2_extract_sas_columns.py \
  --sas "$SAS/<文件>.sas7bdat" \
  --column n_eid --column n_30020_0_0 --column n_30010_0_0 --column n_30040_0_0 \
  --out "$OUT/taskb_extract/traits.tsv" \
  > "$OUT/taskb_extract/traits.log" 2>&1 &
```

若抽出的 id 列名是 `n_eid` 而 covariates 里是 `eid`，在后面用 `--kinship-flag-id-column n_eid` / `--phenotype-id-column n_eid` 指定即可。

## 3. 冻结无亲缘子集的 split（纯标准库）

```bash
python3 scripts/b3_build_unrelated_split.py \
  --covariates "$COV" \
  --covariate-id-column eid --genotype-id-column iid \
  --fam "$FAM" \
  --kinship-flag "$OUT/taskb_extract/kinship_flag.tsv" \
  --kinship-flag-id-column n_eid --kinship-flag-column n_22021_0_0 \
  --qc-flag n_22027_0_0 --qc-flag n_22019_0_0 \
  --strata-column n_22001_0_0 --strata-column n_22006_0_0 \
  --seed 20260917 \
  --out-dir "$OUT/taskb_split_rc0"

cat "$OUT/taskb_split_rc0/SPLIT_SUMMARY.json"
```

**有撤回同意名单一定要加 `--exclude`。** 预期保留约 34 万人，`dropped_related_or_unassessed` 约 14 万。

## 4. 准备 numpy 环境

```bash
cd /home/tyuxiao
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
bash Miniforge3-Linux-x86_64.sh -b -p /home/tyuxiao/mf3
/home/tyuxiao/mf3/bin/conda create -y -n taskb python=3.11 numpy pandas pyreadstat
source /home/tyuxiao/mf3/etc/profile.d/conda.sh && conda activate taskb
```

## 5. 建 SNP panel

先小规模确认能跑通（单染色体、少量样本）：

```bash
cd /home/tyuxiao/Grassmann_model/SNP_model/incoming/20260917_taskb_tq_unrelated_rc0
python3 scripts/b4_build_snp_panel.py \
  --bed /data3/ukb_all/genotype/mnt_gsq/ukb_imp_chr1_v3.bed \
  --split-manifest "$OUT/taskb_split_rc0/tq_split_manifest.tsv" \
  --out-dir "$OUT/taskb_panel_smoke" \
  --thin-bp 500000 --max-samples 20000
```

通过后跑正式 panel。每条染色体给一个 `--bed`，重复副本只用一份：

```bash
nohup python3 scripts/b4_build_snp_panel.py \
  --bed /data3/ukb_all/genotype/mnt_gsq/ukb_imp_chr1_v3.bed \
  --bed /data3/ukb_all/genotype/data2_external_drive_copies/gsq_root/ukb_imp_chr2_v3.bed \
  --bed /data3/ukb_all/genotype/data1_ukb_crosschr_fixed/ukb_imp_chr3_v3.bed \
  --bed /data3/ukb_all/genotype/data2_external_drive_copies/onetouch_root/ukb_imp_chr4_v3.bed \
  --bed /data3/ukb_all/genotype/data1_ukb_crosschr_fixed/ukb_imp_chr5_v3.bed \
  --bed /data3/ukb_all/genotype/mnt_onetouch/ukb_imp_chr6_v3.bed \
  --bed /data3/ukb_all/genotype/mnt_onetouch/ukb_imp_chr7_v3.bed \
  --bed /data3/ukb_all/genotype/external_new_chr8_9_10/gsq/ukb_imp_chr8_v3.bed \
  --bed /data3/ukb_all/genotype/external_new_chr8_9_10/gsq/ukb_imp_chr9_v3.bed \
  --bed /data3/ukb_all/genotype/external_new_chr8_9_10/onetouch/ukb_imp_chr10_v3.bed \
  $(for c in 11 12 13 14 15 16 17 18 19 20 21 22; do \
      echo --bed /data3/ukb_all/genotype/新建文件夹/ukb_imp_chr${c}_v3.bed; done) \
  --split-manifest "$OUT/taskb_split_rc0/tq_split_manifest.tsv" \
  --out-dir "$OUT/taskb_panel_rc0" \
  --thin-bp 50000 --max-samples 120000 \
  > "$OUT/taskb_panel_rc0.log" 2>&1 &
```

`--thin-bp 50000` 全基因组约取 5 万个位点；`--max-samples 120000` 下 panel 约 6 GB。磁盘写 `/home/tyuxiao`（70 TB），不要写 `/data3`（只剩 3 TB）。

## 6. 资格检验：先只看 validation

```bash
PCS=$(for i in $(seq 1 40); do echo --covariate-column n_22009_0_${i}; done)

python3 scripts/b5_task_qualification.py \
  --panel-dir "$OUT/taskb_panel_rc0" \
  --covariates "$COV" --covariate-id-column eid \
  $PCS \
  --covariate-column n_21003_0_0 --square-column n_21003_0_0 \
  --categorical-column n_22001_0_0 --categorical-column n_54_0_0 \
  --phenotype "$OUT/taskb_extract/traits.tsv" --phenotype-id-column n_eid \
  --trait-column n_30020_0_0 \
  --out-dir "$OUT/taskb_tq_hgb"
```

看 `validation_scan` 各阈值的 `delta_r2`。数值合理、负对照方向正确之后，再填 `CONFIG.json`。

## 7. 开 test（每个 trait 只能一次）

先把 `CONFIG.template.json` 填好存成 `CONFIG.json`（trait 列表、SESOI、split 与 panel 的 SHA-256 都要写死），然后：

```bash
python3 scripts/b5_task_qualification.py \
  ... 同上 ... \
  --out-dir "$OUT/taskb_tq_hgb_test" \
  --allow-test
```

输出 `TQ_RESULTS.json` 含 `r2_A`、`r2_B`、`delta_r2`、95% CI、负对照与判定。同目录再开会被拒绝。

## 8. 回传

```bash
cat "$OUT/taskb_fields/FIELD_SCAN.md"
cat "$OUT/taskb_split_rc0/SPLIT_SUMMARY.json"
cat "$OUT/taskb_panel_rc0/PANEL_SUMMARY.json"
cat "$OUT/taskb_tq_hgb_test/TQ_RESULTS.json"
```

都不含参与者标识。**不要贴** `tq_split_manifest.tsv`、`traits.tsv`、`kinship_flag.tsv`。
