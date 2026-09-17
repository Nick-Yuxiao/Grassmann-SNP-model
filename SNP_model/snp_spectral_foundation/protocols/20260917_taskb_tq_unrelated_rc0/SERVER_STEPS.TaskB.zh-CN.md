# Task B 服务器执行步骤

_路径按 `c20-5090` 写。除 `b2` 抽列外，全部只读原始数据。_

---

## 0. 环境与变量

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

`b1`/`b2b`/`b3` 纯标准库，系统 `python3` 即可。`b2`/`b4`/`b5`/`b6` 需要 numpy 与 pyreadstat：

```bash
python3 -m ensurepip --user 2>/dev/null; python3 -m pip install --user numpy pyreadstat
# 失败则装 miniforge（后面 torch 也要）
# bash Miniforge3-Linux-x86_64.sh -b -p /home/tyuxiao/mf3
# /home/tyuxiao/mf3/bin/conda create -y -n taskb python=3.11 numpy pandas pyreadstat
python3 -m unittest discover -s tests     # 12 个全过
```

## 1. 字段扫描（已完成）

`22021` 在 `new` / `ukb2` / `ukb_alldata` 中。最小读取组合：trait 走 `ukb5.sas7bdat`（3.8 GiB），flag 走 `ukb2.sas7bdat`（52.7 GiB）。详见 `FIELD_BINDING_20260917.md`。

装上 pyreadstat 后可以补一次权威扫描（方法会从 `byte_scan_heuristic` 变成 `pyreadstat_metadata`）：

```bash
python3 scripts/b1_scan_sas_fields.py --sas "$SAS/ukb5.sas7bdat" --sas "$SAS/ukb2.sas7bdat" \
  --out-dir "$OUT/taskb_fields_authoritative"
```

## 2. 抽列

先 200 行冒烟确认标识列名（`n_eid` 还是 `eid`）：

```bash
python3 scripts/b2_extract_sas_columns.py --sas "$SAS/ukb5.sas7bdat" \
  --column n_eid --column n_30020_0_0 \
  --out "$OUT/taskb_extract/traits_smoke.tsv" --row-limit 200
head -3 "$OUT/taskb_extract/traits_smoke.tsv"
```

确认后跑全量，两个并行后台：

```bash
nohup python3 scripts/b2_extract_sas_columns.py --sas "$SAS/ukb5.sas7bdat" \
  --column n_eid --column n_30020_0_0 --column n_30010_0_0 --column n_30040_0_0 \
  --column n_30780_0_0 --column n_30760_0_0 --column n_30870_0_0 --column n_30690_0_0 \
  --out "$OUT/taskb_extract/traits.tsv" > "$OUT/taskb_extract/traits.log" 2>&1 &

nohup python3 scripts/b2_extract_sas_columns.py --sas "$SAS/ukb2.sas7bdat" \
  --column n_eid --column n_22021_0_0 --column n_22027_0_0 --column n_22019_0_0 --column n_22006_0_0 \
  --out "$OUT/taskb_extract/kinship_flag.tsv" > "$OUT/taskb_extract/kinship.log" 2>&1 &
```

## 3. 审计 22021 的编码（冻结过滤条件的前提）

**先审计，再过滤。** 不要跳过这一步直接跑 `b3`：

```bash
python3 scripts/b2b_audit_flag_coding.py \
  --table "$OUT/taskb_extract/kinship_flag.tsv" --id-column n_eid \
  --column n_22021_0_0 --column n_22027_0_0 --column n_22019_0_0 \
  --cross-reference "$COV" --cross-id-column eid \
  --out-dir "$OUT/taskb_flag_audit"

cat "$OUT/taskb_flag_audit/FLAG_AUDIT.md"
```

三个必须确认的点（详见 `FIELD_22021_CODING.zh-CN.md`）：

1. 观测取值只有 `{0, 1, 10, -1, <blank>}`，`undocumented_values_present` 为空
2. `keep_share` 落在 0.6–0.75 —— 偏离说明列选错或编码不同
3. `participant_overlap.coverage_of_cross_reference` 接近 1 —— 明显偏小说明 `ukb2` 不覆盖同一套人，应改用 `new.sas7bdat`

## 4. 冻结 split（含 bridge holdout）

```bash
python3 scripts/b3_build_unrelated_split.py \
  --covariates "$COV" --covariate-id-column eid --genotype-id-column iid \
  --fam "$FAM" \
  --kinship-flag "$OUT/taskb_extract/kinship_flag.tsv" \
  --kinship-flag-id-column n_eid --kinship-flag-column n_22021_0_0 \
  --qc-flag n_22027_0_0 --qc-flag n_22019_0_0 \
  --strata-column n_22001_0_0 --strata-column n_22006_0_0 \
  --reserve-bridge-holdout \
  --seed 20260917 --out-dir "$OUT/taskb_split_rc0"

cat "$OUT/taskb_split_rc0/SPLIT_SUMMARY.json"
```

有撤回同意名单必须加 `--exclude`。重点看 `relatedness_axis.dropped_by_value`（逐值account）和 `filtering.dropped_no_kinship_flag`（覆盖缺口）。

## 5. 建 SNP panel

先小规模冒烟：

```bash
python3 scripts/b4_build_snp_panel.py \
  --bed /data3/ukb_all/genotype/mnt_gsq/ukb_imp_chr1_v3.bed \
  --split-manifest "$OUT/taskb_split_rc0/tq_split_manifest.tsv" \
  --out-dir "$OUT/taskb_panel_smoke" --thin-bp 500000 --max-samples 20000
```

正式 panel，每条染色体一个 `--bed`（重复副本只用一份）：

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

panel 写 `/home/tyuxiao`（70 TB），不要写 `/data3`（剩 3 TB）。MAF 与缺失率只用 `train`。

## 6. Trait preflight（只读 train + validation）

```bash
PCS=$(for i in $(seq 1 40); do echo --covariate-column n_22009_0_${i}; done)

python3 scripts/b6_trait_preflight.py \
  --phenotype "$OUT/taskb_extract/traits.tsv" --phenotype-id-column n_eid \
  --trait-column n_30020_0_0 --trait-column n_30010_0_0 --trait-column n_30040_0_0 \
  --trait-column n_30780_0_0 --trait-column n_30760_0_0 --trait-column n_30870_0_0 \
  --trait-column n_30690_0_0 \
  --covariates "$COV" --covariate-id-column eid \
  $PCS --covariate-column n_21003_0_0 --square-column n_21003_0_0 \
  --categorical-column n_22001_0_0 --categorical-column n_54_0_0 \
  --split-manifest "$OUT/taskb_split_rc0/tq_split_manifest.tsv" \
  --out-dir "$OUT/taskb_preflight"

cat "$OUT/taskb_preflight/PREFLIGHT.md"
```

确认 `splits_left_closed` 含 `tq_test` 与 `bridge_holdout`。据此挑 1–3 个 trait，写进 `CONFIG.json`。

## 7. 资格检验：先只看 validation

```bash
python3 scripts/b5_task_qualification.py \
  --panel-dir "$OUT/taskb_panel_rc0" \
  --covariates "$COV" --covariate-id-column eid \
  $PCS --covariate-column n_21003_0_0 --square-column n_21003_0_0 \
  --categorical-column n_22001_0_0 --categorical-column n_54_0_0 \
  --phenotype "$OUT/taskb_extract/traits.tsv" --phenotype-id-column n_eid \
  --trait-column n_30020_0_0 \
  --test-split tq_test \
  --out-dir "$OUT/taskb_tq_hgb"
```

看 `validation_scan` 各阈值的 `delta_r2` 走势是否合理。

## 8. 开 test（每个 trait 只能一次）

`CONFIG.template.json` 填好存成 `CONFIG.json`（trait 列表、split 与 panel 的 SHA-256、eligibility 规则都写死）之后：

```bash
python3 scripts/b5_task_qualification.py \
  ... 同上 ... --out-dir "$OUT/taskb_tq_hgb_test" --allow-test
```

判定：`eligible = ΔR² > 0 且 paired 95% CI 下界 > 0`。结果里应有 `splits.never_touched: ["bridge_holdout"]`。

## 9. 回传

```bash
cat "$OUT/taskb_flag_audit/FLAG_AUDIT.md"
cat "$OUT/taskb_split_rc0/SPLIT_SUMMARY.json"
cat "$OUT/taskb_panel_rc0/PANEL_SUMMARY.json"
cat "$OUT/taskb_preflight/PREFLIGHT.md"
cat "$OUT/taskb_tq_hgb_test/TQ_RESULTS.json"
```

都不含参与者标识。**不要贴** `tq_split_manifest.tsv`、`traits.tsv`、`kinship_flag.tsv`。
