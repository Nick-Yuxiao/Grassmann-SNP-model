# 备选路径：用现有 hard-call 自建 KING kinship

_仅在 UKB field `22021` 不可得、且官方 `ukb_rel_*.dat` 也拿不到时使用。_

---

## ⚠️ 先说折扣

在 **imputed hard-call** 上估计 kinship，精度低于在 array 直接分型上估计：

- imputation 会把罕见基因型平滑掉，`IBS0` 被压低，**亲缘系数系统性偏高**
- 受 reference panel 影响，不同 ancestry 的偏差不一致
- `0.0442` 这个阈值是为 array KING 校准的，用在这里应视为近似

因此自建 kinship 的定位是**保守过滤器**：宁可多划掉一些人，也不要让真实亲缘对跨 split。建议把阈值放宽到 `0.0353`（更激进地剔除），并在协议里登记该偏差。拿到官方文件后应重新冻结 split。

## 💰 成本

`--make-king-table` 的代价随样本数平方增长。487,409 人约 1.19×10¹¹ 对。即使配合 `--king-table-filter` 只输出超阈值的对，全量计算仍是**数十 CPU-小时量级**，需要大内存和充足临时空间。开始前先确认机器规格和排队策略。

先决条件：plink2 是单个静态二进制，不需要管理员权限。

```bash
mkdir -p /home/tyuxiao/bin && cd /home/tyuxiao/bin
curl -L -o plink2.zip "https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64_20241114.zip"
unzip -o plink2.zip && ./plink2 --version
export PATH=/home/tyuxiao/bin:$PATH
```

## 1. 逐染色体 LD pruning

```bash
GENO=/data3/ukb_all/genotype
OUT=/home/tyuxiao/Grassmann_model/e0_runs/king
mkdir -p "$OUT"

for f in $(find "$GENO" -name "ukb_imp_chr*_v3.bed" | sort -u); do
  stem="${f%.bed}"; chr=$(basename "$stem" | sed 's/.*chr\([0-9]*\)_.*/\1/')
  plink2 --bfile "$stem" \
    --snps-only just-acgt --maf 0.05 --geno 0.02 --hwe 1e-6 \
    --exclude-long-range-ld \
    --indep-pairwise 1000kb 0.1 \
    --threads 16 --memory 64000 \
    --out "$OUT/prune_chr$chr"
done
```

同一染色体有重复副本时，`sort -u` 之外还要人工确认只用一份。

## 2. 提取并合并 pruned panel

```bash
for f in ...; do
  plink2 --bfile "$stem" --extract "$OUT/prune_chr$chr.prune.in" \
    --make-bed --threads 16 --out "$OUT/pruned_chr$chr"
done

ls "$OUT"/pruned_chr*.bed | sed 's/\.bed$//' > "$OUT/merge_list.txt"
plink2 --pmerge-list "$OUT/merge_list.txt" bfile \
  --make-bed --threads 16 --out "$OUT/pruned_all"
```

目标是保留约 10 万个近独立位点；`pruned_all` 大约 12 GB。

## 3. 计算 KING 并只导出超阈值的对

```bash
nohup plink2 --bfile "$OUT/pruned_all" \
  --make-king-table --king-table-filter 0.0353 \
  --threads 32 --memory 200000 \
  --out "$OUT/king" \
  > "$OUT/king.stdout.log" 2> "$OUT/king.stderr.log" &
```

输出 `king.kin0` 的列为 `#FID1 IID1 FID2 IID2 NSNP HETHET IBS0 KINSHIP`。

## 4. 接回既有 family 轴代码

`.kin0` 正是 E0 `R2` builder 的输入格式，可以直接冻结 family manifest：

```bash
python3 r1_ukb_binding/scripts/r2_build_family_manifest.py \
  --sample-source <panel>.fam \
  --kinship "$OUT/king.kin0" \
  --kinship-threshold 0.0442 \
  --strata-file <covariates>.tsv --strata-column n_22001_0_0 \
  --seed 20260917 \
  --out-dir "$OUT/family_rc0"

python3 r1_ukb_binding/scripts/r2_validate_family_manifest.py \
  --family-manifest "$OUT/family_rc0/family_manifest.tsv" \
  --kinship "$OUT/king.kin0" \
  --output "$OUT/family_rc0/VALIDATION.json"
```

列名非标准时用 `--kinship-id-columns 1 3 --kinship-value-column 7` 显式指定。

## 5. 必须登记的内容

`CONFIG.json` 中记录：kinship 来源为自建、plink2 版本、pruning 参数、pruned 位点数、`--king-table-filter` 阈值、`king.kin0` 的 SHA-256，以及「在 imputed hard-call 上估计、偏差方向为高估」这一条。拿到官方文件后重新冻结并对比两套 components。
