# R1 服务器执行步骤

_全部只读。任何一步都不训练、不读 phenotype、不打开 test。_

---

## 0. 上传与自检

```bash
cd /path/to/workdir
unzip r1_ukb_binding.zip      # 或 git clone 后 cd 到本目录
cd r1_ukb_binding
python3 --version             # 需要 >= 3.9
python3 -m unittest discover -s tests
```

预期：`OK`，13 个测试通过。失败就先停，把输出贴回来。

设置一个**不在 UKB 目录内**的输出根（正式产物可能几十 GB，放大盘）：

```bash
export OUT=/mnt/bigdisk/e0_runs        # 按你的服务器改
mkdir -p "$OUT"
df -h "$OUT"
```

## 1. 先做一次快速扫描（不算哈希）

先确认路径和检出结果对，再花时间算哈希：

```bash
python3 scripts/r1_inventory.py \
  --root /path/to/ukb/genotype \
  --root /path/to/ukb/kinship \
  --out-dir "$OUT/r1_probe" \
  --skip-hashes

cat "$OUT/r1_probe/INVENTORY_SUMMARY.md"
```

看三件事：

1. `Genotype filesets` 表里有没有出现你期望的 array panel，`primary eligible` 是不是 ✅
2. `build hint` 是不是唯一且符合预期（`GRCh37` 对应 UKB array/haplotype，`GRCh38` 需要单独确认）
3. `Kinship` 那一行有没有解析成功，component 数量是否合理

没找到文件就加 `--max-depth 8`，或者多给几个 `--root`。

## 2. 正式 inventory（带哈希，可后台）

大文件走 head/tail 指纹，小文件走全量 sha256。默认阈值 2 GiB，可调：

```bash
nohup python3 scripts/r1_inventory.py \
  --root /path/to/ukb/genotype \
  --root /path/to/ukb/kinship \
  --root /path/to/genetic_map \
  --out-dir "$OUT/r1_inventory_rc0" \
  --cohort-id UKB \
  --authorization-id <你的 application id> \
  --access-mode local_institutional_copy \
  --full-hash-max-bytes $((4*1024*1024*1024)) \
  > "$OUT/r1_inventory_rc0.stdout.log" 2> "$OUT/r1_inventory_rc0.stderr.log" &
echo $!
```

按协议要求：**不要前台空等**。过一会用 `tail -5 "$OUT/r1_inventory_rc0.stdout.log"` 看一次即可。

## 3. 独立验证 inventory

```bash
python3 scripts/r1_validate_inventory.py \
  --inventory "$OUT/r1_inventory_rc0/INVENTORY.json" \
  --output "$OUT/r1_inventory_rc0/R1_VALIDATION.json"
```

`status` 必须是 `PASS` 且 `r2_allowed` 为 `true` 才继续。常见 `FAIL` 与处理：

| 报错 | 处理 |
| --- | --- |
| `no parsable kinship file` | 用 `--root` 指到 kinship 目录；或该文件列名非标准，第 4 步显式给列号 |
| `kinship scan was truncated` | 加大 `--max-scan-lines` 重跑第 2 步 |
| `no directly genotyped ... fileset` | 确认你指的是 array calls 而不是 imputed；必要时用 `--array-variant-min/--array-variant-max` |
| `inventory output lives inside scanned root` | 换一个 `--out-dir` |
| `forbidden identifier-bearing key` | 不要手动编辑 `INVENTORY.json`，重跑 |

## 4. 冻结 family 轴

最简形式（kinship 文件是 KING `.kin0` 风格、带表头）：

```bash
python3 scripts/r2_build_family_manifest.py \
  --sample-source /path/to/panel.fam \
  --kinship /path/to/ukb_rel_aXXXXX.dat \
  --kinship-threshold 0.0442 \
  --seed 20260916 \
  --out-dir "$OUT/r2_family_rc0"
```

需要补充的常见选项：

- **撤回同意名单必须传**：`--exclude /path/to/w<app>_<date>.csv`，可重复
- 列名非标准时显式给 0-based 列号：`--kinship-id-columns 1 3 --kinship-value-column 7`
- 想在 ancestry/sex/array 内平衡：先准备一个 TSV（`sample_id` + 分层列），再加
  `--strata-file covars.tsv --strata-column ancestry --strata-column array`

`--sample-source` 支持 `.fam` / `.psam` / `.sample` / 每行一个 ID 的文本。

## 5. 独立验证 family 轴

```bash
python3 scripts/r2_validate_family_manifest.py \
  --family-manifest "$OUT/r2_family_rc0/family_manifest.tsv" \
  --kinship /path/to/ukb_rel_aXXXXX.dat \
  --exclude /path/to/w<app>_<date>.csv \
  --output "$OUT/r2_family_rc0/R2_FAMILY_VALIDATION.json"
```

必须同时满足：`status=PASS`、`component_cross_split_count=0`、`related_pairs_crossing_splits=0`、`partition_check.identical=true`。

## 6. 生成 CONTRACT 草稿

```bash
python3 scripts/r1_bind_contract.py \
  --inventory "$OUT/r1_inventory_rc0/INVENTORY.json" \
  --template ../CONTRACT.template.json \
  --family-summary "$OUT/r2_family_rc0/FAMILY_SUMMARY.json" \
  --output "$OUT/r1_inventory_rc0/CONTRACT.R1_DRAFT.json"
```

`genome_build`、`phasing_status`、`kinship_source/version`、`genetic_map_source` 默认留空，需要你确认后显式传参（例如 `--genome-build GRCh37 --phasing-status unphased_calls`）。草稿里的 `r1_unbound_fields` 就是进入 R2 之前还欠的字段清单，`run_authorized` 恒为 `false`。

## 7. 回传

把下面三个文件贴回来（不含参与者标识）：

```bash
cat "$OUT/r1_inventory_rc0/INVENTORY_SUMMARY.md"
cat "$OUT/r1_inventory_rc0/R1_VALIDATION.json"
cat "$OUT/r2_family_rc0/FAMILY_SUMMARY.json"
cat "$OUT/r2_family_rc0/R2_FAMILY_VALIDATION.json"
```

**不要**贴 `family_manifest.tsv`、`HASH_PLAN.tsv` 里的完整路径（如涉敏）或任何 `sample_id`。
