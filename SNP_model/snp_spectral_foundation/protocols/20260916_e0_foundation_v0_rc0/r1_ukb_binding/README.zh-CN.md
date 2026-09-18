# R1/R2-family UKB 绑定执行包

_协议 `20260916_e0_foundation_v0_rc0` · 里程碑 `R1` 与 `R2` 的 family 轴 · 状态 `READ_ONLY_EVIDENCE_ONLY`_

---

## 📋 这个包做什么

`STATUS_20260916.md` 里唯一的阻断是「UKB genotype、kinship、genetic map 与授权路径尚未绑定」。本包就是把这条阻断变成证据：

1. **R1 只读 inventory**：扫描你给的 UKB 目录，报出 fileset、样本数、variant 数、染色体覆盖、build 提示、kinship 结构、磁盘路径、大小与哈希计划，并顺带绑定 CUDA / Java / plink2 / bcftools 等运行时事实。
2. **R1 独立验证**：不信任生成器，重新检查结构完整性、标识符卫生、哈希覆盖率与只读边界，给出 `PASS/FAIL` 与 `r2_allowed`。
3. **R2 family 轴**：用授权的 kinship 文件做 connected components，按 `70/15/15` 在 strata 内冻结 `family_manifest.tsv`。
4. **R2 family 独立验证**：从 kinship 文件**重新**推导 components 并与 manifest 比对，检查跨 split 亲缘对、component_size、排除名单与比例。
5. **CONTRACT 草稿**：把机器能确定的字段填进 `CONTRACT.json` 草稿，人必须决定的字段保持 `null`。

## 🚫 硬边界

这个包在设计上做不到以下事情，验证器也会强制：

- **不写入任何被扫描的目录**；`--out-dir` 落在扫描根内会直接拒绝运行
- **不读 phenotype**，不打开 test 数据，不训练任何模型
- **不输出参与者标识**：所有 JSON/MD 报告只含计数、路径、大小与哈希；验证器会扫描禁用字段名与可疑的长字符串列表并判 `FAIL`
- **不签发 `RUN_AUTHORIZED`**：`run_authorized` 在本包所有输出里恒为 `false`

`family_manifest.tsv` 本身当然含 `sample_id`，它必须留在服务器上，**不要**提交到仓库或贴进聊天。可以回传的是 `FAMILY_SUMMARY.json` 与验证报告。

## ⚠️ 一个必须先确认的科学点：用 array calls，不要用 imputed

E0 的正向任务是 masked genotype prediction，`M1b` 是 Beagle haplotype-HMM。如果把 **imputed dosage** 当作 masking 目标，那么 target 本身就是一个 HMM 的输出，E0 与它的最强基线会变成循环比较，`δ_X` 失去意义。

所以 inventory 会给每个 fileset 打标签：

| role | 判据 | 能否做 E0 primary |
| --- | --- | --- |
| `directly_genotyped_candidate` | variant 数落在 array 量级（默认 100k–3M） | ✅ primary panel |
| `imputed_candidate` | 文件名含 `imp/mfi/dosage`，或 variant 数 > 5M | ❌ 不可做 masking target |
| `phased_haplotype_candidate` | 文件名含 `hap/phase/shapeit` | ❌ primary；但它是 `M1b` reference 的候选 |

阈值可用 `--array-variant-min/--array-variant-max` 调整；分类结果不会静默改变任何 contract 字段。

同理，`data.phasing_status` 只能由你显式写入。协议规定：unphased dosage 只能支持「局部 LD / 联合基因型上下文」结论，phased claim 需要输入和 truth 都是合格 phased haplotypes。

## ▶️ 三步执行

具体命令见 [`SERVER_STEPS.R1.zh-CN.md`](SERVER_STEPS.R1.zh-CN.md)。最短形式：

```bash
python3 scripts/r1_inventory.py \
  --root /path/to/ukb/genotype --root /path/to/ukb/kinship \
  --out-dir "$OUT/r1_inventory_rc0" \
  --cohort-id UKB --authorization-id <你的 application id>

python3 scripts/r1_validate_inventory.py \
  --inventory "$OUT/r1_inventory_rc0/INVENTORY.json" \
  --output "$OUT/r1_inventory_rc0/R1_VALIDATION.json"

python3 scripts/r2_build_family_manifest.py \
  --sample-source <panel>.fam --kinship <kinship file> \
  --seed 20260916 --out-dir "$OUT/r2_family_rc0"
```

只依赖 Python 3.9+ 标准库，没有第三方依赖；torch 只在存在时用于探测 CUDA，缺失不影响 R1。

包内的 `CONTRACT.template.json` 是协议目录同名文件的副本，便于本包在服务器上自包含运行；以协议目录的版本为准。

## 📤 跑完回传什么

把这三个文件贴回来（都不含参与者标识），我据此写 `R1` 状态报告并出 R2 的 block 轴代码：

- `INVENTORY_SUMMARY.md`
- `R1_VALIDATION.json`
- `FAMILY_SUMMARY.json` 与 `R2_FAMILY_VALIDATION.json`

如果路径里有你不想公开的目录名，可以在贴出前替换掉路径字段，其它计数保持原样即可。

## ⏭️ 下一步（R2 block 轴，本包不含）

block 轴必须等 R1 报出真实的染色体覆盖、variant 密度和 genetic map 版本，才能：

1. 只用 train families 或冻结 reference 构造 LD blocks
2. 按 chromosome / 位置 / block size / MAF / LD density 分层抽取**连续 held-out runs**
3. 以 `max(encoder receptive field, 相邻 LD block, 1 cM)` 在 run 外边界设 guard
4. 用 train families 审计跨边界 `r² ≥ 0.1`，不通过就扩 guard 或并成 atomic block

生成的 `block_manifest.tsv` 会与本包的 `family_manifest.tsv` 一起送进既有的 [`scripts/validate_manifests.py`](../scripts/validate_manifests.py)；本包的 family manifest 列名就是按那个验证器的契约写的。

## 🧪 自测

```bash
python3 -m unittest discover -s tests -v
```

13 个合成测试覆盖：fileset 与 build 提示识别、kinship component 计数、只读边界拒绝、不写入扫描根、标识符泄漏必须 `FAIL`、imputed 量级不得进入 E0 primary、contract 草稿保持未绑定字段为 `null`、family split 确定性、排除名单生效、以及故意制造的跨 split 亲缘泄漏必须 `FAIL`。
