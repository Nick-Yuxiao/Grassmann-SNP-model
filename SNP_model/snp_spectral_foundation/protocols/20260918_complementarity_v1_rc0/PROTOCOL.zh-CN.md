# Concept Gate v1 协议（20260918_complementarity_v1_rc0）

## 0. 这一步在整条链路里的位置

Task Qualification（B 包）回答的是：**在协变量之上，raw dosage 还能不能带来样本外增量。**
它不能回答的是：这份增量里，哪一部分是"人群结构"，哪一部分是"功能注释"，两者互相之间是否可替代。

Concept Gate 就只回答这一件事，用最便宜的方式：

```
delta F | P  =  R2(PF) - R2(P)     有了人群结构以后，功能注释还有没有增量
delta P | F  =  R2(PF) - R2(F)     反过来
```

**这一步不训练任何 representation，不加载任何 checkpoint，不涉及 Grassmann。**
P 和 F 都是从已经冻结的 panel 直接算出来的确定性特征。

## 1. 五个 arm

| arm | 内容 | 来源 |
|-----|------|------|
| A | 协变量（sex、age、age²、array、assessment centre、PC1–40） | 与 B5 完全一致 |
| B | A + 冻结阈值下的 marginal-effect score | 阈值从 `TQ_RESULTS.json` 读入，**不重新挑** |
| P | A + 分块局部 PCA 成分 | `c1` |
| F | A + 功能注释 × dosage 聚合 | `c2` |
| PF | A + P + F | — |

## 2. 为什么 P 必须是"分块"PCA

arm A 里**已经含 PC1–40**。如果 P 用全基因组 PCA，它会和 A 高度共线，`delta P | F`
被构造性地压到 0——那不是一个实验结果，那是设计错误。

所以 P 定义为：在冻结 panel 上按**连续区块**（默认 100 个 variant 一块，不跨染色体）
做 PCA，每块取前 3 个成分，**只在 train 上拟合**，再投影到全体。这样 P 携带的是
局部单倍型结构，而全局祖先轴不表示这一层。

`c1` 会把 P 与 PC1–40 的 **canonical correlation** 写进 `P_SUMMARY.json`。这个数必须报告：
- 接近 1 → P 基本重复了 A 已有的东西，`delta P | F` 没有解释力；
- 明显小于 1 → P 确实带了新的结构。

特征宽度 = `ceil(panel_variants / block_size) * pcs_per_block`。2.5 万 SNP、block 100、
每块 3 个 → 750 列。内存不够就把 `--block-size` 调大。

## 3. F 的构造与硬约束

静态注释本身在人与人之间没有差异，所以注释只能当作**每个位点上的权重**作用在
这个人自己的 dosage 上：

```
F[i, a] = sum_j  standardised_dosage[i, j] * annotation[j, a]
```

默认按染色体分组聚合（`--group-by chromosome`），保留一点空间分辨率。

**硬约束：注释必须 phenotype-free。**
- 可以用：ChromHMM / 染色质状态、ATAC/DNase 可及性、跨物种保守性、基因结构。
- **不能用**：GWAS Catalog、同性状 eQTL、已有 PRS/PGS 权重、TWAS、fine-mapping 后验、
  LDSC baseline 中由关联导出的部分。

理由不是洁癖：由关联导出的注释把 outcome 塞进了 predictor，F 的增量就变成了循环论证。
`c2` 对任何没有声明 `phenotype_free=true` 的 track 直接 `SystemExit`；`c0` 对名字或路径
形似关联产物的 `--local` track 直接拒绝。

**join 只用位置，不用等位基因方向。** 本 cohort 的 A1/A2 来自一次没有显式 REF/ALT 模式的
BGEN 转换，allele-aware join 不可靠，positional join 不受影响。panel 是 **GRCh37**；
GRCh38 的注释会 join 到 0 个位点，`c2` 会直接停下来——这是预期的失败，不要绕过它。

## 4. 解码器

所有 arm 共用同一个解码器，差异只能来自特征本身：

- ridge，截距不惩罚；
- 惩罚网格 `1e-4 … 1e5`（10 档）；
- **λ 由 train 内部的 K 折交叉验证选出**（默认 5 折），每个 arm 各自选，
  评估 split 不参与拟合，也不参与选 λ；
- 每列都在 train 上标准化，让同一个 λ 对各列含义一致。

注意：arm A 在 B5 里是普通最小二乘，在这里是 ridge，所以**两边的 R2_A 不必相等**，
这已写进 `GATE_RESULTS.json` 的 `decoder.note`。

## 5. 判定规则

对每个 comparison：

```
通过  ⟺  delta_r2 > 0  且  paired 95% CI 下界 > 0
```

bootstrap 对同一批行同时重采样两个 arm（paired），2000 次，seed 记录在输出里。

verdict：

| 条件 | verdict |
|------|---------|
| `PF-P` 与 `PF-F` 都通过 | `BOTH_CARRIERS_ADD` |
| 只有 `PF-P` 通过 | `ONLY_FUNCTIONAL_ADDS` |
| 只有 `PF-F` 通过 | `ONLY_POPULATION_ADDS` |
| 都不通过 | `NEITHER_CARRIER_ADDS` |

## 6. 这个 verdict 是 developmental，不是 confirmatory

必须写清楚，并且已经写进 `GATE_RESULTS.json` 的 `status` 与 `why_developmental`：

validation 这条 split 到目前为止已经被用了两次——
1. Task Qualification 用它挑 B 的 p-value 阈值；
2. Concept Gate 用它做评估。

所以这里的 PASS 的含义是：**值得为它写一个 confirmatory 协议，并在封存 split 上跑一次**，
而不是"已经证实"。任何把这个结果当作已证实效应对外陈述的写法都是错的。

`c3` 里**没有**打开封存 split 的开关。`tq_test` / `bridge_holdout` / `test` 在
`SEALED_SPLITS` 里，`--eval-split` 指向它们会被直接拒绝，属于这些 split 的个体
在建设计矩阵之前就被剔除，并计入 `participants.dropped.sealed_split`。

## 7. trait 的角色

- **MCV = development / smoke-test trait。** 它在 validation 上 ΔR² 最大，但正因为如此，
  用它当唯一 confirmatory primary 是拿选择结果当证据。它的作用是把管线跑通、
  把量级看清楚。
- confirmatory primary trait 在 Task Qualification 的 `tq_test` 开一次之后，
  从 **QUALIFIED** 的 trait 里按预先写定的顺序选，写进 confirmatory 协议再动。

## 8. 输出边界（数据治理）

`c3` 只写四个文件，全部是汇总量：

```
GATE_RESULTS.json     完整记录：arm、comparison、participants 计数、sha256、seed
PENALTY_TRACE.json    每个 arm 的 λ 网格 CV 误差
arms.csv              trait, arm, r2, pearson, n_features, penalty, n_eval, eval_split
comparisons.csv       trait, comparison, delta_r2, ci_low, ci_high, replicates, seed
```

**没有任何个体级行**，包括假名化的。`--include-predictions` 这类开关在本包里不存在。
个体级 dosage、phenotype、预测值全部留在服务器。带出环境的只有上面这些汇总表。

如果之后确实出现"没有个体级预测就无法诊断"的问题，正确做法是在服务器内复算，
只导出 bootstrap 的分位数、seed、N 和 checksum。

## 9. 前置条件

1. Task B 的 panel 已冻结，`panel.int8.npy` / `panel_variants.tsv` / `panel_samples.tsv` 就位；
2. 每个 trait 有一份 `TQ_RESULTS.json`（validation-only 那次就够，它已含
   `selected_p_threshold`）；
3. 注释文件已下载到服务器（见 `ANNOTATION_SOURCES.zh-CN.md`）；
4. withdrawal 名单问题**不阻塞本步**——Concept Gate 不开封存 split。但它仍然阻塞
   `tq_test` 的那一次打开，这条冻结条款继续有效。

## 10. 不在本协议范围内

- 训练任何 encoder、加载任何 checkpoint；
- Grassmann 几何的任何部分；
- 用 phenotype 结果反选模型、层、秩、seed 或 panel；
- 打开 `tq_test` 或 `bridge_holdout`；
- 导出任何个体级数据。
