# 候选架构的验证合同

> ⚠️ **状态更新（2026-09-15）：** 本文件只保留为未来新架构的 Phase II 工程合同，不再是当前研究主线。当前唯一主线由 `protocols/20260915_0526_stage2_bridge_rc1/PROTOCOL.zh-CN.md` 冻结；在该协议的 phenotype bridge 通过前，不执行这里的 foundation-model attribution 或新架构训练。

本目录实现的是**待否证候选**，不是已被现有实验支持的最终模型。任何全规模训练前，按顺序执行以下 gate；前一层失败时不解释后一层。

## Gate 0：数据与切分

- ALT allele frequency、MAF、LD block、标准化量和变异筛选仅在训练个体上估计；AF 用于 dosage residual，MAF 用于频率分层，二者不可混用。
- 家系/亲缘个体不得跨 train、validation、test；跨 ancestry 与跨 cohort 另报。
- phase-aware 任务使用 phased haplotype 输入；dosage 输入不能证明模型利用了 phase。
- 所有模型共享 SNP panel、切分、mask、优化预算和调参预算。

## Gate 1：local genotype pretraining

冻结 local encoder，用新 mask 在 held-out people 上比较：

1. masked-genotype loss / accuracy；
2. train-only locus-mode 或 HWE marginal baseline；
3. 保持 MAF、破坏 LD 后的性能差；
4. block 内 context shuffle 后的性能差。

必须先证明相对 marginal baseline 有预先规定的实际增量，才进入 spectral memory。

## Gate 2：结构归因，而不是结构命名

对同一冻结 encoder 做以下干预，并测预训练、embedding 与下游性能保留率：

| 干预 | 主要保留 | 主要破坏 |
|---|---|---|
| 每 SNP 跨个体独立置换 | MAF / genotype count | LD、haplotype |
| 整 block 跨个体置换 | block 内 local LD | block 间依赖 |
| dosage 不变的 phase scrambling | dosage | phase / cis-trans |
| ancestry 内置换与 MAF matching | 指定边际 | population contrast 的其他部分 |
| 10 kb–1 Mb context destruction | 指定半径内结构 | 更长程结构 |

干预应做成成对、固定随机种子的 factorial ablation；报告效应量和区间，不能只看显著性。

## Gate 3：memory 公平比较

在相同输入、相同 memory 维度/字节、相同 global mixer、相同 decoder 下比较：

- `spectral_only`：Grassmann projector + eigenvalues + MAF；
- `residual_only`：有符号、位置对齐的 dosage residual sidecar；
- `hybrid`：上述两者融合；
- `aligned_control`：保留位置与方向的强对照。

主要问题不是 hybrid 是否胜过一个弱 baseline，而是它是否在新染色体/新 cohort 上胜过 `aligned_control`，且增益能归因到 spectral branch。

## Gate 4：global mixer 与 phenotype conditioning

- 用 block shuffle、距离分层 context destruction 比较 local-only 与 global slots。
- conditioning 只能使用 trait identity 或部署时可得的 context，绝不能把真实 outcome 输入模型。
- phenotype-conditioned effect 表示若用于推断/解释，必须 cross-fit；产生某人的表示时，projector 未见该人的 outcome。
- prediction、association/interaction characterization 分开评价。

## 最小消融矩阵

这是一个 `2 × 2 × 2` 析因设计：spectral branch（开/关）× residual sidecar（开/关）× global mixer（开/关）。另加 `aligned_control` 作为主动对照。至少使用多个数据种子与训练种子，并以 person/family 为独立重复层级。不要把 SNP 数或重复 mask 当独立样本。

停止条件：local pretraining 不胜 marginal baseline；hybrid 不胜 aligned control；global branch 在长程破坏实验中无增量；或 cross-cohort 增量不复现。任何一个成立，都应收缩模型而不是继续扩参。
