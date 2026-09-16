# 0526 Stage 2 phenotype bridge：开发 Pilot rc2

_冻结日期：2026-09-15｜状态：`FROZEN_PILOT_STARTED`｜取代过度防御的 rc1，但不删除 rc1_

## 唯一主问题

在固定的真实个体级 genotype–molecular-phenotype 数据上，genotype-only 预训练得到的隐藏表示 `H`，是否相对于同一 SNP panel 的简单加性 dosage ridge 增加 held-out predictive R²？

主 estimand 只有：

\[
\Delta R^2_{H}=R^2(C)-R^2(B).
\]

`A` 是非遗传参照；`D/E` 是机制性次要比较。它们都不是 `C−B` 的先决 Gate。

## 数据与声明边界

- 数据：Bioconductor `GeuvadisTranscriptExpr 1.40.0`，GEUVADIS CEU、chr19、个体级 genotype 与 transcript expression。
- 角色：真实分子表型 **development Pilot**。
- 不能声称：外部复制、临床表型预测、全基因组有效、Grassmann 优势、因果或独立 biological effect。
- GEUVADIS 与 1000 Genomes 共享参与者；这不妨碍开发 Pilot，但禁止称为相对 1000G 的独立验证。

## 冻结执行顺序

1. 校验源包哈希、版本与内容；只保留 chr19 双等位 genotype/表达共同个体。
2. 以 seed `5262026` 建立一次固定的 `70/15/15` individual-disjoint split。
3. 只用 train 个体拟合缺失填充、AF、标准化、PCA 与 trait 选择。
4. trait 选择不查看任何 `H` 或 test 结果：在可用 transcript 中选择 train 方差最大的有限表达 transcript；并报告其 ID。
5. 只用 train genotype 做 masked-genotype pretraining；训练结束后保存 checkpoint 与 SHA-256，随后冻结 encoder。
6. 对 train/validation/test 各提取一次 `H`；不按 phenotype 选择 layer、rank 或 seed。
7. 所有臂使用相同 split、ridge family、alpha grid、validation MSE 和 test 打开时点。
8. test 只评价一次；所有预先指定结果都报告，不以失败触发改标准。

## 比较臂

| Arm | 输入 | 角色 |
| --- | --- | --- |
| A | train-fitted genotype PCs（若数据包有 sex 则再加 sex） | 非遗传/结构参照 |
| B | A + 同一 panel 的 additive ALT dosage | 最强简单主基线 |
| C | A + 冻结 encoder 的完整、按 SNP 顺序对齐的 `H` | 主受试表示 |
| D | A + train-fitted PCA 压缩的 aligned `H` | 次要压缩参照 |
| E | A + 每个 block 的 `H` 子空间投影矩阵上三角元素与谱值 | 次要 Grassmann 表示 |

`D/E` 使用同一冻结 rank；rank 由实现可容纳的 `min(4, hidden_width, block_size-1)` 确定，不根据 phenotype 表现搜索。

## Decoder、指标与判定

- Decoder：Ridge；alpha grid 为 `10^-6 ... 10^6`；validation MSE 最小，平局取更大 alpha。
- outcome：train mean/SD 标准化，validation/test 原样应用。
- test predictive R²：`1 - SSE / sum((y_test - mean(y_train))^2)`。
- 不确定性：固定预测上的 5,000 次 paired individual bootstrap；样本量小导致区间宽是 Pilot 结果，不是新增 Gate 的理由。
- **支持 bridge：** `C−B > 0` 且 paired 95% CI 下界 `> 0`。
- **不支持 bridge：** 上述条件不成立。不得事后加入/删除样本、改 trait、改 seed、改 alpha grid 或另找阈值挽救本次结果。
- 不设 `0.005` SESOI，不做 Holm，不要求 `B−A` 先通过。

## 失败只停止什么

- `C−B` 未通过：只说明本数据、panel、encoder 与样本量下未建立 `H` 的 phenotype 增量；不否定 genotype pretraining 或未来独立数据上的 phenotype 工作。
- `E` 不如 `D`：只停止本设置的 Grassmann 默认压缩主张；不否定 C/D 或 foundation route。
- 运行/数据错误：结果标为 invalid 并保留日志；修复只能产生新版本，不能覆盖本次输出。

## 后续授权

只有 `C−B` 通过后，才在新版本和新 test 上考虑 biology prior 或 CropARNet-style phenotype weighting。当前 Pilot 禁止动态 phenotype weighting。
