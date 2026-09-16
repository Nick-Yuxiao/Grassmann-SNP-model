# GEUVADIS chr19 development Pilot rc1 结果

## 冻结主结论

本次结果**不支持 0526 Stage 2 bridge**。冻结主比较 `C-B` 的
ΔR² = **-0.006803**，paired 95% bootstrap CI =
[-0.059185, 0.044373]。判定严格使用运行前规则：
点估计大于 0 且区间下界大于 0；没有追加 0.005 门槛、P1 前置门或多重校正门。

## 数据与模型

- 数据：Bioconductor `GeuvadisTranscriptExpr 1.40.0`，CEU chr19，真实个体级 genotype–transcript expression。
- 样本：train/validation/test = 63/13/15，家系隔离。
- 性状：`ENST00000450764.1` / `ENSG00000105576.9`，仅由 train 选择。
- SNP：640 个，20 blocks；缺失填补、AF、MAF 和 PCA 只由 train 拟合。
- Encoder：40 epochs masked-genotype pretraining；checkpoint 在提取 H 和 phenotype decoder 前保存并冻结。

## 各臂 held-out predictive R²

| Arm | R² | paired-bootstrap marginal 95% CI | alpha |
| --- | ---: | ---: | ---: |
| A | 0.000003 | [-0.000009, 0.000016] | 1e+06 |
| B | 0.022519 | [-0.001348, 0.044659] | 10000 |
| C | 0.015716 | [-0.047531, 0.072397] | 100000 |
| D | 0.033765 | [-0.100780, 0.147438] | 10000 |
| E | 0.017183 | [-0.108832, 0.101826] | 10000 |

## 解释边界

这是小样本真实分子表型 development Pilot，不是外部复制、全基因组有效性、临床预测或因果证据。
`E` 是 Grassmann 候选机制臂；无论它输赢，都不能改写唯一主问题 `C-B`。后续是否进入新数据或 biology-prior
阶段，只按本文件已经给出的冻结结果决定，不事后修改本次标准。
