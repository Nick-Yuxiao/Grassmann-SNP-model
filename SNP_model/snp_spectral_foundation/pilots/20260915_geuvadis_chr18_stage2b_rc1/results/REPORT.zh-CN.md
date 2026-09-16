# 0526 Stage 2B GEUVADIS chr18 rc1 结果

## 冻结结论

最终状态：**TASK-INELIGIBLE**。

Task `B-A` = **0.010475**，paired 95% CI = **[-0.03848765353856709, 0.059945794164239376]**；判定 **TASK-INELIGIBLE**。

Bridge phenotype 未打开。

## 审计摘要

- Development / Task / Bridge / RC2 excluded = 150 / 50 / 230 / 15；家庭零交叉。
- traits：ENSG00000263753.7, ENSG00000242550.5, ENSG00000263006.6, ENSG00000270112.3, ENSG00000141384.12, ENSG00000179981.10, ENSG00000198796.6, ENSG00000266850.1。
- 设计灵敏度：预设 ΔR²=0.02，Monte Carlo power=0.8159，要求 ≥0.80。
- A-E 使用每个 trait 完全相同的 256 个 cis SNP、样本、covariates、decoder family 和一次性 test opening。
- RC2 原结论保持：`0526 Stage 2 bridge, RC2: NO-GO in this setting.`

## 解释边界

本结果只裁决本冻结 setting 的 **Task eligibility**；由于 Task Gate 未通过，它没有裁决 full-H bridge。它不把 MAF/LD/haplotype 学习判为失败，也不裁决整个 SNP foundation-model programme。发布的 normalized expression 在分组前已由上游处理，因此这不是 raw-expression preprocessing-sealed replication。
