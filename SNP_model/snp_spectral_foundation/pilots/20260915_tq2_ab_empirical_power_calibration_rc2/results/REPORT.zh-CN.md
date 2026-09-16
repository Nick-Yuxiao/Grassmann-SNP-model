# TQ2 A/B empirical power calibration 结果

## 裁决

最终状态：**AB_DESIGN_UNDERPOWERED**。停止在当前 GEUVADIS 规模继续修改 Task Gate；不启动真实 TQ2，不打开新 Bridge。

## 为什么做这项校准

Stage 2B 的规划功效 `0.8159` 没有把高维 decoder 拟合、alpha 选择和训练样本波动纳入计算。TQ2 使用真实 chr18 genotype 并注入已知 additive signal，直接检验完整 A/B pipeline 的实际错误率和功效。全过程未读取真实 phenotype，也未使用 encoder、H 或 Grassmann。

## rc1：原 120/30/230 设计

| 注入 ΔR² | A/B Gate 通过率 |
|---:|---:|
| 0.00 | 0.37 |
| 0.01 | 0.36 |
| 0.02 | 0.39 |
| 0.05 | 0.35 |
| 0.10 | 0.62 |

零效应通过率 37%，说明 individual test bootstrap 没有覆盖 decoder fitting/validation selection 的波动，并且 B 可以通过 dosage 重复表达 A 中的结构。它既反校准，也没有达到 `ΔR²=0.02` 的 80% power。

## rc2：family cross-fit + residual dosage

| 注入 ΔR² | oracle mean ΔR² | A/B Gate 通过率 |
|---:|---:|---:|
| 0.00 | 0.0000 | 0.00 |
| 0.01 | 0.0103 | 0.00 |
| 0.02 | 0.0198 | 0.00 |
| 0.05 | 0.0495 | 0.00 |
| 0.10 | 0.0996 | 0.01 |

rc2 消除了零效应假阳性，但在 445 人、每 trait 256 dosage features 的条件下，residual ridge 无法稳定学习即使 oracle 已存在的增量。它不是合格的替代 Gate。

## 综合解释

1. E0 PASS 仍然成立：local encoder 能利用真实 genotype context。
2. TQ1 的 `B-A=-0.01745` 仍是合格的阴性结果，但不能被升级为“expression 没有 genotype signal”或“encoder 无用”。
3. 更关键的新发现是：当前 GEUVADIS 样本规模无法为高维 A/B Task Qualification 同时提供可靠的零效应校准和 `ΔR²=0.02` 检出功效。
4. 继续换 bootstrap、alpha 或 residual 规则会进入 outcome-independent 但仍然无止境的 protocol tuning，因此在 rc2 后停止。

## 下一步约束

下一次 Task Qualification 必须绑定更大的 individual-level genotype + quantitative phenotype cohort，并在 outcome 访问前完成：样本权限、family/relatedness、assembly/allele、协变量、train/validation/test 或 cross-fit 设计以及 empirical null/power calibration。当前服务器 UKB genotype inventory 存在，但冻结 manifest 明确写有 `phenotype_access_permitted=false`，因此在合法 phenotype 授权与路径绑定前不能运行 UKB TQ。

新 Bridge 仍需满足：新 TQ PASS + E0/对应新 encoder genotype Gate PASS；Grassmann只作为可选 secondary arm。

