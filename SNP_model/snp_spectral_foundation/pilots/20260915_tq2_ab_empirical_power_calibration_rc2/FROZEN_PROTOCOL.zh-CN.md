# TQ2 A/B empirical power calibration rc2（冻结）

_冻结日期：2026-09-15；rc1 已永久保留为 `AB_DESIGN_UNDERPOWERED`。本实验不读取真实 phenotype。_

## rc1 暴露的问题

旧 120/30/230 pipeline 在零 genotype increment 下的经验通过率为 0.37，在注入 `ΔR²=0.02` 时仅为 0.39。individual test bootstrap 没有覆盖 decoder fitting 与 30-person validation alpha selection 的不稳定性；同时 B 中的 dosage 可重复表达 A 的 ancestry/covariate structure。

## rc2 仅修正 Task Qualification 设计

- 使用全部 445 人的 genotype，按 family 分组、population 分层做 5-fold out-of-fold prediction；真实 TQ 使用此设计后整批数据被消耗，Bridge 必须使用外部新 cohort。
- 每个 outer training fold 内重新估计 AF/imputation、10 PCs 和所有标准化。
- A 使用含截距的 OLS：sex + population one-hot + 10 train-only genotype PCs。
- dosage 的每一列先在 outer train 中对 A 做线性投影，B 只允许读取该 residual dosage。
- B prediction = A prediction + residual-dosage ridge prediction。
- ridge alpha 由 outer-train generalized leave-one-out criterion 选择，不建立小型 30-person validation split。
- 其余 generator 保持 rc1：8 traits、每 trait 4 causal SNP、covariate variance 0.20、100 phenotype seeds、增量网格 `0/0.01/0.02/0.05/0.10`。
- 每个 seed 使用 1000 次 family-cluster bootstrap；PASS 为 macro `B-A>0` 且 CI95 lower >0。

## 冻结裁决

- null false-positive rate ≤0.05；
- `ΔR²=0.02` empirical power ≥0.80；
- 同时满足才是 `AB_DESIGN_POWER_READY`。

本修正只服务于“dosage 是否在协变量之外提供可泛化 phenotype signal”这一 estimand，不改变 E0，不评价 representation 或 Grassmann。若仍失败，停止在 GEUVADIS 规模继续修 Task Gate，转向更大 individual-level cohort。

