# TQ2 A/B empirical power calibration rc1（冻结）

_冻结日期：2026-09-15；本实验不读取任何真实 phenotype。_

## 问题

在 Stage 2B/TQ1 的真实 genotype、固定 8×256 SNP panels、120/30/230 family-disjoint split、A/B features、ridge alpha grid 和 Task PASS rule 下，当真实 additive genotype increment 已知时，A/B pipeline 是否具有足够的经验检出功效？

这是一项设计校准，不是新的 Task Qualification，不验证 encoder、hidden state 或 Grassmann。

## 固定设计

- A：sex + population one-hot + 10 个 train-only genotype PCs。
- B：A + 256 additive dosages。
- 每个 trait/seed 使用 4 个 panel 内 causal SNP；权重、covariate score 和 noise 由固定 seed 生成。
- genetic score 在 train 中对 A residualize，避免把 ancestry/covariate signal 误计为 genotype increment。
- covariate variance fraction：0.20。
- 注入增量网格：`0, 0.01, 0.02, 0.05, 0.10`。
- 每个增量 100 个独立 phenotype seeds；每个 seed 含 8 个 traits。
- decoder、validation alpha selection、test R² denominator 与 Stage 2B 相同。
- 每个 seed 用 1000 次 paired individual bootstrap；PASS rule 仍为 macro `B-A>0` 且 CI95 lower > 0。

## 设计裁决

- null false-positive rate 必须 ≤0.05；
- 最小科学增量 `0.02` 的 empirical power 必须 ≥0.80；
- 两者同时满足才记为 `AB_DESIGN_POWER_READY`。

若失败，说明现有 120-person decoder training design 不能可靠检验原先声明的最小增量；它不会改写 TQ1 的阴性结果，但会阻止在同样本规模上继续寻找 phenotype task。下一版必须优先增大 decoder training N、减少无根据的 feature burden，或使用外部冻结效应，而不是调 encoder/Grassmann。

