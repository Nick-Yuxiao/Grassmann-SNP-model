# TQ1 / E0 parallel readiness rc1（冻结）

_冻结日期：2026-09-15。本文在读取原 Stage 2B `bridge_test` phenotype、或计算新的 E0 holdout 指标之前写入。_

## 决策结构

```text
H:TQ1：仅 A/B phenotype task qualification ─┐
                                             ├─ 两者 PASS → 允许设计全新 Bridge
E0：仅 genotype masked pretraining audit ────┘
```

两条证据线相互独立。TQ1 不读取 checkpoint、hidden state、C/D/E 或 Grassmann；E0 不读取 phenotype。任一失败都只停止“用当前 task + 当前 local encoder family 进入新 Bridge”，不停止 SNP foundation modelling，也不把学习 MAF、LD、haplotype 或 population structure 定义为失败。

## H:TQ1

- 目标 estimand：固定 8 个 chr18 expression traits 上，230 名原 `bridge_test` 个体的 macro `R²(B)-R²(A)`。
- A：sex + population one-hot + 10 个仅由 `dev_train` genotype 拟合的 PCs。
- B：A + 同一 trait 已冻结的 256 个 additive dosages。
- trait、SNP、样本、ridge alpha grid、development decoder、R² denominator 与旧 Stage 2B 完全相同。
- primary PASS：`B-A > 0` 且 paired individual-bootstrap 95% CI lower > 0；保持旧 Task Gate 标准，不事后提高门槛。
- family-cluster bootstrap 仅为预先指定的敏感性分析，不改变 primary 判定。
- 读取这 230 人 phenotype 后，该 pool 永久标记为 `CONSUMED_BY_TQ1`，其任何既存 C/D/E prediction 均不得用于 Bridge 推断。

TQ1 的 positive task 是 cis-expression prediction；最强简单 baseline 是 A；保留 additive dosage，控制 sex、population 和 genotype PCs。TQ1 FAIL 仅说明该固定 task panel 仍不够稳定，不能说明 encoder 或 Grassmann 失败。

## E0

- 目标 estimand：固定的 8 个 Stage 2B local encoders，在 development 之外 280 名个体上是否利用同一个人的局部 genotype context 改善 masked-genotype prediction。
- 训练和模型选择：禁止；只加载既存 40-epoch checkpoints。
- 输入：同一 8×256 SNP panels；每 panel 使用 5 个固定 20% masks。
- 最强 marginal baseline：仅由 `dev_train` 估计、Jeffreys smoothing (`+0.5`) 的逐位点 empirical genotype probabilities。
- context destruction：在 population 内固定循环错配 donor individual，保留 locus identity、AF、mask 与 ancestry composition，但用 donor 的非遮蔽 genotype context 预测原 individual 的 masked targets。
- co-primary estimands：
  1. `CE(empirical marginal)-CE(real-context model)`；
  2. `CE(donor context)-CE(real context)`。
- E0 PASS：两个 macro contrast 均 > 0，且 family-cluster paired-bootstrap 95% CI lower 均 > 0。
- accuracy、逐 trait 数值和旧 validation contextual lift 全部报告，但不替代 co-primary 判定。

E0 只验证这 8 个 phenotype-blind local encoder/checkpoint 在这些 locus panels 上的 genotype objective；它不证明 genome-wide foundation model、phenotype relevance、Grassmann utility 或跨 cohort transport。

## 合并裁决

- `READY_FOR_NEW_BRIDGE`：TQ1 PASS 且 E0 PASS。
- `NOT_READY_TASK`：TQ1 FAIL，不论 E0 结果。
- `NOT_READY_ENCODER`：TQ1 PASS 但 E0 FAIL。
- `NOT_READY_BOTH`：二者均 FAIL。

即使 `READY_FOR_NEW_BRIDGE`，也只能新建设计并绑定未参与本轮的 phenotype individuals/cohort；不得重新开启已经被 TQ1 消耗的 GEUVADIS 230 人 C/D/E。

