# 0526 Stage 2B Bridge Replication Gate rc1（冻结）

_冻结时间：2026-09-15；冻结发生在任何 chr18 phenotype 筛选或 Task/Bridge outcome 分析之前。_

## 唯一科学问题

在独立 Task Gate 已证明 `genotype → expression` 可检测之后，冻结的 genotype-pretrained full hidden state `H` 是否比同一 SNP panel 的 additive dosage 提供额外 held-out predictive R²？

- Task estimand：8 个 trait 的 macro mean `R²(B)-R²(A)`。
- Primary bridge estimand：8 个 trait 的 macro mean `R²(C)-R²(B)`。
- RC2 固定保留为 `0526 Stage 2 bridge, RC2: NO-GO in this setting.`
- 本次失败最多终止 generic `H > dosage` bridge claim；不终止 SNP foundation modelling，不把 LD/MAF 学习称为失败。

## 数据、防火墙与坐标

- 数据：Broad TensorQTL 官方 example 中的 GEUVADIS 445-sample normalized expression BED 与 GRCh38 chr18 PGEN/PVAR/PSAM。
- expression BED：0-based half-open；每行 `[start,end)` 为一个 TSS point，分析用 `TSS_1based=start+1`。
- PVAR：GRCh38、`chr18`、VCF/PVAR 1-based `POS`；cis 定义为 `|POS-TSS_1based| <= 1,000,000`。
- RC2 15 个 test individuals 永久排除。家庭是不可拆分分配单元。
- 冻结样本数：development 150（dev_train 120 / dev_validation 30）、Task Gate 50、Bridge Test 230、RC2 excluded 15；五个人群均分层进入三个新分区。
- phenotype-derived TensorQTL covariates不进入 A；A 只含 sex、population one-hot 与按 trait panel 在 dev_train 拟合的 10 个 genotype PCs。
- supplied normalized expression 是发布前已在完整 445 cohort 上处理的 outcome；因此本研究是独立个体 holdout replication，但不是 raw-expression preprocessing-sealed replication。这个限制不能事后删除。

## Development-only trait 与 SNP panel 冻结算法

1. 只查看 development phenotype；Task 与 Bridge phenotype 在此阶段不进入内存。
2. 候选为 chr18、有限值且 development 方差大于 0 的 gene rows。
3. 每个 gene 的 ±1 Mb cis SNP 在 dev_train 上要求 MAF ≥ 0.05、missingness ≤ 0.10。
4. 对每个候选 gene，以 development 中最大单 SNP `|Pearson r|` 为完全确定的筛选分数。
5. 按分数降序、gene ID 升序贪心选择 8 个 trait；任意两者 TSS 至少相隔 2 Mb，且 development 中经 sex/population 回归后的 phenotype residual `|r| ≤ 0.20`。
6. 每个 trait 选择 development 单 SNP `|r|` 最大的 256 个 cis SNP，再按坐标排序为 8 × 32 的 local blocks。B–E 对该 trait 使用完全相同的 SNP identity 与顺序。

这些步骤可以提高正向 task 的信号，但不能制造确认性结果，因为 Task Gate 与 Bridge Test 个体均未参与选择。

## 冻结模型

| Arm | 输入 | 地位 |
| --- | --- | --- |
| A | sex + population + 10 train-fitted genotype PCs | covariate/structure reference |
| B | A + 256 additive dosages | 最强简单 task baseline |
| C | A + 256 × 16 coordinate-aligned full `H` | 唯一 primary bridge arm |
| D | A + rank-4 aligned low-rank `H` | RC2 登记的 secondary hypothesis |
| E | A + 每 block 的 rank-4 Grassmann projector + spectrum | secondary mechanism arm |

Encoder 继承 RC2：d_model 16、2 local layers、4 heads、FF 64、dropout 0、20% masked genotype、AdamW、lr 0.002、weight decay 0.0001、batch 16、40 epochs、最终 epoch checkpoint。Encoder 只见 dev_train/dev_validation genotype，不见 phenotype。Decoder 均为同一 dual ridge alpha grid；alpha 只由 dev_validation MSE 选择，tie 取较大 alpha。

## Task Gate 与 Bridge Gate

预测 R² 的 denominator 使用 dev_train phenotype mean。5000 次 paired individual bootstrap 中，每次对同一批 individuals 重采样，并对全部 traits/arms 共用该索引；trait 不伪装成独立样本。

- Task PASS：macro `B-A > 0` 且 paired bootstrap 95% CI lower > 0。
- Task FAIL：记录 `TASK-INELIGIBLE`，Bridge phenotype 不打开；不作 bridge NO-GO。
- Task PASS 后只允许打开一次 Bridge Test。
- Bridge PASS：macro `C-B > 0` 且 paired bootstrap 95% CI lower > 0。
- Bridge FAIL：`Stage 2B generic full-H bridge: NO-GO under this high-signal replication setting`。
- D-B、E-B 与各 trait 数值必须报告，但不能改写 C-B。

运行中不得改变阈值、trait 数、SNP 数、rank、decoder、alpha grid、bootstrap 或 representation。任何异常只能修代码并产生新审计记录，不能用 outcome 调参。

## 设计灵敏度

运行前最小科学效应固定为 macro ΔR² = 0.02；目标单侧 α=0.025、power ≥0.80。Development-only trait 选定后，用冻结的 K=8、Bridge n=230、B R²=0.10、trait residual exchangeable correlation 取选中 traits 的观测上界（至少 0.20）进行 20,000 次 Gaussian Monte Carlo。若 estimated power <0.80，Task/Bridge 均不打开，状态为 `DESIGN-INELIGIBLE`；不得事后把最小效应放大。

## 结果解释边界

这是一项 biological-effect bridge 的分子表型检验，不是因果推断、临床预测、全基因组证明或跨研究 expression replication。MAF/LD/haplotype 是 encoder 的合法 genotype 信号；ancestry 在这里由分层分配与 population/PC covariates控制，以免把 cohort structure 当 phenotype gain。

