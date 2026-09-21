# G2PT 论文规格抽取（bioRxiv v2，2025-04-11）

来源：`paper/G2PT_biorxiv_v2_20250411.pdf`（用户提供），DOI `10.1101/2024.10.23.619940`，23 页。
抽取范围：Materials and Methods（p13–17）、Results 中的数量（p5–p8）、Abstract（p2）。
**所有条目均已从原文读出，不再是检索摘要。**

> ⚠️ **版本警告：本文件抽自 v2，最新是 v3（2026-01-15）。** v2 的实验只有 **TG/HDL 一个性状**；v3 摘要另称覆盖 LDL、T2D 与跨人群外推。**若复现目标包含 LDL/T2D/跨人群，v2 不够用，必须补 v3。**

---

## 1. 数据与队列

| 项 | 原文取值 |
|----|---------|
| 队列 | UK Biobank，**Caucasian ancestry** 参与者 |
| UKB 申请号 | **51436** 与 **26041**（作者的；我们必须用自己的） |
| 样本量 | **423,888** 名参与者（有 TG/HDL 与配套基因型） |
| 基因型覆盖 | **203,126** 个 SNP（SNP array） |
| 基因型编码 | 0 = 纯合主等位（参考），1 = 杂合，2 = 纯合次等位 |
| 目标性状 | **log₂(TG/HDL)**，TG 与 HDL 单位 mmol/L |
| 协变量 | **sex, age, PC1–PC10**（m = 12） |

## 2. SNP 筛选（这是 G2PT 的输入特征定义）

- 方法：**BOLT-LMM** 做关联分析，因变量 log₂(TG/HDL)，协变量 sex + age + top 10 PC。
  （原文写作 "Bayesian logistic regression analysis using BOLT-LMM"，但因变量是连续量，实际应为 BOLT-LMM 的线性混合模型；照抄原文用词即可，实现按连续性状走。）
- 阈值：**p = 10⁻⁵ 到 10⁻⁸** 的一系列阈值，**阈值本身是被网格搜索的超参**。
- 规模：**p=10⁻⁸ → 2,429 个 SNP；p=10⁻⁵ → 5,510 个 SNP**。
- **关键：SNP 筛选只在 training set 上做**（见 §5 嵌套交叉验证），不能在全量上先筛再分割，否则信息泄漏。

## 3. SNP → 基因映射

三条证据取**并集**（一个 SNP 可映射到多个基因）：

1. **cS2G** —— 2023 年 11 月取自 `https://alkesgroup.broadinstitute.org/cS2G`
2. **GTEx v7 eQTL** —— `https://www.gtexportal.org/home/downloads/adult-gtex/qtl`，限定 7 个组织：
   adipose subcutaneous、adipose visceral omentum、liver、pancreas、adrenal gland、muscle skeletal、uterus
3. **最近基因**（基因组坐标，**GRCh37 / hg19**）

## 4. 本体（GO 层次）构建 —— B3 至此解除

- 基础：**Gene Ontology Biological Process，version 2023-07-27**
  下载地址原文给定：`https://release.geneontology.org/2023-07-27/index.html`
- 剪枝：用 **DDOT** 包，只保留与"有显著 SNP 映射的基因"相关的 term
- 两条排除规则：
  1. 映射基因数 **< 5** 的 system 剔除
  2. 注释基因集合与某个子系统**完全相同**的 system 剔除
- 最终图：三类节点（SNP / gene / system），三类有向边（SNP→gene、gene→system、subsystem→supersystem，后者取 GO 的 `is_a` 与 `part_of`）

> **这条直接否定了上游 notebook 的示范路径。** `Collapse_Gene_Ontology_Based_on_GWAS_results.ipynb` 里下载的是 GWAS Catalog 的 `GCST90257283`，而论文用的是**作者自己在 UKB 上跑的 BOLT-LMM 结果**。notebook 只是流程示范，不是论文路径。

## 5. 模型结构

| 项 | 原文取值 |
|----|---------|
| 嵌入维度 d | **64** |
| 传播阶段注意力头数 | **n = 4** |
| 翻译阶段（participant 节点） | **Differential Attention，n = 1 头**（为了可解释性） |
| 前向传播 | (1) SNP→gene，(2) gene→system，(3) subsystem→supersystem |
| 反向传播 | (4) system→subsystem，(5) system→gene |
| 合子性 | 按 z ∈ {hetero, homominor} 使用**不同的 HiGT 模块**；纯合主等位视为参考态，不参与 |
| 每步权重 | 五个步骤各自是独立层，**权重分开学习** |
| 协变量通路 | Pcov=(sex, age, PC1..PC10) 经 MLP 投到 d 维得 Pembed，再用 HiGT 以 genes+systems 为邻居更新 |
| 输出 | 各模块嵌入拼接后过最终层，预测 log₂(TG/HDL) |

**与代码默认值的对照（重要）：** 论文 d=64 对应 `--hidden-dims 64`，即上游 `train_model.sh` 里那个值，**不是**代码默认的 256。`--n-heads 4` 与代码默认一致。

## 6. 训练与评估

| 项 | 原文取值 |
|----|---------|
| 评估框架 | **嵌套交叉验证**（nested CV） |
| 划分比例 | train : val : test = **3 : 1 : 1** |
| 各集用途 | train 用于**筛 SNP + 拟合参数**；val 用于**网格搜索调超参**；test 为留出集 |
| 损失 | **MSE**（真实 vs 预测 TG/HDL） |
| 优化器 | **AdamW**（decoupled weight decay + L2 正则） |
| 网格搜索范围（G2PT） | **p-value 阈值** 与 **训练 epoch 数** |
| 对照模型 | ElasticNet、XGBoost、PRS C+T、PRS T、**LDPred2**、**Lassosum** |
| 对照模型的网格 | XGBoost: n_estimators / max_depth / subsample / lr；ElasticNet: alpha / L1 ratio / tol；LDPred2 与 Lassosum 用 training 群体的 LD 信息校正效应量 |

> **对 B9 的重要修正：论文的模型选择靠"对 epoch 做网格搜索"，不是靠代码里的 EarlyStopping。**
> 这意味着上游那个失效的早停/`.best` 逻辑**大概率没有影响论文结果**——作者按验证集在若干 epoch 检查点里挑。
> 但对我们的要求反而更明确：**必须保留逐 epoch 检查点（`{out}.pt.N`）并自己在验证集上选，绝不能用 `.best`。**

## 7. 可解释性与上位效应

**重要性打分：** 用训练好的模型，计算每个 gene/system 的个体注意力值与预测 TG/HDL 的 **Pearson 相关**，**男女分开算**，再取**绝对值的平均**作为重要性分数。
（这与上游 `sys_importance.csv` 的列名 `corr / corr_female / corr_male / corr_mean_abs` 完全对得上——可作为实现一致性的交叉验证点。）

**上位效应搜索（Fig. 4a）：**
1. 取重要性 **top 20** 的 system；
2. 对每个 system，取所有映射到该 system 及其子节点注释基因的 SNP；
3. **卡方检验**筛选"对该 system 注意力高 vs 低（按 system attention 排序的 **top 10% vs bottom 10%**）的个体间频率有显著差异"的 SNP；
4. 对 SNP 对做组合线性模型，检验交互项 α₁₂ ≠ 0，**Wald test + Benjamini-Hochberg 校正，调整后 p < 0.05**；
5. **排除相距 < 1 Mb 的 SNP 对**（避免 LD 混杂）。

## 8. v2 的主要数量结果（复现判定靶子）

- 注意力显著覆盖 **20 个 multigenic system**，共 **1,395 个 SNP**、**253 个基因**
- 1,395 个 SNP 中 **252 个**存在多重 SNP→基因映射，且 G2PT 把注意力显著集中在其中一个基因上（>2 倍差距）
- 检出 **40 个上位交互**，含 **APOA4–CETP**（phospholipid transfer）
- 在 10⁻⁸ ~ 10⁻⁵ 全部阈值上，R² 显著高于 PRS 与正则化回归

## 9. 仍未解除的条目

| 编号 | 仍缺什么 |
|------|---------|
| B6 残留 | **学习率、weight decay 具体值、batch size、dropout、epoch 网格的具体取值范围、nested CV 的折数**——Methods 只说"grid search"，未给网格点。需查 Supplementary，或按 val 集自行搜索并记为偏离 |
| B9 残留 | 论文是否传 `--target-phenotype` 仍不可知；但因论文按 epoch 网格选模型，此问题降级为"我们自己不要用 `.best`" |
| B2 | 无随机种子问题不受影响，且因 nested CV × 网格搜索而被放大 |
| 版本 | v3 新增的 LDL / T2D / 跨人群结果，v2 没有 |
