# TQ-G1：annotation-free Gene-layer Readout Pilot rc0（待冻结）

_起草：2026-09-17｜冻结：2026-09-17｜状态：`FROZEN_AWAITING_RUN`_

> 冻结发生在任何 chr18 outcome 被读取之前。被冻结文件的 SHA-256 见同目录 `FREEZE.sha256`；运行器会把实际使用的 protocol/config/core 哈希写进 `results/RUN_BINDING.json`，两者必须一致。
> 冻结之后不得修改本文件、`CONFIG_FROZEN.json` 或任何分析代码。需要改动就开 rc1，不覆盖 rc0。

## 唯一主问题

在 GEUVADIS chr18 上，把一个**共享的、冻结的 genotype-only encoder** 的 hidden state 按 **gene 为单位**、**只用基因组距离**池化之后，相对同一 cis 窗口的 additive dosage ridge，是否增加 held-out predictive R²？

主 estimand：

\[
\Delta R^2 = \frac{1}{G}\sum_{g=1}^{G} R^2_g(\texttt{C\_gene}) - \frac{1}{G}\sum_{g=1}^{G} R^2_g(\texttt{B})
\]

**gene 是复制单位。** 这是本 Pilot 相对 rc1/rc2 的核心方法学改动：0526 Stage 2 rc2 用 1 个 transcript、chr18 Stage 2B rc1 用 8 个 trait，两次的 CI 宽度都由 trait 数量而非个体数量主导。本次以数百个 gene 做 macro 平均，并对 individual × gene 做 two-way clustered bootstrap。

## 本 Pilot 不回答什么

- **不是 burden test。** 本目录不读取任何 WES、LoF、missense 或 rare-variant 数据，也不产生 causal effect。
- 不使用 eQTL、promoter/enhancer、chromatin contact 或任何 functional annotation 做 SNP→gene routing。routing 权重是 **TSS 距离的确定性函数**，理由见下节。
- 不检验 Grassmann；`E` 臂不在本 Pilot 中。既有 `NO_GO` 不因本 Pilot 改变。
- 不是外部复制，不是临床表型预测，不是全基因组结论。

## 为什么 routing 必须 annotation-free

若 `a_jg` 取自 eQTL / chromatin / functional annotation，则由模型导出的 gene effect score 与未来的 WES burden `β̂_g` 之间的任何相关，都可以被第三变量「该基因已被既有遗传学证据标注为重要」解释。那样检验的是注释，不是 encoder。

因此本 Pilot 只允许几何 routing：

\[
a_{jg} \propto \exp\!\left(-\frac{|\mathrm{pos}_j - \mathrm{TSS}_g|}{L}\right),\qquad L = 100{,}000\ \mathrm{bp}
\]

annotation-injected routing 是**另一个**预注册对照，只有在本 Pilot 通过后才允许写协议。

## 数据、角色与防火墙

- 资产与 `20260915_geuvadis_chr18_stage2b_rc1` 完全相同：GRCh38 chr18 PGEN/PVAR/PSAM（445 人、367,759 变异）与 445-sample normalized expression BED；哈希在 `PRE_RUN_BINDING.json` 中强制。
- **sample split 复用 rc1 的 `SAMPLE_MANIFEST.tsv`，不重新生成**，哈希 `7136c1c4…9875b`。家庭是不可拆分单位。
- `dev_train` 120：QC、AF、PC、encoder 预训练、ridge 拟合。
- `dev_validation` 30：ridge alpha 选择、encoder validation、目标量 `T_g`。
- `task_gate` 50：本 Pilot 的 held-out 评价集。
- `bridge_test` 230：**保持密封。** 运行器在配置指向密封角色时硬失败，除非同时给出 `--open-sealed-test` 与书面 `OPEN_AUTHORIZATION.json`。
- RC2 的 15 个 individuals 保持永久排除。
- expression 只按 `dev_train ∪ dev_validation ∪ task_gate` 的样本列读取；密封角色的 outcome 从不进入内存。

### 必须写明的限制

`task_gate` 的 50 个个体在 Stage 2B rc1 中已被打开过一次（用于 `B−A`）。因此它们**不是 virgin test set**。本 Pilot 的 estimand 不同（gene-macro `C_gene−B`），但这一点不能事后删除，必须与结果一同报告。

## 冻结的选择规则（都不看 phenotype）

1. QC：`dev_train` 上 MAF ≥ 0.05、missingness ≤ 0.10。
2. 候选 gene：chr18、expression 有限值、development 方差 > 0、cis ±1 Mb 内 QC 变异 ≥ 256。
3. cis SNP：在 cis 窗口内按坐标**均匀抽稀**到 256 个，再切成 8 × 32 的 local blocks。
   - 这是相对 rc1 的第二个有意改动：rc1 取 development 中单 SNP `|r|` 最大的 256 个，那是 phenotype 驱动的，会同时抬高 `B` 并让 `C−B` 更难解释。本 Pilot 用纯几何抽稀。
4. gene 上限：对合格 gene 按 `sha256("tqg1-rc0:" + gene_id)` 排序取前 200 个。这是确定性的、与信号无关的截断。
5. 合格 gene 少于 60 个则判 `DESIGN-INELIGIBLE` 并终止，不放宽规则。

## 共享 encoder

| 项 | 值 |
| --- | --- |
| 架构 | `LocalGenotypeEncoder`，d_model 16、2 层、4 heads、FFN 64、dropout 0 |
| 位置编码 | `relative_continuous`（整数 bp 导出的相对位置与 marker gap） |
| block lookup | **禁止**（`use_local_block_embedding=false`） |
| 预训练 panel | chr18 QC 变异均匀抽稀出的 512 blocks × 32 SNP = 16,384 SNP |
| 目标 | 20% masked genotype，三分类 CE |
| 优化 | AdamW、lr 0.002、weight decay 1e-4、batch 16、40 epochs |
| 输入 | 只有 `dev_train` genotype；不含任何 phenotype |

**一个 encoder 服务全部 gene**，不是 rc1 的每 trait 一个 encoder。没有 block lookup、没有 locus-specific 参数，这是「共享 local encoder」这一说法在本仓库 E0 审计后的唯一合法实现方式。

## 比较臂

| Arm | 输入 | 地位 |
| --- | --- | --- |
| A | sex + population one-hot + 10 个 `dev_train` 拟合的 genotype PC | 结构参照 |
| B | A + 256 个 additive cis dosage | **最强简单基线，主对照** |
| C_gene | A + 冻结 encoder 的 annotation-free gene 池化读出 | **唯一 primary 臂** |
| C_full | A + 256 × 16 坐标对齐 full `H` | rc1 参照臂 |

`C_gene` 的池化特征为三个几何视图拼接：cis 均匀均值、TSS 距离核加权均值、逐 block 均值。

Decoder 统一为 dual ridge；alpha grid `10^-6 … 10^6`，只由 `dev_validation` MSE 选择，平局取更大 alpha。R² 分母统一为 `Σ(y_eval - mean(y_train))²`。

## 判定

- **支持：** `ΔR² > 0` 且 two-way clustered bootstrap 95% CI 下界 `> 0`。
- **不支持：** 上述不成立。不得事后增删 gene、改 routing 核宽、改 alpha grid 或改 bootstrap 单位。
- Bootstrap：5,000 次，同时对 evaluation individual 与 gene 有放回重抽。
- 次要对比：`C_full−B`、`C_gene−C_full`、`B−A`，均报告，均不作为 primary 的前置 gate。

## Gene effect score：burden test 的**预演**，不是 burden test

对每个 gene 定义扰动分数

\[
s_g = \mathrm{sd}_{i \in \text{eval}}\left[\hat y_g(x_i) - \hat y_g(x_i^{\text{cis-}g\ \text{removed}})\right]
\]

- `B` 臂：cis dosage 列替换为 `dev_train` 均值 `2\hat p`。
- `C_gene` 臂：cis genotype token 全部替换为 encoder 的 MASK token 后重新编码。

两种干预不是同一个数学操作（一个在 dosage 空间，一个在 token 空间），这一点必须与结果同报。共同点是「把该 gene 的 cis genotype 信息拿掉」。协变量在两次预测中完全相同，故其贡献在差分中抵消。

目标量取 `T_g = dev_validation 上 B 臂的 held-out R²`，即该 gene 的 cis 可预测性。它与 `s_g` 使用**不相交的个体**。

预演统计量：

\[
\rho_{\text{model}} = \mathrm{Spearman}(s_g^{\texttt{C\_gene}}, T_g),\qquad
\rho_{\text{linear}} = \mathrm{Spearman}(s_g^{\texttt{B}}, T_g)
\]

**决定性的量是差** `ρ_model − ρ_linear`，以 gene 为单位做 paired bootstrap。若它不显著大于 0，则 foundation encoder 在 gene effect map 上相对纯 additive ridge 零增量，后续 UKB + WES burden 的大计划不应启动。

Matched null：在 `cis SNP 数十分位 × 平均 MAF 十分位 × development 表达方差十分位` 构成的层内置换 `s_g`，5,000 次，报告单侧与双侧 p。

### 迁移到真正的 burden test 时只换一处

把 `T_g` 从「dev_validation 的 cis 可预测性」换成独立 WES 的 `β̂_g`，估计量、matched null 与 `ρ_model − ρ_linear` 的判据全部不变。本 Pilot 的作用就是先证明这套估计量在一个已知有信号的分子表型上能工作。

## 已知限制（不可事后删除）

1. 不是 burden test，不产生 causal effect。
2. `task_gate` 个体在 rc1 中已被打开一次。
3. 供应的 expression 是在完整 445 cohort 上预处理过的 outcome，因此是 individual holdout，不是 preprocessing-sealed replication。
4. `T_g` 的 alpha 由 `dev_validation` 选出，故 `T_g` 含轻微乐观偏差；该偏差对 model 与 linear 两个 `ρ` 对称，因此差值仍可解释。
5. 两臂的扰动干预不是同一个操作。
6. chr18 单染色体、单一细胞类型（LCL），不外推全基因组或其他组织。

## 失败只停止什么

`C_gene−B` 不支持，只停止「在本设置下、annotation-free 几何 routing 的 gene 池化优于同窗口 additive dosage」这一条。它不停止：SNP foundation modelling、E0 masked-genotype gate、annotation-injected gene layer、或在更大 cohort 上的重做。
