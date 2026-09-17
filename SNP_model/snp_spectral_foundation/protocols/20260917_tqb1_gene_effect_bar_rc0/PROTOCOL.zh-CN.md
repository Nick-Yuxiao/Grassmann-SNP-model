# TQ-B1：additive gene-effect bar（rc0,待冻结）

_起草：2026-09-17｜状态：`DRAFT_AWAITING_FREEZE`_

## 唯一主问题

在一个真实 cohort 上,**纯加性** cis 模型的 gene 扰动分数 `s_g`,对独立 WES burden 效应 `|β̂_g|` 的秩预测能力有多强?

\[
\mathrm{bar} = \mathrm{Spearman}(s_g,\ |\hat\beta_g|)
\]

这个数就是杆。任何 foundation model 在 gene effect map 上的分数,必须超过它才算有贡献。

**先量杆,再造东西。** 这是 TQ-G1 rc0 的直接教训。

## 不做什么

- 不训练任何神经网络,不使用 GPU,不依赖 torch。
- 不计算 burden 统计量。`β̂_g` 必须由**独立的 WES 分析**提供,且该分析不得使用本协议所用的芯片数据。
- 不做因果声明。`s_g` 是幅度,不是有符号效应,因此对照 `|β̂_g|`。
- 亲缘关系只通过用户提供的 `group` 列处理;不提供则视全体无关。

## 冻结的两阶段顺序

### 阶段 1：可检测性 Gate（在 validation 上）

\[
\frac{1}{G}\sum_g \left[R^2_g(B) - R^2_g(A)\right] > 0,\quad \text{95\% CI 下界} > 0
\]

- `A`：仅协变量。
- `B`：协变量 + 该基因 cis 窗口的 256 个 additive dosage。
- Bootstrap 同时重抽 individual 与 gene。

**不过的 trait,其 evaluation split 永不打开。** 这是硬 Gate,不是参考指标。

> TQ-G1 rc0 把同一个量降级为"次要对比",结果整轮不可解释。此处不重复该错误。

### 阶段 2：杆（在 evaluation 上,仅限通过 Gate 的 trait）

\[
s_g = \mathrm{sd}_{i\in\text{eval}}\left[\hat y(x_i) - \hat y(x_i^{\text{cis-}g\ \text{置为训练均值}})\right]
\]

协变量在两次预测中完全相同,其贡献在差分中精确抵消。

Matched null：在 `cis SNP 数十分位 × 平均 MAF 十分位 × burden se 十分位` 构成的层内置换 `s_g`。
加入 burden se 是必要的——它是 carrier 数的代理,否则"罕见变异多的基因"会同时驱动两侧。

## 两个必须保持的估计量约束

**1. 协变量不受惩罚。**

若协变量与基因型共用一个 alpha,加入无信息的 dosage 列会逼高 alpha、压扁协变量系数,使 `B` 结构性地劣于 `A`——与有无遗传信号无关。本包对协变量块施加零惩罚,于是 `alpha → ∞` 时 `B` 精确退化为 `A`,`B` 永远不会结构性输给 `A`。

该缺陷存在于 TQ-G1 rc0,见 `../../pilots/20260917_tqg1_gene_layer_readout_rc0/IMPLEMENTATION_INCIDENT_02.md`。

**2. 输出全精度。**

`GENE_TABLE.tsv` 用 `repr` 级精度写出。6 位小数会在秩相关上制造并列,使独立复算失败,见该目录的 `IMPLEMENTATION_INCIDENT_01.md`。validator 会显式检查浮点列的去重计数等于基因数。

## 不看 phenotype 的选择规则

1. QC：训练集上 MAF ≥ 0.01、missingness ≤ 0.05。
2. 合格基因：cis ±1 Mb 内 QC 变异 ≥ 256。
3. cis SNP：窗口内按坐标**均匀抽稀**到 256 个。不按相关性挑。
4. 基因上限：按 `sha256("tqb1-rc0:" + gene_id)` 排序截断。
5. 合格基因少于 200 判 `DESIGN-INELIGIBLE`。

## 判读

| bar | 含义 |
| --- | --- |
| ≈ 0.8 | 加性方法已还原 burden 排序;foundation model 空间极小 |
| ≈ 0.3 | 加性方法在 gene 层分辨率有限;有真实空间 |
| matched-null p 不显著 | 该排序可被 cis SNP 数/MAF/burden 精度单独解释,杆不成立 |

## 失败只停止什么

Gate 不过,只说明**该 trait 在该 cohort、该 cis 半径、该协变量集合下不可检测**。不停止 foundation modelling,不构成对 gene layer 的任何结论。
