# Grassmann 表型信号压缩判决 · Gate 0A 设计包 rc2

_状态：`DESIGN_ONLY_RUN_NOT_AUTHORIZED`；日期：2026-09-04；取代 rc1;只冻结设计,不授权运行。_

---

## 📋 rc2 相对 rc1 改了什么(以及为什么)

rc1 的主测试臂是 raw genotype 上的 RBF kernel-PCA——它**不是**你的预训练表示,也**不是** Grassmann,却把 kill 判据写成"决定 Stage 1 + Stage 2 / Grassmann"。这是过度归因。rc2 修正范围:

- 明确降级为 **Gate 0A:非线性压缩 headroom**,只回答"raw genotype 上有没有比 PCA 更保信号的廉价非线性映射";
- **删除** rc1 的 Grassmann 几何特异 kill(本 gate 没有 Grassmann 臂,判不了);
- 冻结六个 DGP regime,分开判,不做跨 regime 平均;
- 加 MAF 标准化臂、total-budget 口径、成本轴、`R²_genetic`、replicate-bootstrap 推断、detectability gate;
- 记录 0A→0B→1→2 gate ladder(仅 0A 在本包设计,其余不授权)。

## 🎯 这个 gate 判什么

> 在共享基因型 panel 上,**逐 DGP regime**,有没有一种廉价的非线性低维映射,在固定表示预算下,比匹配的线性重建 null(block PCA)保留更多**遗传**信号?

它**不能**判预训练 encoder(Gate 0B),**不能**判 Grassmann(Gate 1)。三者是不同命题,见 `GATE_LADDER_ROADMAP.md`。

## 🎯 为什么值得跑(机制风险)

无监督重建优化的是**基因型方差**,而表型因果方向没有理由与基因型 top 方差对齐——基因型 top PC 主要装的是祖源与 LD block 结构,不是某性状的因果方向。所以有很强的先验:除非非线性专门重排信号,否则谱压缩会把表型信号丢掉。Gate 0A 用 `spectral_tail` 这个对抗 regime 直接测这个失败模式——它是**主判决 regime**。

## 🎯 null 为什么是 block PCA(以及它的确切地位)

目标是无监督重建,PCA 是线性**重建最优**,所以是匹配的线性对照——但它**不是**"所有线性表型保留方法的上界"(监督式 PLS 是另一回事,且本包禁用)。主 null `B_pca_z` 与 kernel encoder 用同一套 MAF 标准化输入,消除预处理混淆;`B_pca_raw` 作预处理敏感性次 null。

## 🧪 对照臂(同 total-budget,输入/split/性状逐字节相同)

| 臂 | 表示 | 角色 |
| --- | --- | --- |
| `A_test`(主) | RBF kernel-PCA top-`k` | 主测试臂 |
| `A_test`(次) | 浅层 per-block autoencoder | 仅稳健性 |
| `B_pca_z` | MAF 标准化上的 block PCA top-`k` | **主线性 null** |
| `B_pca_raw` | raw-centered 上的 block PCA top-`k` | 预处理敏感性 null |
| `B_rand_bilinear` | 同维随机双线性 | **仅探索,无 kill** |
| `B_rand_proj` | 高斯随机投影 | 平凡参考 |
| `B_ldprune` | 预算匹配的 LD-pruned 子集 | 朴素选择参考 |
| `B_blockmean` | block 均值 / haplotype dosage | 地板 |

`B_rand_bilinear` 在 Gate 0A 只作探索,不携带 Grassmann kill;它只有到 Gate 1 才成为匹配对照。

## 📐 是 budget,不是 "k"

per-block `k` 的总维是 `D_total(k)=Σ_b min(k,m_b)`,不是 `k`。判决横轴是 **total representation budget** 与压缩率,per-block `k` 并列报告。主预算点是 per-block `k∈{8,16}`(运行前定死,避免 `argmax_k` 选择偏差);`{1,2,4,8,16,32,64}` 为探索。

## 🧬 六个 DGP regime(分开判)

见 `config/DGP_REGIMES.json`,只做定量性状,`h²∈{0.05,0.1,0.2,0.4}`,各带多基因背景:

1. `major_LD_aligned`——加性,因果在 leading-PC/强 LD(PCA 有利);
2. `spectral_tail`——加性,因果在低方差谱尾(**对抗;主判决 regime**);
3. `low_MAF`——加性,因果在低频变异;
4. `within_block_interaction`——block 内 `γ·G_i·G_j`,匹配交互 head;
5. `between_block_additive`——远 block `β_A G_A+β_B G_B`,**分布式多基因,非** long-range 交互;
6. `between_block_interaction`——跨 block `γ·G_A·G_B`,**真正的 long-range 测试**,匹配交互 head。

交互 regime 的下游是对表示做二次多项式展开后的 ridge,**所有臂完全相同**——线性 ridge 读不出任何表示里的交互,所以必须给交互能力,且相同以免偏袒某臂。丢掉交互因子的表示会在这里失败,这正是想要的信号。

## 📊 指标

- **主判决指标:** `R²_genetic(k)=R²(g, ĝ)`,`ĝ=下游(Z_k)` 拟合到**已知**遗传值 `g`;它把表示损失与表型噪声分开,并绕开 `p≫n` 的全基因型上限问题;
- **次(现实)指标:** `R²_pheno(k)=R²(y, ŷ)`;
- 第一轮**删掉** binary/AUC(未定义 penetrance/liability/prevalence/scale)。

## 🔬 推断(基于 replicate,不是 fold)

推断单元是**模拟 replicate(seed)**,不是 CV fold。每个 replicate `r`:`Δ_r=R²_genetic(A_test,r)−R²_genetic(B_pca_z,r)`;95% CI 用 replicate 上的 percentile bootstrap(10000 次)。CV fold 只用于惩罚项选择与拟合。

## 🔪 kill 判据(仅 Gate 0A)

在**每个** DGP regime **分别**、在主预算点上:若 `A_test` 对 `B_pca_z` 的 paired-mean `ΔR²_genetic` 的 replicate-bootstrap 95% CI **不排除零**,则该 regime 无非线性压缩 headroom。

> The paired-mean criterion is evaluated separately within each pre-registered DGP regime; no pooled cross-regime average is used for the primary decision.

程序含义:若连 `spectral_tail` 与 `between_block_interaction` 都失败,则无监督谱压缩前提错误,该路径停;通过的 regime 定义带往 Gate 0B 的 headroom。本 gate **没有** Grassmann kill,也**没有**笼统的"Stage 1+2 不存在"kill。

## 💰 成本轴(恢复)

记录 bytes/sample、峰值内存、encode 时间、下游 fit/predict 时间、复杂度 vs SNP 数。同等 retention 但预算/成本更低是赢;retention 只高一丁点却贵很多不算赢。恢复"retention × cost"二维判决。

## 🔍 设计包自检(不执行实验)

```bash
python scripts/validate_design_package.py
```

只校验设计包自洽与 manifest,输出 `DESIGN_PACKAGE_INTEGRITY_VALID`——**仅**代表文件包自洽,不代表科学设计已 valid。

## 🚧 运行前 gate(缺一不可)

1. 修好并验证**中心化与标准化**(仅用训练折统计量);
2. 绑定数据契约 hash;
3. 跑 **detectability gate** 定死 margin 与 replicate 数(h²=0.05、N≈2247 下 20 seeds 很可能不够);
4. 跑 compute smoke(3 block × 2 k × 1 fold),或预注册 randomized-SVD/Nyström/Lanczos 回退;
5. 取得项目负责人运行授权——然后才写运行代码。

## 🚫 不授权事项

不执行本包任何运行;不训练 Transformer;不设计或运行 Gate 0B/1/2;不做 phenotype gating;不做 GPU 或 v7/A1-R;不读 HGDP;不产生 binary/真实表型 claim;不产生架构 GO/NO-GO;不产生本 gate 的 Grassmann 特异 claim;主判决不跨 regime 平均。
