# 20260918_complementarity_v1_rc0 — Concept Gate v1

一句话：**在已经冻结的 panel 上，用最便宜的确定性特征，测"功能注释"和"人群结构"
互相之间还有没有增量。不训练任何模型，不打开任何封存 split，不导出任何个体级数据。**

## 回答什么

```
delta F | P = R2(PF) - R2(P)    有了人群结构以后，功能注释还有没有增量
delta P | F = R2(PF) - R2(F)    反过来
```

## 不回答什么

- self-supervised pretraining 能不能降低 phenotype sample complexity（这需要 encoder，本包没有）；
- 任何 learned representation 是否胜过 raw dosage；
- Grassmann 几何的任何部分；
- 任何 confirmatory 结论（见下）。

## 五个 arm

`A` 协变量 · `B` A + 冻结阈值的 dosage score · `P` A + 分块局部 PCA ·
`F` A + 注释×dosage 聚合 · `PF` A + P + F

同一个 ridge 解码器，同一套 λ 网格，λ 一律由 **train 内部 K 折 CV** 选出。

## 两条设计红线

1. **P 必须是分块 PCA，不能是全局 PCA。** arm A 里已经含 PC1–40，全局 PCA 会和它
   共线，`delta P | F` 被构造性地压到 0。`c1` 会把 P 与 PC1–40 的 canonical
   correlation 报出来，让重叠是被量出来的而不是被假设的。
2. **F 的注释必须 phenotype-free。** ChromHMM / ATAC / 保守性 / 基因结构可以；
   GWAS Catalog、同性状 eQTL、已有 PRS 权重不行。`c2` 对未声明 `phenotype_free=true`
   的 track 直接退出，`c0` 对形似关联产物的名字直接拒绝。

## 结果是 developmental，不是 confirmatory

validation 这条 split 已经被用了两次：Task Qualification 用它挑 B 的阈值，
Concept Gate 用它做评估。所以 PASS = "值得写 confirmatory 协议"，不是 "已证实"。
`c3` 里没有打开封存 split 的开关，`tq_test` / `bridge_holdout` 在建设计矩阵之前
就被剔除。

## 文件

```
PROTOCOL.zh-CN.md              完整协议：arm 定义、解码器、判定规则、输出边界
ANNOTATION_SOURCES.zh-CN.md    F 的注释从哪来、为什么 phenotype-free、怎么下
SERVER_STEPS.zh-CN.md          服务器上从头到尾的命令
scripts/lib_arms.py            ridge、内层 CV、paired bootstrap、canonical correlation
scripts/c0_prepare_annotations.py   下注释、归一化到 hg19 BED、写 c2 的 manifest
scripts/c1_build_population_features.py   P
scripts/c2_build_functional_features.py   F
scripts/c3_complementarity_gate.py        A/B/P/F/PF 与两个 delta
tests/test_complementarity.py  20 个测试
```

## 前置

Task B 的 panel 已冻结，每个 trait 有一份 `TQ_RESULTS.json`（validation-only 那次即可）。
注释文件需要外网取一次。withdrawal 名单问题不阻塞本步，但仍然阻塞 `tq_test` 的打开。

## 输出边界

只写 `GATE_RESULTS.json` / `PENALTY_TRACE.json` / `arms.csv` / `comparisons.csv`，
全部是汇总量。没有个体级行，包括假名化的。个体级数据全部留在服务器。
