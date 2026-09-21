# 注释覆盖结果（c0，2026-09-21，c20-5090）

panel: 20801 个位点，GRCh37。三个 UCSC 源全部 `200 OK`，走在线路径，无 lift-over。

## 正确性证据：15 个 ChromHMM 状态合计 = 0.9988

这些状态互斥且铺满基因组，所以按位点统计的覆盖率**必须**合计到 ≈1。实测 99.88%，
缺的 0.12% 落在 segmentation 的缝隙里。这一条同时验证三件事：基因组版本对（hg19，
不是 GRCh38）、BED 半开区间的边界处理对、没有重复计数。任何一条错了这个和都不会接近 1。

## 覆盖率

| track | coverage | 读法 |
|---|---|---|
| `13_heterochrom_lo` | 0.6517 | thin 过的位点大多落在基因间区，符合预期 |
| `11_weak_txn` | 0.1276 | |
| `12_repressed` | 0.0854 | |
| `10_txn_elongation` | 0.0433 | |
| `7_weak_enhancer` | 0.0272 | |
| `9_txn_transition` | 0.0139 | |
| `4_strong_enhancer` | 0.0089 | |
| `6_weak_enhancer` | 0.0088 | |
| `8_insulator` | 0.0075 | |
| `1_active_promoter` | 0.0070 | |
| `2_weak_promoter` | 0.0058 | |
| `5_strong_enhancer` | 0.0057 | |
| `14_repetitive_cnv` | 0.0032 | 勉强过 `--min-coverage 0.002` |
| `3_poised_promoter` | 0.0014 | |
| `15_repetitive_cnv` | 0.0013 | **低于阈值，已剔出 manifest** |
| `phastcons100way` | 0.0342 | 保守元件约占基因组 5%，稀疏采样后略低，对得上 |
| `tss_proximity` | 0.7645 | **见下，这个数不能按字面读** |

活跃调控状态合计约 7%（启动子 0.70% + 强增强子 1.46% + 弱增强子 3.60% + 绝缘子 0.75%）。

进入 manifest 的是 **15 条 track**。`--group-by chromosome` 下 F = 15 × 20 条染色体
= 300 列（chr18、chr21 不在 panel 里）。

## `tss_proximity` 的覆盖率是个假数，不要引用

`c0` 统计 coverage 的口径是 "score > 1e-6"。这条 track 的 score = `exp(-d/10kb)`，
反解得 d < 138 kb——等于在问"离最近的转录起点 138 kb 以内吗"，基因组上绝大部分位置
都满足，所以 0.7645 不代表"76% 的位点靠近 TSS"。

**这是 `c0` 的口径缺陷，不是数据问题**：那个覆盖率指标是为区间型 track 设计的，
套到连续型打分上就失去意义。track 本身可用——权重按指数衰减，实际有贡献的是
20–30 kb 以内那些，远处权重接近 0。

对外写报告时：区间型 track 可以引 coverage，`tss_proximity` 不要引，要引就引
score 的分位数。

## 结论

可以进 `c2`。基因组版本、坐标系、区间语义都已验证，没有需要回头的地方。
