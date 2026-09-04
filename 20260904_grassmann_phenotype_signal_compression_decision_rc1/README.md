# Grassmann 表型信号压缩判决实验 · 设计包 rc1

_状态：`DESIGN_ONLY_RUN_NOT_AUTHORIZED`；日期：2026-09-04；只冻结设计，不授权运行。_

---

## 📋 这个包在回答什么

在继续任何 Stage-1(谱压缩)、Stage-2(Grassmann mixing)、A1-R 或 phenotype gating 之前,先把整座架构赖以成立的**上游那一件事**测清楚:

> 对一个 LD block,学到的非线性重建表示的 top-`k`,是否比同 `k`、信息匹配的最佳线性基线,保留更多的表型信号?

若答案是否定的,Stage 1 + Stage 2 没有存在理由,GPU / A1-R 继续停。这就是四层设计跳过的"地基承重测试"。

**本实验不训练任何 Transformer、不用 GPU**,全部闭式 / CPU,目标是一台机器一周内跑完。

## 🎯 为什么 null 是 block PCA(而不是 PLS)

压缩目标是**无监督重建**。在重建(方差最大化)目标下,线性 encoder 的 top-`k` 子空间**就是** block PCA 子空间;因此接近线性的 encoder 在数学上退化成 block PCA,block PCA 就是测试臂必须打赢的精确线性上界。

监督式 null(block-PLS)被**排除**:它会让测试臂靠 encoder 从没见过的标签信息"赢",把"非线性 vs 线性"偷换成"用了标签 vs 没用"。null 必须只看到 encoder 看到的信息——重建 encoder 对应的就是 block PCA。

## 🧪 对照臂(同 `k`,输入/split/性状逐字节相同)

| 臂 | 表示 | 角色 |
| --- | --- | --- |
| `A_test`(主) | RBF kernel-PCA top-`k`(闭式) | 被测的非线性重建压缩器 |
| `A_test`(次) | 浅层 per-block autoencoder 瓶颈 `k` | 仅稳健性,不参与判决 |
| `B_pca` | block PCA top-`k` | **主线性 null** |
| `B_rand_bilinear` | 同维随机双线性/二次特征 | **恢复的消融**:几何 vs 一般非线性 |
| `B_rand_proj` | 高斯随机投影到 `k` 维 | 平凡下参考 |
| `B_ldprune` | LD-pruned SNP 子集,数目匹配 `k` | 朴素 SNP 选择参考 |
| `B_blockmean` | block 均值 / 简单 haplotype dosage | 地板参考 |

`B_rand_bilinear` 是上一版设计删掉的那条,和 `B_pca` 不重复:PCA 回答"非线性有没有超过线性最优",随机双线性回答"Grassmann 几何比一个同维二次特征多做了什么"。两个问题必须画在同一条曲线上。

## 📊 下游评价

固定为岭回归(连续)/逻辑回归(二分类),惩罚项由训练折内嵌套 CV 选,5×5 嵌套 CV。每臂每个 `k` 报 R²/AUC 随 `k` 的曲线,单性状与多性状聚合都报。

## 🧬 表型模拟器(把两个读数统一)

该 panel 是 1KGP 基因型,没有真实表型,所以表型信号读数全部在**模拟**性状上进行(已知遗传结构):LD 取自冻结的 v7 chr22 panel、多基因背景、已知因果位置、h² 网格 `{0.05,0.1,0.2,0.4}`,因果效应在 **within-block** 与 **between-block** 两个 regime 分别注入、分别报告。

R²/AUC-vs-`k` 曲线和 causal-retention-vs-压缩率曲线是**同一批**模拟运行的两个读数。`causal_signal_retention` = 下游线性模型从 `k` 维表示能恢复的注入因果方差比例(相对全基因型上限)。**between-block** 那条曲线是"LD-preserving / long-range 压缩"的唯一直接证据——也就是 Grassmann 的真正卖点。

## 🔪 现在就写死的判据(运行前冻结)

1. **Stage-1/2 kill:** 若 `A_test` 在 within-block 与 between-block **两个** regime、跨性状面板都打不过 `B_pca`(margin 的 CV 不确定区间需排除零),则 Stage 1 + Stage 2 无存在理由。
2. **几何特异性 kill:** 若 `A_test` 打赢 `B_pca` 却打不赢 `B_rand_bilinear`,则赢的是一个一般二次特征,Grassmann 特异 claim 失败——哪怕非线性确实有用。

margin 阈值与 CV 不确定规则在运行前冻结,看到臂间结果后不得调整。

## 🚧 运行前 gate(缺一不可)

1. 修好并验证此前指出的**中心化 bug**——每条重建/PCA 臂都依赖正确的 per-block 去均值,未中心化的 PCA null 不是合法上界;
2. 把数据契约绑定到冻结的 v7 chr22 panel manifest hash 与 block-version hash(见 `DATA_CONTRACT.json`);
3. 冻结 margin 与不确定规则;
4. 取得项目负责人的运行授权——然后才写运行代码。

## 🔍 设计包自检(不执行实验)

```bash
python scripts/validate_design_package.py
```

该脚本只校验设计包自身完整性与 manifest,不训练、不碰数据、不执行任何一条对照臂。

## 🚫 不授权事项

- 不执行本包的任何运行(设计包)
- 不训练任何 Transformer,不做 phenotype gating
- 不做 GPU 或 v7 / A1-R 工作
- 不读 HGDP holdout
- 不产生真实表型预测 claim
- 不产生超出上述 kill 逻辑的架构 GO/NO-GO
