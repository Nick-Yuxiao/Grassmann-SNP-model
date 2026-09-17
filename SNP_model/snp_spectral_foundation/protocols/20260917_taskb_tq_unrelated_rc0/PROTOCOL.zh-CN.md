# Task B 资格检验协议：unrelated-subset 上的 A vs B

_版本 `rc0` · 2026-09-17 · 状态 `STRUCTURE_FROZEN_PENDING_22021` · 不覆盖 E0 协议，也不是 Stage 2 bridge_

---

## 📋 这一步问什么

只问一个问题：

> 在 relatedness-disjoint 的 held-out 参与者上，**raw additive dosage 相对纯 covariate 模型，是否增加了样本外的 phenotype 信息**？

这是 bridge 的**前置资格检验**，不是 bridge 本身。通过只说明「该 trait 在该 panel、该 split 上存在可检测的 genotype 信号」，也就是后面问「学到的表示是否胜过 raw dosage」这个问题才有意义。

不回答：任何表示是否优于 raw dosage、生物特异性、因果性、Grassmann 是否有价值。

## 🧬 为什么 Task B 不被 imputed-only 阻断

E0 被 imputed 阻断，是因为它的**预测目标**是 genotype 本身——用 IMPUTE4 的输出当 target，而最强基线 `M1b` 又是同类 HMM，构成循环。

Task B 的预测目标是 **phenotype**，genotype 只作为**输入**。imputed dosage 作为输入是标准做法（全部 UKB GWAS/PRS 都这么做），不产生循环。因此：

| 资产 | E0 | Task B |
| --- | --- | --- |
| array 直接分型 | **必需** | 不必需 |
| phased haplotype | phase claim 必需 | 不必需 |
| 可靠 relatedness 轴 | **必需** | **必需** |
| phenotype | 禁止 | 必需 |

Task B 真正的前提只有三条：**可靠 relatedness 轴、phenotype 字段确认、A/B 独立资格检验**。

## 👥 Relatedness 轴

### 首选：UKB field 22021

保留 `22021 == 0` 的参与者，即 UKB 自身的 KING 推断未发现三度以内亲属者。每个保留个体因此是一个 singleton relatedness component，**个体级划分即 relatedness-disjoint**。

必须同时承认代价：

- `22021` 是**每人一个汇总值**，不是配对边；无法重建 connected components
- 结论范围从「新家系」收窄为「**新的无亲缘个体**」
- relatedness-disjoint 的强度取决于 UKB 自身 kinship 推断的质量
- 约损失 30% 样本，剩余约 34 万人，对本资格检验绰绰有余

### 备选：自建 KING

若 `22021` 不可得，用现有 chr1–22 hard-call 做 LD-pruned KING，生成 `.kin0` 后交给 E0 的 `r2_build_family_manifest.py` 冻结 family 轴。成本与精度折扣见 [`KING_FALLBACK.zh-CN.md`](KING_FALLBACK.zh-CN.md)。

### 明确禁止

**不得用随机个体划分开 bridge。** UKB 中约三成参与者有三度以内亲属；随机划分会让跨 split 的亲缘对把 genotype→phenotype 的家系共享成分泄漏进 test，测到的是泄漏而不是泛化。

## 📊 Trait 选择

优先连续、测量稳定、遗传信号强的表型，按此顺序：

| 优先级 | UKB field | trait |
| --- | --- | --- |
| 1 | `30020` | haemoglobin concentration |
| 2 | `30010` | red blood cell count |
| 3 | `30040` | mean corpuscular volume |
| 4 | `30780` | LDL direct |
| 5 | `30760` | HDL cholesterol |
| 6 | `30870` | triglycerides |

本轮最多评估 **3 个 trait**，在开 test 前写入 `CONFIG.json`。每个 trait 的 test 只开一次。

## ⚙️ 两个臂

| 臂 | 内容 |
| --- | --- |
| **A** | intercept、genetic sex、age、age²、array、assessment centre、PC 1–40 |
| **B** | A 的全部 + raw additive dosage |

B 臂中 dosage 通过一个**只在 train 上拟合**的线性预测子进入：先用 train 系数把 phenotype 与每个 SNP 都对 covariates 残差化，在 train 上估计每个 SNP 的边际效应，再按 p 值阈值加权求和成一个分数。这就是同 panel 的 raw/additive genotype 基线，与 E0 协议中 phenotype 路线的 baseline 定义一致。

阈值网格 `{1, 0.5, 0.1, 0.05, 0.01, 1e-3, 1e-4, 1e-5, 1e-6}` **只在 validation 上选一个**，选定后冻结。

## 📐 判定

主终点是 test 上的样本外增量决定系数

\[
\Delta R^2 = R^2_B - R^2_A,
\]

其中 \(R^2 = 1 - \sum (y-\hat y)^2 / \sum (y-\bar y_{train})^2\)，`train` 均值作参照以避免泄漏。

| 结果 | 判定 |
| --- | --- |
| `ΔR² ≥ 0.005` 且 bootstrap 95% CI 下界 `> 0` | `QUALIFIED` |
| CI 下界 `> 0` 但 `ΔR² < 0.005` | `STATISTICALLY_POSITIVE_PRACTICALLY_TIED` |
| 其余 | `NOT_QUALIFIED` |

`0.005` 是本轮预注册的最小有意义增量。2,000 次个体级 bootstrap；因为设计上全部为 unrelated singleton，个体即独立重采样单位。

**负对照**：把 train 拟合的分数在 test 内打乱后重算。真实信号会让该值明显为负；若它接近或超过 `ΔR²`，说明分数与结局错位，该次运行作废。

## 🔒 Firewall

- Test 目录默认关闭；`b5` 不加 `--allow-test` 只出 validation 结果
- 开过 test 会写 `TEST_OPENED.marker`，同目录再开会被拒绝
- AF、缺失填补、covariate 系数、SNP 效应、阈值选择全部只用 train/validation
- Trait 列只在 Task B 使用；**E0 不得读取本包产出的 phenotype 文件**
- 所有 summary 只含计数与哈希，不含参与者标识

## ⚠️ 已知限制

1. **Allele 方向**：本 cohort 的 `.bed` 由 BGEN 转换时未指定 `--bgen ref-first`，A1/A2 与官方 REF/ALT 的对应未验证。panel 因此把 dosage 定义为「`.bim` A1 等位基因计数」，且全流程不与任何外部 reference 比对。对本检验无害，但任何跨 cohort 断言前必须先做方向审计。
2. **Imputed 输入**：dosage 来自 imputation v3 hard-call。作为**输入**合法，但 panel 的有效独立位点数低于同名义数量的直接分型位点，`ΔR²` 不应被解释为可捕获遗传度的上界。
3. **Panel 是稀疏抽样**：按物理距离抽薄，不是全基因组 PRS。`ΔR²` 是资格信号的下界，不是最优预测性能。
4. **22021 的推断质量**：relatedness-disjoint 的强度继承自 UKB 自身的 kinship 推断，本包无法独立验证。

## ▶️ 通过之后

`QUALIFIED` 只解锁下一个问题：在同一 split、同一 panel、同一 covariate 集合下，**冻结的 genotype representation 是否胜过 raw dosage**。那需要单独的预注册、单独的 test firewall，且不得用本轮的 test 结果反选表示、层数、rank 或 checkpoint。
