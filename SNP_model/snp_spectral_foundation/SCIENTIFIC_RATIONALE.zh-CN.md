# 粘贴方案、旧方案与新候选模型的科学关系

> ⚠️ **路线更新（2026-09-15）：** 本文此前建议将 foundation-model attribution audit 作为立即主线，该建议已被 0526 Stage 2 冻结协议取代。当前主线是不再重做 Stage 1，而是冻结一个 genotype representation，直接进行 `Full H / aligned low-rank / Grassmann → phenotype` 五臂验证。权威协议见 `protocols/20260915_0526_stage2_bridge_rc1/PROTOCOL.zh-CN.md`。

## 总裁决

我认同粘贴方案作为**当前主研究策略**：先审计公开 genotype foundation model 的能力由哪些 population-genetic structures 支撑，再决定新模型需要保存什么。它比直接把 Grassmann 扩成 backbone/memory 更符合已有证据，也更容易产生可否证、可复用的结果。

我不认同把审计方案直接解释为“已经证明应该训练下面这套新大模型”。本目录的模型是 audit 后的 Phase II 候选；在 audit 结果出来前，只能搭建接口、对照和 gate，不能预设 `Grassmann U` 有增量价值。

## 与旧方案的本质区别

旧方案（第 13 次决策）的核心链条是：

`genotype → context h_j → cross-fitted phenotype projector → effect e_j → small-set joint geometry → prediction / characterization`

粘贴方案的核心链条是：

`frozen public foundation model → MAF/LD/haplotype/ancestry 干预审计 → 可解释部分与 residual → 再问 phenotype/effect geometry`

| 维度 | 旧方案 | 粘贴方案 |
|---|---|---|
| 第一研究对象 | phenotype-conditioned effect representation | 已训练 foundation model 的 capability / representation |
| 第一科学问题 | joint-effect geometry 是否可估计、是否有用 | 模型能力由哪些遗传结构解释 |
| Grassmann 的位置 | small candidate set 的候选几何 | residual 出现后才进入的候选假说 |
| 对 phenotype 的依赖 | 早，需要严格 cross-fit | 前三类 audit 可 phenotype-free，最后才接 phenotype |
| 主要风险 | 小样本 outcome overfit、effect estimand 不清 | 干预不正交、把结构破坏误读为唯一贡献 |
| 资源使用 | 训练 projector/预测模型 | 先冻结公开 checkpoint，主要做干预和归因 |
| 对旧实验的复用 | 复用 effect-space 与公平表示比较 | 复用全部 perturbation/invariance/residual/context-range 工具 |

二者不是互斥方案。合理顺序是：**粘贴方案先做 representation audit；旧方案作为 residual 通过后的 phenotype bridge。**

## 为什么不能把归因写成简单的“贡献百分比”

MAF、LD、haplotype 与 ancestry 不是正交因素。比如逐 SNP 置换保留边际频率，但会同时破坏 local LD、long-range LD、haplotype 与部分 ancestry signal。一次性能下降不能唯一归给 LD。应使用多干预、负对照和顺序敏感性分析，并把结论写为“该干预破坏的一组结构对性能是必要的”，除非额外设计识别出独立成分。

同样，embedding residual `r = h - f(MAF,LD,H,A)` 依赖表示坐标和回归器容量。需要：

- 在训练折拟合解释器、在未见个体计算 residual；
- 用线性与非线性解释器形成容量曲线；
- 同时报告 residual norm、跨 seed/cross-cohort 稳定性和 downstream incremental utility；
- 对旋转不定的 hidden states 使用 CKA/CCA/Procrustes 或可预测性，而非逐坐标相减。

## 对用户给出流水线的修订

建议把它理解为四层，而不是一条不可回退的流水线：

1. **SNP reader**：raw genotype、AF/MAF、位置输入，masked-genotype local pretraining；必须胜过 train-only marginal baseline。
2. **可审计 block memory**：Grassmann projector、eigenvalues、MAF summary 与独立 signed residual sidecar；必须与 aligned/oriented control 同预算比较。
3. **可扩展 global context**：latent slots 在 block 数上线性扩展；必须通过 context-range destruction 才能声称利用 long-range signal。
4. **trait effect 与预测**：trait identity conditioning、effect projection、spectral/residual fusion；用于解释时必须 cross-fit。

本实现特意没有直接传递任意符号的 eigenvector 坐标。Grassmann 分支使用 projector probe `aᵀUUᵀa`，满足 `U → UQ` 不变性；eigenvalues 独立保存。由于旧实验已证明纯 quotient 可能丢掉任务所需的 signed/aligned coefficients，另设位置对齐 residual sidecar，并提供 `aligned_control`。

## 决策顺序

1. 立即主线：SNPBag/GenoBERT（以及相位对照）冻结审计。
2. 同步工程：维护本目录的小规模候选实现与真实数据 adapter，但不全规模训练。
3. audit 得到清晰设计原则后：冻结模型配置和消融合同，在未使用染色体/外部 cohort 上训练。
4. 只有 hybrid 稳定胜过 aligned control 且增量来自 spectral branch，才保留 Grassmann；否则删除该分支，保留 local encoder + residual + scalable global mixer。
