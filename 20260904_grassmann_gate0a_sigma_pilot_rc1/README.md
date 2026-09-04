# Gate 0A σ-pilot · rc1(variance calibration)

_状态：`VARIANCE_CALIBRATION_ONLY`;日期：2026-09-04;只在**已绑定**的冻结契约上运行,只估方差,不判胜负。_

---

## 📋 干什么

估计 Gate 0A 决策程序里配对差 `Δ_r=R²_genetic(A_test)−R²_genetic(B_pca_z)` 的 replicate 级 SD `σ`,再据此定所需 replicate 数 `R`。

## 🔒 顺序(绝不反)

绑定 panel hash → 绑定 LD-block 版本 hash → 冻结/验证预处理(`BINDING.json → BOUND`)→ 导出 panel → 跑 pilot → 估 `σ` → 映射到 `R` → 冻结 `R_formal`。run 脚本在契约未 `BOUND` 时直接 `RUN_BLOCKED`,所以在未冻结 SNP 集/block 版本/split 上测到的 `σ` 永远进不了 `R`。

## 🎯 范围:只估方差

只估 `SD(ΔR²_genetic)`,**不**判谁赢。结果只能定 **replicate 数与算力可行性**——不能改 margin、不能决定 DGP 是否保留、不能改 `k`/arm/success threshold、不能重定哪个 regime 是 primary(这些在 rc2 已冻)。

## 🧮 为什么不需要 h² 网格

主指标 `R²_genetic` 把下游拟合到**无噪声**遗传值 `g`,评 `R²(g,ĝ)`;`g` 不含表型噪声,故 `ΔR²_genetic` 的 replicate SD **与 h² 无关**,pilot 不需 h² 网格(h² 只进 Gate 0A 报告用的次指标 `R²_pheno`)。

## 🧪 设计

- regime(代表性):`spectral_tail_adversarial`(必须,主门控)、`major_LD_aligned`(易参照)、`between_block_interaction`(高方差参照,双线性 head);
- per-block `k∈{8,16}`;`R_pilot=30`;外层 5-fold CV;ridge λ∈{0.01,0.1,1,10};3 个因果方向;
- 臂只在训练折拟合;MAF-z 用训练折统计量;主臂 KPCA-RBF(训练折 median 带宽),null `B_pca_z`;
- 加性 regime 读有界 block 集(≤40,记录);交互 regime 用涉及 block 对的匹配双线性 head。

## 📐 σ→R 与 R_formal

每个 cell 的 `σ̂` 用与 detectability gate 完全相同的判据(percentile-bootstrap CI 下界 `>0`)映射到 `{20,50,100,200,400,800}` 中达到 `Δ=0.005` 下 ≥0.90 power 的最小 `R`。然后

```
R_formal = 主相关 cell 上 required R 的最大值(绝不取平均)
```

若无候选 `R` 够用 → 加 `R`,绝不降 margin。

## 🧬 参考实现

`src/gate0a` 实现臂(`pca_z`/`kpca_rbf`)、逐 regime DGP(`g` 落在 block 谱尾/主方向/跨 block 乘积)、匹配双线性 head、以及嵌套 CV 下的 `R²_genetic`。它是 Gate 0A 将复用的参考实现,纯 numpy,已在合成 panel 上单测通过(`tests/test_gate0a.py`);真实 `σ` 来自服务器上对已绑定 panel 的运行。

## 🔍 复现(机器逻辑)

```bash
python -m unittest discover -s tests -p 'test_*.py'   # 合成数据自测
python scripts/validate_package.py
# 真实运行见 server_ops/SERVER_STEPS.sigma_pilot.md(需 BINDING==BOUND)
```

## 🚫 不授权

- 未 `BOUND` 不跑;长度不符的 panel 拒跑;
- 不判胜负、不改 margin/DGP/k/arm/threshold/primary;
- 不做 GPU/v7。
