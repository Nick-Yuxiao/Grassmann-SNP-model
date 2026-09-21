# G2PT 基线规格：证据账本与冻结决定

**日期：** 2026-09-21　**状态：** rc0
**配套文件：** `../configs/g2pt/CONFIG_FROZEN.json`、`../experiments/B10_B11_spec_discrepancy_matrix.yaml`、`../PAPER_SPEC.zh-CN.md`

本文件回答一个问题：**当"论文怎么说"和"代码怎么写"对不上时，我们的基线以哪个为准？**
答案是：**不选，两个都建，并在报告增量时写明是相对哪一个。**

---

## 1. 证据账本

```text
                        G2PT evidence ledger
                                │
        ┌───────────────────────┼───────────────────────┐
        ▼                       ▼                       ▼
  Paper-confirmed         Code-confirmed            Unresolved
  （v2 原文直读）          （fc13abd 直读）          （两者冲突或都没说）
                                │
  d = 64                  SwiGLU FFN               init: default vs xavier   (B10)
  heads = 4               默认 nn.Embedding         FFN:  swiglu vs gelu      (B11)
  translation head = 1    只对 q 做 Pre-LN          batch size
  60/20/20                无任何随机种子 (B2)        dropout
  MSE + AdamW             epoch 网格未固化           epoch / early-stop 设置
  GO 2023-07-27                                     5 个传播步骤 ↔ CLI 开关映射
  grid: p 阈值 × epoch
                                │
        ┌───────────────────────┴───────────────────────┐
        ▼                                               ▼
  Secondary（二手，待原文确认）
  lr = 1e-5 (grid 1e-3/1e-4/1e-5)
  wd = 1e-3 (grid 0/1e-5/1e-3)
  cv = 5-fold
                                │
                                ▼
                         Frozen baseline
                                │
                ┌───────────────┴───────────────┐
                ▼                               ▼
        G2PT-code-faithful              G2PT-paper-faithful
        (arm A: default + SwiGLU)       (arm D: xavier + GeLU)
                └───────────────┬───────────────┘
                                ▼
                           your method
```

---

## 2. 三条冻结规则

**R1｜每个字段必须带证据等级。** `CONFIG_FROZEN.json` 里没有裸数值，只有 `{value, status, source}`。四个等级：`CONFIRMED_PRIMARY` / `SECONDARY_PENDING_PRIMARY` / `CODE_OBSERVED` / `UNRESOLVED`。训练代码消费任何字段时，必须把该字段的 `status` 一并写进 run record——否则三个月后没人说得清某个数字是读来的还是猜的。

**R2｜`UNRESOLVED` 字段不许有默认值。** 它必须作为析因实验的一个轴出现。悄悄沿用代码行为，等于替论文做了一个它没做的选择。

**R3｜增量必须指明基线。** 任何"相对 G2PT 提升 X"的表述，写明是 `vs G2PT-code-faithful` 还是 `vs G2PT-paper-faithful`。若最终确认论文与代码确实不一致，**两个都报**。

---

## 3. B10 与 B11 分开的理由

上一版把两者合在一起是错的，它们的可处理性完全不同：

| | B10（初始化） | B11（FFN） |
|---|---|---|
| 性质 | 优化设置 | **函数类** |
| 能否事后补救 | 能——敏感性分析 | **不能**，SwiGLU 与 GeLU FFN 是两个不同的模型 |
| 对"增量是否真实"的影响 | 改变基线的强弱 | 改变基线**是什么** |
| 优先级 | 次 | **主** |

B11 更硬。但 B10 也不是纯训练细节——见下节的实测。

---

## 4. B10 的机制：一条撤回，一条换上

### 撤回

早先写过"默认初始化下 attention logits 一开始就接近 one-hot"。**这是错的，已被实测证伪。**

step-0 探针（`../experiments/probe_init_step0.py`，配置 d=64/h=4）：

| Lk | init | logit std | entropy / ln(Lk) | max prob |
|----|------|-----------|------------------|----------|
| 8 | default | 0.3299 | **0.978** | 0.1916 |
| 8 | xavier | 0.0036 | **1.000** | 0.1256 |

两种初始化在 step 0 **都近似均匀**。原因也不是"Pre-LN 抹掉了尺度"——恰恰相反，见下——而是 `nn.Linear` 的默认初始化 `kaiming_uniform_(a=√5)` 本身足够小，即使 `k_std=1.0`，logit std 也只有 0.33，离 one-hot 很远。

### 换上（这条是新的，而且更硬）

**结构发现：HiGT block 是 Pre-LN，但只对 query 做。**

```python
q_norm = self.norm_attn(q)                                   # 只有 q
attn_output = self.drop(self._xattn(q_norm, k, v, mask, ...))  # k, v 原样
x = q + attn_output
```

`k` 和 `v` 未经归一化直接进入投影。所以嵌入的初始尺度**确实**原样乘进 attention logits——logit std 比值 91×，与嵌入尺度比 (1.0 / 0.011) 完全吻合。

真正的后果不在 logits 的绝对大小，而在于：**`xavier_uniform_` 的 bound 依赖词表大小，于是它在三类节点之间造出一个尺度分层，而默认初始化让三类全部等于 1.0。**

| 节点 | 词表 V | xavier std |
|------|-------|-----------|
| SNP | 16,532 | 0.0110 |
| gene | 253 | 0.0794（SNP 的 7.2 倍）|
| system | 20 | 0.1543（SNP 的 **14 倍**）|

由于残差是 `x = q + attn_output` 且 `v` 未归一化，这个分层直接改变每个传播步骤的混合比：

| 传播步骤 | default | xavier | 倍数 |
|---------|---------|--------|------|
| SNP → gene（前向）| 0.1233 | 0.0163 | **弱 7.6×** |
| gene → system（前向）| 0.0884 | 0.0426 | 弱 2.1× |
| system → gene（**反向**）| 0.1413 | 0.2633 | **强 1.9×** |

**xavier 不是均匀地缩小一切，而是重新配平了前向与反向传播的相对强度**：前向被压低、反向被抬高，净差约 14 倍。

这一点为什么要紧：论文的生物学叙事核心之一，正是"系统状态反向调制其子系统与基因"（v2 Methods p15 的 reverse propagation）。**这条通路的初始强度，直接由一个论文与代码不一致的初始化选择支配。** 所以 B10 有架构级后果——但它仍可用敏感性分析处理，这一点与 B11 不同。

---

## 5. 执行顺序

1. ✅ 落盘：本文件 + `CONFIG_FROZEN.json` + `B10_B11_spec_discrepancy_matrix.yaml` + `probe_init_step0.py`
2. ✅ B6 → `PARTIALLY_CLOSED`
3. ✅ B10（init）与 B11（FFN）拆分建档
4. ⬜ 在复现分支实现 `--init {default,xavier}`、`--ffn {swiglu,gelu}`、`--seed`
5. ✅ step-0 诊断（阶段 1，已完成，结论见 §4）
6. ⬜ 2×2 短程优化诊断（阶段 2，A/B/C/D × 3 seeds × 5 个检查点）
7. ⬜ 拿到 v3 后，核对优先级：**initialization 原话 → FFN 原话 → grid-search 表 → batch/dropout/epoch**
8. ⬜ 据此决定哪一臂进入正式 Stage B 基线（可能是两臂）

---

## 6. 一句话

B10/B11 不是复现工作的边角料，是**"增量是否真实"的校准层**。现在花两天把它钉死，比 Stage B 跑完之后发现基线规格漂了便宜得多。
