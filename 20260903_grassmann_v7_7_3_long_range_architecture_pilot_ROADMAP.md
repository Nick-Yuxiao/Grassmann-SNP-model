# v7.7.3 自我规划说明（架构 pilot 合约，无 GPU）

## 到目前为止的链条

- v7.6.4：注册结论 `A1R_LD_REGIME_DEPENDENT`（不可变），没有任何 regime 展示 +0.010 nats 的 Grassmann 增益。
- v7.7.0：签署 long-range Gate 合约（合成远端 XOR、三级 Gate、2×2 阶乘、主对比与 GO 规则）。
- v7.7.1：实现 CPU-only 生成器 / 静态不变量 / 解析探针。
- v7.7.2：执行 CPU 任务有效性 → **`LONG_RANGE_TASK_VALIDITY_PASS`**。
  LR0（local ≈ marginal，90% CI ⊂ ±0.010）、oracle 成功、LR1 常规全局 MLP 明显成功
  （95% LCB≈0.694）、target-shuffled 仍等价。**结论：任务能区分"局部不足"与"全局可解"，
  但这不是 Grassmann 胜利。**

## 本次交付：v7.7.3（唯一被授权的下一步 = 起草架构 pilot，不授权 GPU）

按既有 patch 约定生成，全部本地验证通过（`validate` PASS，8/8 单测通过，builder 端到端跑通）。
它是一个**合约草案**：不训练 Grassmann、不碰 GPU、不给任何 arm 出结论。冻结的内容：

1. **完整 2×2 阶乘**：Grassmann(缺/有) × 共享全局 router(缺/有)，四个 decision-eligible 格
   （LR_A00 / LR_A10 / LR_A01 / LR_A11）+ oracle / shuffled 两个控制。
2. **共享全局 router 的关键澄清**：与 v7.7.2 任务有效性正控制不同，架构 router **不被告知源位置**，
   感受野 = 全 8192-token 序列、输入 = 全 token 序列；两个 router-present 臂（A01/A11）的实现、
   感受野、输入、优化调度、训练数据**完全相同**。这一点是主对比可解释的前提。
3. **主估计量**：`delta_primary = NLL(LR_A01) − NLL(LR_A11)`（router 条件下的增量 Grassmann 收益），
   实践边界 0.010 nats/target，未来 GO 需配对单侧 95% LCB > 0.010 且 validity/shuffled/fairness 全过。
   Grassmann 主效应与 Grassmann×router 交互为次要，不能覆盖主对比。
4. **公平性**：matched-parameter / matched-compute 标签在 non-identity 审计 + realized-compute 审计
   通过前一律禁用。
5. **复制与盲法方差 pilot**：重复单元 = 真值 seed；正式 seed 数**此处不定**，留给下一步的
   **CPU-only 盲法方差 pilot**（缩比 proxy、隐藏 arm 均值、只释放离散度与冻结后的 N，
   pilot 数据禁止进入正式分析）。

**终态**：`LONG_RANGE_ARCHITECTURE_PILOT_CONTRACT_SIGNED_NO_GPU`。
**唯一转移**：`IMPLEMENT_V7_7_4_CPU_BLINDED_VARIANCE_PILOT_NO_GPU`。

## 之后的路（规划，非授权）

- **v7.7.4**：实现 CPU 盲法方差 pilot（缩比 proxy 跑四臂，估 `delta_primary` 跨 seed 的 SD，
  按功效冻结正式 N；不释放 arm 均值、不出架构结论、不碰 GPU）。
- **v7.7.5**：把 non-identity + realized-compute 公平性审计实现并通过（解锁 matched-* 标签）。
- **v7.8.0（需全新、单独的 GPU 授权）**：在冻结的 N 与公平审计通过后，运行 GPU 2×2 阶乘，
  按主对比 GO 规则判定；若比较足够精确而无实践收益 → **关闭 v7 Grassmann-primary 路线**，
  而不是再找新的 outcome-driven 任务。

## 待你拍板

1. 是否继续实现 **v7.7.4 CPU 盲法方差 pilot**（仍无 GPU）？
2. 缩比 proxy 的规模（序列长度是否降采样、router 参数预算、samples/seed）——会写死进 v7.7.4。
3. 功效目标：以多大把握在 N 个 seed 上把 0.010 的边界与 95% LCB 区分开（决定冻结的 N）。
