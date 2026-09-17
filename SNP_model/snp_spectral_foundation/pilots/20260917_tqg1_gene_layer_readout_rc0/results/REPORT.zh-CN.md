# TQ-G1 rc0 结果报告：`TASK-INELIGIBLE`

_运行：2026-09-17｜评价角色 `task_gate`（50 人）｜`bridge_test`（230 人）仍密封_

## 判定

**本轮没有回答预注册的主问题,因为前提不成立。**

不能宣称 "annotation-free gene 池化劣于 additive dosage"。正确的判定是
`TASK-INELIGIBLE`：在本设置下 cis 遗传信号不可检测,主对比因此没有解释力。

## 数据

| Arm | 输入 | macro R² |
| --- | --- | --- |
| A | sex + population + 10 PC | **+0.0263** |
| B | A + 256 additive cis dosage | **+0.0119** |
| C_gene | A + annotation-free gene 池化读出 | **−0.0196** |
| C_full | A + 256 × 16 对齐 full `H` | **+0.0036** |

| 对比 | ΔR² | 95% CI | 判定 |
| --- | --- | --- | --- |
| `C_gene − B`（primary） | −0.0315 | [−0.0540, −0.0124] | 不支持 |
| `C_full − B` | −0.0083 | [−0.0199, +0.0050] | 跨 0 |
| `C_gene − C_full` | −0.0232 | [−0.0459, −0.0052] | 池化劣于 full `H` |
| **`B − A`** | **−0.0144** | **[−0.0386, +0.0122]** | **跨 0 → 任务不可检测** |

301 个合格基因,按冻结的 hash 规则取 200 个分析。

## 为什么是 TASK-INELIGIBLE 而不是负结果

`B − A` 跨 0 意味着 additive dosage 相对协变量没有增量。四个臂的 macro R² 全部
落在 0 附近（−0.02 至 +0.03）。在没有可检测遗传信号的回归上比较表示,测到的是
**哪个臂过拟合更少**,不是哪个表示携带更多信息：dev_train 只有 120 人,B 加 256
维、C_full 加 4096 维、C_gene 加 160 维,ridge 收缩后净贡献为负是预期行为。

因此 `C_gene − B = −0.0315` 的统计显著性是真实的算术事实,但它不支持任何关于
gene 单位的科学结论。

## 协议缺陷（必须记录,不得事后删除）

rc0 把 `B − A` 写成"次要对比,不作为 primary 的前置 Gate"。这是设计错误。

Stage 2B rc1 曾把同一个量作为 Task Gate,并在其不过时判 `TASK-INELIGIBLE`、
拒绝打开 bridge test（当时 `B−A = +0.0105`, CI `[−0.038, +0.060]`）。rc0 取消这
道 Gate 的理由是"gene 作为复制单位提供足够功效"——该推理混淆了两件事：

- 检测**两臂之差**的功效：200 个基因确实提供了远高于 8 个 trait 的功效；
- **存在可测信号**这一前提：与基因数无关,由个体数与基因选择方式决定。

rc0 同时移除了 rc1 的信号富集（基因按 hash 随机取、cis SNP 按几何抽稀,均不看
phenotype）。该选择在方法学上仍然正确——phenotype 驱动的选择会抬高 B 并让
`C−B` 难以解释——但它与 `n_eval = 50` 结合后使任务整体不可检测。

## 一个不能被这个负结果解释掉的正面观察

共享 encoder 的 masked-genotype validation：

```
model_accuracy 0.7209 / marginal(HWE) 0.6494 / contextual_lift 0.0715
masked_cross_entropy 0.6292
```

encoder 相对 AF/HWE 基线有 7.15 个点的 lift,所以本轮结果**不能**归因于
"encoder 未学到任何东西"。但该 lift 只针对 E0 层级中最弱的 M0b 基线,对是否
胜过 M1a 局部 multinomial 与 M1b haplotype-HMM 不构成证据。

## 本轮停止了什么

只停止一件事：**在 `n_eval = 50`、无信号富集的随机 chr18 基因样本上,用本设置
检验 gene 单位增量**。

不停止：gene layer 这一想法本身、annotation-free routing、SNP foundation
modelling、E0 masked-genotype gate,也不构成 `bridge_test` 的任何结论。

## 资产状态

- `bridge_test` 230 人：**仍然密封,未被读取**。
- `task_gate` 50 人：本轮打开一次；该组个体在 Stage 2B rc1 中已被打开过,本来
  就不是 virgin test set。

## 已知缺陷

见 `../IMPLEMENTATION_INCIDENT_01.md`：`GENE_TABLE.tsv` 的浮点精度不足以让独立
validator 逐位重现两个 Spearman（偏差 `1.3e-4` / `3.0e-5`）。不影响 primary
estimand 与任何结论,修复进入 rc1。
