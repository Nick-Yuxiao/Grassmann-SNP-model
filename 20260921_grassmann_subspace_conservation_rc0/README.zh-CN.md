# Cross-context effect-subspace conservation · rc0

_日期：2026-09-21｜状态：`PARKED_INCOMPLETE`｜Grassmann 主线已结案，本包为封存产物_

---

## ⚠️ 这个包是什么，不是什么

这是一条**未完成**的路线的封存代码。Grassmann 主线在 `89 → 86 → 8 → 4 → 0` 的真实
interaction 复现链之后已结案；本包在结案前写到一半，保留下来是因为里面的几何层和
估计量经过测试、可复用，**不是**因为这条路线还活着。

**本包不含任何真实数据结论。** 全部输出来自合成世界，只用于验证估计量本身。

## 📦 已完成（测试通过）

| 组件 | 状态 |
| --- | --- |
| `src/subspace_conservation/geometry.py` | 完成，11 项测试 |
| `src/subspace_conservation/nulls.py` | 完成 |
| `src/subspace_conservation/estimators.py` | 完成 |
| `src/subspace_conservation/arms.py` | 完成，表示臂与 transport 基线 |
| `src/subspace_conservation/simulate.py` | 完成，四个合成世界 |
| `scripts/l1_implementation_gate.py` | 完成，`PASS` |
| `scripts/l0_calibration.py` | 完成，4/4 世界正确判别 |
| **`scripts/l2_phase_proxy.py`** | **未写** |
| **`scripts/l0_real.py`（真实数据 runner）** | **未写** |
| **`PROTOCOL.zh-CN.md`** | **未写** |

```bash
python3 -m unittest discover -s tests -p 'test_*.py'      # 34 tests
python3 scripts/l1_implementation_gate.py --trials 200
python3 scripts/l0_calibration.py --pairs 150 --skip-detectability
```
仅需 numpy / scipy；纯 CPU，无数据，无 GPU。

## 🧮 两个值得留下的技术结论

**1. nuisance 群是 `GL(2)` 不是 `O(2)`。** 若 `E = B Λᵀ`，换 ancestry 即换 `Λ`，
在 `K=2`、`Λ` 可逆时这正是 `E → E M`，`M ∈ GL(2)`，而 `span(E M) = span(E)`。
关键推论：右乘正交矩阵下 `E Eᵀ` **精确不变**且比 projector 更富信息，所以一个
"只做旋转"的正对照会把分数送给 spectrum 臂而不是 Grassmann 臂。
`l1_implementation_gate.py` 实测 `O(2)` 偏差 `8.9e-16`、`GL(2)` 偏差 `0.61`。

**2. per-pair `GL(2)` 对齐残差恒等于 chordal subspace distance。**
`min_M ||Q_B − Q_A M||_F` 就是子空间距离（`test_geometry.py` 断言到 1e-10）。
即 Grassmann = "已经做过 per-pair GL(2) 自适应的 aligned 表示，且自适应成本为零"；
任何在测试 pair 自己的 target 数据上拟合 `M` 的所谓 baseline 都不是 baseline，
是 Grassmann 臂加泄漏。

## 🔍 估计量：四路判决

原始设计只区分"守恒 / 不守恒"，标定时被合成世界 `N` 证伪——**什么都没守恒时，
subspace 度量比 coordinate 度量更早饱和，于是仍然呈现 `I_sub < I_coord`，被误读成守恒**。
修法是引入 permutation（错配 pair）作为"完全无关"上界，把两个度量都归一到同一
绝对尺度 `C = 1 − (D_ctx − D_null)/(D_perm − D_null)`，并拆成四个互斥判决：

| 合成世界 | C_sub | C_coord | 判决 |
| --- | ---: | ---: | --- |
| `S` 子空间守恒 | +1.09 | −0.11 | `PROCEED_TO_LEVEL_1` |
| `C` 坐标也守恒 | +0.95 | +0.98 | `NO_CONTEXT_DISTORTION_ROUTE_UNMOTIVATED` |
| `N` 什么都不守恒 | +0.04 | −0.03 | `NO_SUBSPACE_CONSERVATION_STOP` |
| `D` 装饰性守恒 | +0.99 | −0.05 | `CONSERVED_BUT_DECORATIVE_STOP` |

`S` 与 `D` 的**单 SNP 边际完全相同**（每列的角位移和模长都匹配），只在两列是否
一起移动上不同。估计量能把它们分开，是这条路线唯一真正被验证过的东西。

## 🚫 不授权

- 不把任何合成数字当作关于遗传学的证据；
- 不在没有真实数据 Level 0 的情况下重启这条路线；
- 不把 `l1_implementation_gate.py` 的 `PASS` 当作科学结论——它只验证实现的不变性。
