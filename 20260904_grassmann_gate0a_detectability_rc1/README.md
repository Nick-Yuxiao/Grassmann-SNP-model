# Gate 0A Detectability Gate · rc1(planning proxy)

_状态：`PLANNING_PROXY`;日期：2026-09-04;不用真实数据,允许运行,但结果永不作 Gate 0A 正式证据。_

---

## 📋 干什么

在正式跑 Gate 0A 之前,先定死一件事:

> 在**固定**主 margin `ΔR²_genetic=0.005` 下,需要多少 simulation replicate `R` 才能可靠检测到它?决策程序在 `Δ=0` 的假阳性率是多少?

margin 先验定死;这个 gate 只定 `R`,**绝不**根据结果选 margin。

## 🎯 程序保真

完全复刻 Gate 0A 的决策:推断单元 = replicate;统计量 = 每个 replicate 的配对差 `Δ_r=R²_genetic(A)−R²_genetic(B)`;不确定性 = replicate 上的 percentile bootstrap CI;判据 = CI 下界 `>0`(单侧)。

## 🧮 shift/scale 恒等式(为什么又快又准)

把配对差建成 `Δ_r~Normal(Δ_true, σ²)`,`σ` 是配对 `R²_genetic` 差的 replicate 级 SD。由于 `L(Δ,σ)=Δ+σ·L₀`,判据 `L>0` 恰等于 `L₀>−Δ/σ`,所以每个 `R` 只需模拟一次单位下界 `L₀`,其余 `(Δ,σ)` 全闭式。

## 📊 冻结网格

- `Δ∈{0,0.002,0.005,0.010}`——假阳性、敏感性、**主 power**、强效应 power;
- `R∈{20,50,100,200}`;
- `σ∈{0.005,0.01,0.02,0.03,0.05}`——**规划轴**,真实值由真实 panel 小 pilot 定;
- Monte Carlo:`M=3000`、`B=1500`、seed `20260904`。

## 📈 本次结果(planning proxy)

见 `results/detectability_table.md`。达到 `Δ=0.005` 下 ≥0.90 power 所需最小 `R`:

| σ | 0.005 | 0.01 | 0.02 | 0.03 | 0.05 |
| --- | ---: | ---: | ---: | ---: | ---: |
| min R | 20 | 50 | 200 | >200 | >200 |

`FPR@0≈0.03`(单侧,符合预期)。**可行性警告:** 若配对 `R²_genetic` 的 replicate SD `σ≳0.03`,即使 200 个 replicate 也够不到 0.005 的 90% power——那时 Gate 0A 要么加 replicate,要么重新考虑 `N`/`h²`。这正是"真实 panel 的 σ pilot"必须作为运行前 gate 的原因。

## 🔍 复现

```bash
python -m unittest discover -s tests -p 'test_*.py'
python scripts/run_detectability.py
python scripts/validate_package.py
```

依赖 numpy。运行是纯 CPU、无数据。

## 🚫 不授权

- 不把这个 `R` 当作 Gate 0A 正式证据;
- 不在真实 panel 的 σ pilot 之前替真实 Gate 0A 定死 `R`;
- 不做任何真实基因型/表型工作;
- 不根据结果选 margin;不做 GPU/v7。
