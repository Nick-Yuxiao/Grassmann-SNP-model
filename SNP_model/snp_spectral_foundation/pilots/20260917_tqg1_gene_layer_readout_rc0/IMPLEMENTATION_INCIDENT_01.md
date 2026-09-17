# 事故 01：GENE_TABLE 精度不足导致 validator 的 Spearman 复算失败

_发现：2026-09-17,rc0 正式运行之后｜影响：报告精度,不影响任何结论_

## 现象

`validate_tqg1.py` 的 23 条检查中有 2 条 FAIL：

```
rehearsal spearman[model]  -> reported 0.5724818120453011  recomputed 0.5726164979951035
rehearsal spearman[linear] -> reported 0.8926978174454361  recomputed 0.8926674141904176
```

差值分别为 `1.3e-4` 与 `3.0e-5`。

## 原因

`run_tqg1.py` 写 `GENE_TABLE.tsv` 时把 `score_model_C_gene`、`score_linear_B` 与
`target_T_g` 格式化为 `:.6f`,而 `validate_tqg1.py` 的容差是 `1e-9`。

`s_g` 是扰动前后预测差的标准差；在本次运行的强收缩 ridge 下其量级足够小,6 位
小数会丢掉有效数字并制造并列值。Spearman 依赖秩次,少量并列即可让复算结果在第
4 位小数上偏离 runner 内部的全精度计算。

确认方式：`GENE_TABLE.tsv` 中该列的去重计数小于基因数即为并列。

## 影响范围

- **不影响** primary estimand。`macro_r2[A/B/C_gene/C_full]`、`macro_delta_r2`
  与 pass rule 四项复算全部 PASS,因为 R² 的量级在 6 位小数下无有效数字损失。
- **不影响** 任何科学结论。两个 Spearman 的差异在第 4 位小数,不改变
  `linear > model` 的方向,也不改变 `rehearsal_model_beats_linear = false`。
- **只影响** 独立复算的可验证性：当前 `GENE_TABLE.tsv` 不足以精确重建 runner 的
  Spearman。

## 处置

rc0 已冻结,**不修改其代码、配置或产物**。`FREEZE.sha256` 保持原值。

修复进入 rc1：

1. `GENE_TABLE.tsv` 的全部浮点列改用 `repr()` 级精度写出,让独立复算可以逐位重现。
2. validator 对 R² 与 Spearman 分别设置容差,不再用单一的 `1e-9`。
3. validator 增加一条检查：浮点列的去重计数必须等于基因数,否则直接报精度不足。

## 教训

冻结产物的数值精度必须按**下游最敏感的统计量**来定,而不是按肉眼可读性。
秩相关对并列极其敏感,6 位小数对 R² 够用,对 `s_g` 不够。
