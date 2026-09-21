# Field 22021 的编码语义与缺失值处理

_冻结过滤条件前必须核对的内容。字段存在 ≠ 过滤条件正确。_

---

## 📋 待核对的官方编码

UKB field `22021`（Genetic kinship to other participants）是**每个参与者一个汇总值**，不是配对边。待核对的取值与处置：

| 取值 | 含义 | 亲缘状态 | 处置 |
| --- | --- | --- | --- |
| `0` | 未发现亲缘关系 | **已知无亲缘** | **保留** |
| `1` | 至少识别出一名亲属 | 已知有亲缘 | 剔除 |
| `10` | 识别出十名及以上三度亲属 | 已知有亲缘 | 剔除 |
| `-1` | 该参与者被排除在 kinship 推断流程之外 | **未知** | **剔除** |
| 空值 | 无记录 | **未知** | **剔除** |
| 其它未记载取值 | — | **未知** | **剔除** |

上表必须与 UKB Showcase 上 field `22021` 的 data-coding 页面逐条核对后才算冻结。本包**不以记忆中的编码为准**：`b2b_audit_flag_coding.py` 打印实际观测到的取值直方图，冻结依据是那张表。

## ⚠️ 最关键的一条：`-1` 与空值必须剔除，不能当作无亲缘

`-1` 的含义是「未参与 kinship 推断」，也就是**这个人的亲缘状态从未被评估**。把它并入 `0` 等于把「不知道有没有亲属」当成「没有亲属」，而未知亲缘个体恰恰是本设计要排除的那一类泄漏源——一旦他们的亲属落在另一个 split，genotype→phenotype 的家系共享成分就会跨过 split 边界。

空值同理。任何未记载的取值也按未知处理。

代价是样本量：`-1` 与空值的比例直接从可用样本里扣除。`SPLIT_SUMMARY.json` 的 `relatedness_axis.dropped_by_value` 会逐值给出这笔账。

## 🔍 冻结流程

```bash
# 1. 抽列（b2）之后，先只做审计，不做任何过滤
python3 scripts/b2b_audit_flag_coding.py \
  --table "$OUT/taskb_extract/kinship_flag.tsv" --id-column n_eid \
  --column n_22021_0_0 --column n_22027_0_0 --column n_22019_0_0 \
  --cross-reference /data3/ukb_all/phenotype/derived/ukb_covariates_chr1fam_aligned.tsv \
  --cross-id-column eid \
  --out-dir "$OUT/taskb_flag_audit"

cat "$OUT/taskb_flag_audit/FLAG_AUDIT.md"
```

看三件事：

1. **观测到的取值集合**是否只有 `{0, 1, 10, -1, <blank>}`。出现未记载取值时先停下来查 Showcase，不要直接跑 `b3`。
2. **`keep_known_unrelated` 的比例**是否符合预期量级。UKB 中有三度以内亲属的参与者约占三成，因此 `0` 的占比应在 0.6–0.75 之间。明显偏离说明列选错了或编码不同。
3. **`participant_overlap`**：`ukb2` 与 covariates 是否覆盖同一套人。`coverage_of_cross_reference` 明显小于 1 时，`ukb2` 不能作为 flag 来源，应改用 `new.sas7bdat` 或 `ukb_alldata.sas7bdat`。

三条都确认后，才把取值集合写进 `CONFIG.json` 的 `relatedness_axis.field_22021_kept_values` 并运行 `b3`。

## 🧾 b3 中的对应实现

- `--keep-flag-value` 默认只有 `0`；要改必须显式传参，改动会记入 summary
- 标识不在 flag 文件中 → 计入 `dropped_no_kinship_flag`
- 取值不在保留集合中 → 计入 `dropped_related_or_unassessed`，并在 `dropped_by_value` 中逐值拆分
- summary 的 `relatedness_axis` 同时保存 `observed_value_counts`、`dropped_by_value`、`coding_reference` 与 `unknown_relatedness_is_dropped_not_kept: true`

## 🧪 对应测试

`tests/test_taskb.py` 的合成 cohort 覆盖 `0 / 1 / 10 / -1 / 空值` 五种取值，并断言：

- 只有 `0` 被保留
- `-1`、空值、`10` 的剔除计数逐值正确
- `b2b` 把「已知有亲缘」与「亲缘未知」分别计数，且未记载取值列表为空
