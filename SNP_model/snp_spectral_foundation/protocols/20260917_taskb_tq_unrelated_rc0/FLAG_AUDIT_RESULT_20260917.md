# Field 22021 编码审计结果

_2026-09-17 · `c20-5090` · 来源 `ukb2.sas7bdat` · 状态 `FILTER_FROZEN`_

---

## ✅ 三项检查全部通过

| 检查 | 判据 | 实测 | 结论 |
| --- | --- | --- | --- |
| 取值集合 | 只含 `{0, 1, 10, -1, <blank>}` | 完全一致，`undocumented_values_present` 为空 | ✅ |
| 保留比例 | `keep_share` 落在 0.60–0.75 | **0.675763** | ✅ |
| 参与者覆盖 | `coverage_of_cross_reference` 接近 1 | **0.999239** | ✅ |

## 📊 观测直方图

`n_22021_0_0`，502,238 行，重复标识行 0：

| 取值 | 人数 | 占比 | 含义 | 处置 |
| --- | --- | --- | --- | --- |
| `0` | 339,394 | 0.6758 | 未发现亲缘关系 | **保留** |
| `1` | 147,447 | 0.2936 | 至少一名亲属 | 剔除 |
| `<blank>` | 14,232 | 0.0283 | 无记录 | 剔除（未知） |
| `-1` | 977 | 0.0019 | 未参与 kinship 推断 | 剔除（未知） |
| `10` | 188 | 0.0004 | 十名及以上三度亲属 | 剔除 |

合计：保留 339,394；已知有亲缘 147,635；亲缘未知 15,209。

QC flags：`n_22027_0_0`（杂合度/缺失率离群）968 人标记为 `1`；`n_22019_0_0`（性染色体非整倍体）651 人标记为 `1`。两者均在 `b3` 中以 `--qc-flag` 剔除。

## 🔬 独立佐证

已知有亲缘者合计 **147,635 人（29.4%）**，与 UKB 官方公布的「约 14.7 万参与者存在三度或更近亲属」一致。该吻合是我们读到的确为 field 22021 而非同名字段的旁证。

## 📐 「未知一律剔除」的实际代价

`flag` 文件 502,238 人，covariates（按 chr1 fam 对齐）487,409 人，交集 487,038，仅 flag 有的 15,200 人。

15,209 个「亲缘未知」中约 15,200 人**本就未做基因分型**，会被 `.fam` 取交集自然剔除。因此该规则在 genotyped 人群中的净损失只有数百人量级，而不是 15,209。这条规则的保守性几乎不花成本。

## 🔒 据此冻结的过滤条件

```text
keep  ⟺  n_22021_0_0 == 0
         且 n_22027_0_0 ∈ {空, 0}
         且 n_22019_0_0 ∈ {空, 0}
         且 存在于 genotype .fam
         且 不在撤回同意名单中
```

写入 `CONFIG.json`：

- `relatedness_axis.source = "ukb_field_22021_unrelated_subset"`
- `relatedness_axis.field_22021_kept_values = ["0"]`
- `relatedness_axis.flag_source_file = "ukb2.sas7bdat"`
- `relatedness_axis.observed_value_counts` = 上表
- `relatedness_axis.claim_scope = "new unrelated individual, not new family"`

## ⚠️ 仍未核对的一项

上表的取值含义取自本包 `FIELD_22021_CODING.zh-CN.md` 中待核对的编码表。观测分布与之完全自洽（五个取值、比例吻合官方统计），但**与 UKB Showcase data-coding 页面的逐条核对尚未由人完成**。在正式报告中引用这些含义前需补上该核对。

---

## 🚨 撤回同意名单缺失（2026-09-17 确认）

冻结的过滤条件包含「不在撤回同意名单中」这一条，但**该条件目前无法验证**。

在 `/data3` 与 `/home/tyuxiao` 的 6 层深度内搜索 `*44430*`、`w[0-9]*.csv`、`*withdraw*`、`*exclus*`、`*.ukbkey`，**全部无匹配**。`/data3/ukb_all/logs/` 只有迁移脚本与 rsync 日志；`/data3/ukb_all` 下的 csv 全是前一位研究者的分析输出。因此 `b3` 运行时 `explicit_exclusions_supplied = 0`。

### 当前的部分覆盖及其边界

`b3` 剔除了 111 个负号标识（`dropped_negative_identifier`）。负号 ID 是 UKB 在刷新数据中标记撤回参与者的惯例写法，因此这是真实的部分覆盖。

但存在时间差：genotype 来自 **2020 年 1 月**的 bgen 转换（见 `ukb_imp_chr1_v3.log` 的时间戳）。**2020 年之后撤回的参与者不会体现在那批负号 ID 中。** 这 111 人因此不能当作完整的撤回处理。

### 建议的处置顺序

不阻断 validation 那一跑——该跑不产生正式判定，只确认实现与数值稳定性。

**但应在打开 `tq_test` 之前解决。** 理由既是合规也是工程：若撤回名单后续到位并改变样本集合，`b3` 必须重跑，split 哈希随之变化，已经开过的 `tq_test` 就落在一个被取代的划分上。test 是一次性的，不应花在可能重做的 split 上。

### `CONFIG.json` 中必须如实记录

```json
"withdrawal_list_supplied": false,
"withdrawal_partial_coverage": "111 negative-sign identifiers dropped by b3",
"withdrawal_coverage_gap": "genotype converted 2020-01; post-2020 withdrawals are not reflected in those identifiers",
"withdrawal_resolution_required_before": "tq_test opening"
```

### 待向数据管理方确认

1. 是否存在最新的参与者撤回名单（通常为 `w44430_<日期>.csv`）
2. genotype 自 2020 年 1 月转换以来是否有过刷新
