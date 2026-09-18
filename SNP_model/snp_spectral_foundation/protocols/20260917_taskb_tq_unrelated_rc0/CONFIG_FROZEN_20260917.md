# CONFIG 冻结凭据

_2026-09-17 · 状态 `FROZEN_PENDING_TEST` · `ready_to_open_test: true`_

---

## 🔑 冻结哈希

```text
CONFIG.json sha256 = 964e64715c9e32a2ade0e51ca0f1aa7c568d6216876db57b436e045951d470ef
```

该哈希在 `tq_test` 被打开**之前**生成。`unbound_fields` 为空，即所有必须在见到 test 结果前定死的字段都已绑定。

## 📌 已绑定的内容

| 字段 | 值 |
| --- | --- |
| `split.manifest_sha256` | `7548dd7fcabe7a9fd6e67ddeae2bf095c7f24027c7c083e3650efdacf371e6cd` |
| `split.retained` | 338,945 |
| `split.sample_counts` | train 186,420 / validation 50,841 / tq_test 50,841 / bridge_holdout 50,843 |
| `panel.shape` | 20,801 variants × 338,945 samples |
| `panel.chromosomes_included` | chr1–17, 19, 20, 22（20 条） |
| `panel.chromosomes_excluded` | chr18, chr21 |
| `traits.selected` | `n_30020_0_0`、`n_30010_0_0`、`n_30040_0_0` |
| `traits.transform` | `zscore`（train 拟合的 ±5 SD 截尾） |
| `analysis.eligibility_rule` | `delta_r2 > 0 and paired 95% CI lower bound > 0` |
| `analysis.sesoi_delta_r2` | 0.005（与 eligibility 分离） |
| `analysis.p_threshold_grid` | 9 个阈值，只在 validation 上选一个 |
| `firewall.test_opened` | `false` |
| `split.bridge_holdout_opened` | `false` |

## ⚠️ 冻结时仍存在的一个缺口

`cohort.withdrawal_list_supplied = false`。本服务器上不存在参与者撤回名单，当前仅有 `b3` 剔除的 111 个负号标识作为部分覆盖，且 genotype 转换于 2020 年 1 月，其后的撤回不体现在那批标识中。

`cohort.withdrawal_resolution_required_before = "tq_test opening"`。**该缺口必须在打开 `tq_test` 之前解决**：若撤回名单后续到位并改变样本集合，`b3` 需重跑、split 哈希改变，已花掉的一次性 test 就落在被取代的划分上。详见 [`FLAG_AUDIT_RESULT_20260917.md`](FLAG_AUDIT_RESULT_20260917.md)。

## ▶️ 冻结之后允许做的事

只有 validation-only 的 `b5`（不加 `--allow-test`）。它在 `train` 上拟合、在 `validation` 上评估，并按预注册网格选定唯一的 p 值阈值。该阈值选择是协议内动作，不是「看结果调参」。

`tq_test` 与 `bridge_holdout` 在此阶段均不打开。
