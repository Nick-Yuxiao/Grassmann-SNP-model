# Task B 字段绑定记录

_2026-09-17 · `b1_scan_sas_fields.py` 于 `c20-5090` 的扫描结果 · 状态 `RELATEDNESS_AXIS_AVAILABLE`_

---

## ✅ 关键结论

UKB field `22021`（genetic kinship to other participants）**存在**。Task B 的 relatedness 轴因此可用，无亲缘子集方案可以冻结，`KING_FALLBACK` 路径不启用。

## 📂 各 SAS 文件的字段覆盖

| 文件 | 大小 | 字段数 | `22021` | `22006` | `22027` | 血常规 `30020/30010/30040` | 血脂 `30780/30760/30870/30690` |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `ukb_alldata.sas7bdat` | 114.9 GiB | 7,385 | ✅ | ✅ | ✅ | ✅ | ✅ |
| `new.sas7bdat` | 88.1 GiB | 4,754 | ✅ | ✅ | ✅ | ✅ | ✅ |
| `ukb3.sas7bdat` | 57.5 GiB | 3,373 | ❌ | ❌ | ❌ | ❌ | ✅ |
| `ukb2.sas7bdat` | 52.7 GiB | 1,551 | ✅ | ✅ | ✅ | ✅ | ❌ |
| `ukb1.sas7bdat` | 28.3 GiB | 1,530 | ❌ | ❌ | ❌ | ✅ | ✅ |
| `ukb5.sas7bdat` | 3.8 GiB | 348 | ❌ | ❌ | ❌ | ✅ | ✅ |
| `data1.sas7bdat` | 1.9 GiB | 253 | ❌ | ❌ | ❌ | ❌ | ❌ |
| `ukb4.sas7bdat` | 1.9 GiB | 252 | ❌ | ❌ | ❌ | ❌ | ❌ |
| `data.sas7bdat` | 1.2 GiB | 64 | ❌ | ❌ | ❌ | ❌ | ❌ |
| `ukb_event.sas7bdat` | 3.5 GiB | 16 | ❌ | ❌ | ❌ | ❌ | ❌ |

## 📉 最小读取量的抽取计划

| 抽取内容 | 选定文件 | 大小 | 理由 |
| --- | --- | --- | --- |
| Trait（7 个候选全覆盖） | `ukb5.sas7bdat` | **3.8 GiB** | 含全部 7 个候选 trait 的最小文件 |
| `22021` + QC flags | `ukb2.sas7bdat` | **52.7 GiB** | 含 `22021`/`22006`/`22027` 的最小文件 |

相比统一从 `ukb_alldata.sas7bdat`（114.9 GiB）抽取，读取量降低约一半，trait 一侧降低约 30 倍。

## ⚠️ 扫描方法的边界

本次扫描在无 pyreadstat 的环境下运行，方法为 `byte_scan_heuristic`：只读每个文件的前 300 MB 并匹配列名文本。

- **`✅` 可信**：列名确实出现在文件元数据中
- **`❌` 不构成不存在的证明**：只说明该字段未出现在被扫描的前 300 MB 内

因此上表的「最小文件」选择是基于阳性发现的，安全；但不能据此断言某文件「没有」某字段。装上 pyreadstat 后重跑 `b1` 会自动改用 `pyreadstat_metadata` 方法给出权威列表。

## ❓ 抽取前仍需确认

1. **标识列名**：`b1` 只匹配 `n_<field>_<instance>_<array>` 形式，未覆盖标识列。需确认 SAS 中是 `n_eid`、`eid` 还是其它写法。
2. **参与者覆盖**：`ukb5` 与 `ukb2` 是否覆盖同一套 487,409 人，尚未验证。`b3` 的 `dropped_no_kinship_flag` 计数会直接暴露覆盖缺口。
3. **instance 选择**：trait 列可能有多个 instance（`_0_0`、`_1_0`）。本轮固定使用 baseline instance `_0_0`，写入 `CONFIG.json`。

## ▶️ 下一步

1. 装 numpy 与 pyreadstat（`b4`/`b5` 与 `b2` 的依赖）
2. 用 `--row-limit 200` 冒烟确认标识列名
3. 从 `ukb5` 抽 trait、从 `ukb2` 抽 `22021` 与 QC flags
4. `b3` 冻结无亲缘子集，检查 `dropped_no_kinship_flag` 与 `dropped_related_or_unassessed`
