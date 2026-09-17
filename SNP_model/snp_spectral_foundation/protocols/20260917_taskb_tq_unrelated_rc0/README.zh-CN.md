# Task B 资格检验执行包

_协议 `20260917_taskb_tq_unrelated_rc0` · 只做 Task Qualification，不开 bridge_

---

## 📋 五个脚本

| 脚本 | 作用 | 依赖 |
| --- | --- | --- |
| `b1_scan_sas_fields.py` | 找出哪个 SAS 文件含哪些 UKB 字段，**重点是 22021** | 标准库（有 pyreadstat 则用它） |
| `b2_extract_sas_columns.py` | 只把点名的列抽成 TSV | pyreadstat 或 pandas |
| `b3_build_unrelated_split.py` | `22021 == 0` 的无亲缘子集 → 分层 70/15/15 | 标准库 |
| `b4_build_snp_panel.py` | 从 `.bed` 抽薄出 SNP panel，train-only MAF/缺失过滤 | numpy |
| `b5_task_qualification.py` | A（covariates）vs B（A + raw dosage），出 `ΔR²` 与判定 | numpy |

`b1` 和 `b3` 纯标准库，你现在的系统 `python3` 就能跑。`b4`/`b5` 只要 numpy。

## 🔀 决策点在 b1

```text
b1 报告 22021 存在 ──► b3 冻结无亲缘子集 ──► b4 ──► b5
                 │
                 └─► 不存在 ──► 不要开 bridge
                                ├─ 首选：向数据管理方要官方 ukb_rel_*.dat（application 44430）
                                └─ 备选：KING_FALLBACK.zh-CN.md 自建 kinship，再走 E0 的 R2 builder
```

**随机个体划分不是选项。** UKB 约三成参与者有三度以内亲属，随机划分会把家系共享的 genotype→phenotype 成分泄漏进 test，测到的是泄漏不是泛化。

## ⚖️ 为什么 Task B 不被 imputed 阻断

E0 卡在 imputed，是因为它的**预测目标**就是 genotype——用 IMPUTE4 的输出当 target，而最强基线又是同类 HMM，构成循环。

Task B 的目标是 **phenotype**，genotype 只是**输入**。imputed dosage 作为输入是全领域标准做法。所以 array 直接分型和 phased haplotype 这一轮都不需要。

## 🔒 Firewall

- 不加 `--allow-test` 时 `b5` 只出 validation 结果
- 开过 test 会写 `TEST_OPENED.marker`，同目录再开直接拒绝
- covariate 系数、SNP 效应、缺失填补、阈值选择全部只用 train/validation
- 本包产出的 phenotype 文件**只属于 Task B**，E0 不得读取
- summary 只含计数与哈希，不含参与者标识；`tq_split_manifest.tsv` 和抽取出的 phenotype TSV 含标识，必须留在服务器上

## ⚠️ 三条必须随结果一起报告的限制

1. **Allele 方向未验证** — `.bed` 由 BGEN 转换时没指定 `--bgen ref-first`，所以 dosage 定义为「`.bim` A1 计数」，全流程不与外部 reference 比对
2. **Panel 是稀疏抽样** — 按物理距离抽薄，`ΔR²` 是资格信号的下界，不是最优预测性能，也不是可捕获遗传度的上界
3. **relatedness 强度继承自 UKB** — `22021` 是每人一个汇总值，不是配对边；结论范围是「新的无亲缘个体」而非「新家系」

## 🧪 自测

```bash
python3 -m unittest discover -s tests -v
```

7 个合成测试：用真实 PLINK 二进制格式模拟一个有已知遗传成分的小队列，跑通 `b3 → b4 → b5`，验证亲缘个体被剔除、split 全为 singleton 且比例正确、`.bed` 字节解码正确、**有信号的 trait 必须 QUALIFIED**、**纯噪声 trait 必须不 QUALIFIED**、不给 `--allow-test` 时 test 保持关闭、同目录重开 test 被拒绝。

## 📄 文档

- [`PROTOCOL.zh-CN.md`](PROTOCOL.zh-CN.md) — 冻结的估计量、两个臂、SESOI、判定与 firewall
- [`SERVER_STEPS.TaskB.zh-CN.md`](SERVER_STEPS.TaskB.zh-CN.md) — 逐条命令
- [`KING_FALLBACK.zh-CN.md`](KING_FALLBACK.zh-CN.md) — 22021 不存在时的自建 kinship 路径
- [`CONFIG.template.json`](CONFIG.template.json) — 开 test 前必须填的预注册字段
