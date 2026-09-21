# 20260921_g2pt_reproduction_feasibility_rc0

G2PT（*A genotype-phenotype transformer to assess and explain polygenic risk*，bioRxiv `10.1101/2024.10.23.619940`）的**复现可行性审计**与**基线规格冻结**。

## 文件

| 路径 | 内容 |
|------|------|
| `FEASIBILITY_AUDIT.zh-CN.md` | 主审计：结论、证据等级、资产盘点、B1–B11 阻塞项、三阶段路径 |
| `PAPER_SPEC.zh-CN.md` | 论文规格逐条抽取（v2 原文 §1–§9；v3 二手超参 §10）|
| `ENVIRONMENT_NOTES.zh-CN.md` | 环境重建实测：四个未文档化硬依赖、torch/xformers 版本互锁、Stage A 冒烟结果 |
| `FETCH_LIST.zh-CN.md` | 待下载清单（本会话出口被拦，需在能上网的机器上取）|
| `UPSTREAM_INVENTORY.json` | 机读清单：上游 commit、实测数字、超参三方对照、全部阻塞项 |
| `docs/G2PT_BASELINE_SPEC_DECISION.md` | **证据账本与冻结决定**：论文与代码冲突时怎么办 |
| `docs/DEEP_BRANCH_DESIGN_MEMO.zh-CN.md` | Deep branch 设计推演：参数预算、功效表、C1/C2/C3 × baseline_variant |
| `configs/g2pt/CONFIG_FROZEN.json` | 冻结基线配置，**每个字段带证据等级** |
| `experiments/B10_B11_spec_discrepancy_matrix.yaml` | B10/B11 的 2×2 析因矩阵与三阶段诊断计划 |
| `experiments/probe_init_step0.py` | step-0 诊断探针（阶段 1，已执行）|

## 结论摘要

- 上游代码（`idekerlab/G2PT` @ `fc13abd`，MIT）完整且质量高；合成数据链路**已在纯 CPU 上跑通**（21 epoch，exit 0），前提是补一行漏写的 CUDA 守卫（B8）。
- **训练脚本没有随机种子**（B2），UKB 主结果的数值在原理上不可 bit-exact 复现；口径须改为多种子统计等价。
- 论文规格已从 v2 原文读全（数据、SNP 筛选、本体、架构、训练、上位效应）；lr/wd/dim/heads/CV 由 v3 二手来源补齐（B6 `PARTIALLY_CLOSED`）。
- **论文描述与公开代码在两处不一致**：初始化（B10）与 FFN（B11）。已拆分建档，做成 2×2 析因，不预设哪一边为准。

## 状态

- 未读到 v2/v3 的 Supplementary，也未读到 v3 正文（出口策略拦截 bioRxiv / PMC / Zenodo）。
- 上游仓库在本会话做只读克隆用于审计，**未纳入本仓库版本管理**。
- 论文 PDF 不入库（预印本声明 All rights reserved），`paper/` 仅存说明与命名约定。
