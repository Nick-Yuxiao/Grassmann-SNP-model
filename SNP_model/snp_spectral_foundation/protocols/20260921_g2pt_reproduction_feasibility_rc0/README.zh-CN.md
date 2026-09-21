# 20260921_g2pt_reproduction_feasibility_rc0

G2PT（*A genotype-phenotype transformer to assess and explain polygenic risk*，bioRxiv `10.1101/2024.10.23.619940`）**精准复现的可行性审计**。本目录只做判定，不含训练代码。

## 文件

| 文件 | 内容 |
|------|------|
| `FEASIBILITY_AUDIT.zh-CN.md` | 主审计报告：结论、证据等级、资产盘点、六项阻塞、三阶段路径 |
| `UPSTREAM_INVENTORY.json` | 机读清单：上游 commit、实测数字、超参三方对照、阻塞项、待你确认的问题 |

## 一句话结论

上游代码（`idekerlab/G2PT` @ `fc13abd`，MIT）完整且质量高，合成数据端到端链路可立即跑通；但**训练脚本没有随机种子**（已在源码中实测确认），因此 UKB 主结果的数值在原理上不可 bit-exact 复现，复现口径须改为多种子统计等价。UKB 数据准入是日历时间上的最大不确定项。

## 状态

- 本次审计**未能读到论文原文**：本会话出口策略拦截 biorxiv.org / pmc.ncbi.nlm.nih.gov / zenodo.org。凡报告中标记 `[待核实]` 的（样本量、SNP 数、论文超参、精度数字）均需 PDF 才能解除。
- 上游仓库已在本会话做只读浅克隆用于审计，**未纳入本仓库版本管理**。
