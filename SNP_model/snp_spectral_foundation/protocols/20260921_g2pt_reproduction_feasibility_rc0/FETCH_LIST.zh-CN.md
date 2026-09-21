# 待下载清单

本会话的网络出口策略拦截 biorxiv.org / PMC / zenodo.org，以下文件**必须由你在能正常上网的机器上取**，再放进仓库或服务器。

## A. 论文相关

| 文件 | 地址 | 放置路径 | 优先级 |
|------|------|---------|--------|
| **v2 Supplementary**（含消融实验：完整 G2PT vs 去掉 system / 去掉 system+gene / 去掉层次连接 / 去掉 SNP→gene 映射 的 R² 对照）| `https://www.biorxiv.org/content/biorxiv/early/2025/04/11/2024.10.23.619940/DC1/embed/media-1.pdf?download=true` | `paper/G2PT_v2_supplementary.pdf` | **高** |
| v3 全文 PDF（2026-01-15） | `https://www.biorxiv.org/content/biorxiv/early/2026/01/15/2024.10.23.619940.full.pdf` | `paper/G2PT_biorxiv_v3_20260115.pdf` | 中 |
| v3 网页版（若 PDF 直链失效，从这里点 Download PDF） | `https://www.biorxiv.org/content/10.1101/2024.10.23.619940v3.full` | — | 备用 |
| v3 Supplementary | **从 v3 网页的 "Supplementary Material" 按钮进入**。把 DC1 路径的日期换成 `2026/01/15` 只是**未经验证的猜测**，可以试，但别当成确定地址 | `paper/G2PT_v3_supplementary.pdf` | 中 |

| **Review Commons 作者回复 / Sciety**（已公开 Supplementary Fig. 8 的超参数值）| `https://sciety.org/articles/activity/10.1101/2024.10.23.619940` | 摘录已入 `PAPER_SPEC.zh-CN.md` §10 | **高** |

> **v2 Supplementary 优先级高于 v3 全文。** 我们现在缺的是网格搜索的具体取值（lr / wd / batch size / dropout / epoch 网格 / CV 折数），这些属于补充材料的内容；而 v3 相对 v2 主要是**扩了性状范围**（加 LDL、T2D、跨人群），对 TG/HDL 这条主线的规格影响有限。

## B. 论文指定的外部资源（Data Availability 原文给出）

| 资源 | 地址 | 说明 |
|------|------|------|
| **Gene Ontology 2023-07-27** | `https://release.geneontology.org/2023-07-27/index.html` | 取 `go.obo` 与当期 `goa_human.gaf.gz`，按 `go_file.ipynb` 逻辑重建 `GO_BP_full.txt` |
| **cS2G** | `https://alkesgroup.broadinstitute.org/cS2G` | 论文取的是 **2023 年 11 月**版本；若该站点只提供最新版，需记为偏离项 |
| **GTEx v7 eQTL** | `https://www.gtexportal.org/home/downloads/adult-gtex/qtl` | 只需 7 个组织：adipose subcutaneous、adipose visceral omentum、liver、pancreas、adrenal gland、muscle skeletal、uterus |
| **Ensembl GRCh37 release 87 GTF** | Ensembl GRCh37 归档区 release-87 的 `Homo_sapiens.GRCh37.87.gtf.gz`（路径请在站点上核对，不要照抄我拼的 URL）| 上游 notebook 引用但仓库内缺失 |

## C. 工具

| 工具 | 用途 | 检查点 |
|------|------|--------|
| **BOLT-LMM** | SNP 筛选（Stage C 的第一步） | 服务器上是否已装 / 能否装 |
| **PLINK** | PRS 对照组（上游 `samples/run_plink_score.sh` 直接用） | 同上 |
| **DDOT** | GO 本体剪枝 | Python 包，需确认与我们的 Python 版本兼容 |

## 取到之后

1. `paper/` 下的 PDF **不入库**（`.gitignore` 已忽略 `paper/*.pdf`），预印本页脚写明 All rights reserved。
2. 把 Supplementary 给我，我补齐 `batch size / dropout / epoch 设置`，B6 才算彻底关闭（lr / wd / dim / heads / CV 已由 v3 二手来源锁定，见 `PAPER_SPEC.zh-CN.md` §10）。
2b. **特别需要核对的**：v3 正文关于 HiGT initialization 与 FFN 的原话——公开代码与转述规格冲突，见 B10。
3. GO 重建后与 `FEASIBILITY_AUDIT.zh-CN.md` B4 里记的指纹比对（55,030 / 134,564 / 27,596 / 17,775）。
