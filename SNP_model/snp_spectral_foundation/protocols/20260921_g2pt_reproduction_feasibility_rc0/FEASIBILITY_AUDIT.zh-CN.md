# G2PT 精准复现：可行性审计（rc0）

**日期：** 2026-09-21
**目标论文：** *A genotype-phenotype transformer to assess and explain polygenic risk*，bioRxiv `10.1101/2024.10.23.619940`
**上游代码：** <https://github.com/idekerlab/G2PT>，commit `fc13abd59d31b9567b418b152cee189de358b095`（2026-04-23），MIT License
**审计范围：** 只回答"能不能精准复现、卡在哪"，不写训练代码。

---

## 0. 结论先行

**代码级复现可行；数字级（bit-exact）复现不可行。** 复现目标必须先拆成三层，因为这三层的阻塞完全不同：

| 层 | 内容 | 判定 | 主要阻塞 |
|----|------|------|----------|
| **S1 冒烟层** | `samples/` 合成数据端到端跑通训练+预测+注意力导出 | **已验证可跑**（CPU 上需一行补丁） | 环境重建（B5/B7）+ 一处漏写的 CUDA 守卫（B8） |
| **S2 模拟层** | 上位效应（epistasis）模拟与检出评估 | **可行，且有固定随机种子** | 论文所用模拟参数未知（B6） |
| **S3 主结果层** | UK Biobank 上 TG/HDL、LDL、T2D 的预测精度、跨人群外推、系统级注意力解释 | **当前不可行** | UKB 数据准入（B1）+ 训练无种子（B2）+ 本体非常量（B3） |

对"精准"二字的直接回答：**即使拿到 UK Biobank 数据，也无法复现出论文的同一组数字**，因为上游训练脚本根本没有设随机种子（见 B2，已在源码中实测确认）。可达到的最强口径是**统计等价复现**——多种子重复训练，报告均值与置信区间，检验论文结论是否落在区间内。这个口径的变更建议在动手前就写进复现协议，否则后期一定会被"为什么 R² 差 0.003"这类问题反复消耗。

---

## 1. 本次审计的证据等级（重要）

本会话的网络出口策略**拦截了 biorxiv.org、pmc.ncbi.nlm.nih.gov、zenodo.org**，因此：

- **[已核实]** 标记的事实 = 直接读上游源码/文件得出，可逐条复查。
- **[待核实]** 标记的事实 = 来自检索摘要，**未读过论文原文**。凡涉及样本量、SNP 数、超参、精度数字的，一律为待核实。

补齐待核实项只需一件事：把论文 PDF（建议 v3，2026-01-15）放进本仓库，我逐条对齐 Methods。

---

## 2. 上游代码资产盘点 [已核实]

```
src/            12,035 LOC   模型、层次化 Transformer、数据集、训练器、可视化、上位效应分析
顶层脚本         1,847 LOC   train_snp2p_model.py / predict_attention.py /
                            train_greedy_phenotype_selection_model.py / evaluate_epistasis_retrieval.py
GO_files/GO_BP_full.txt      189,594 行（全量 GO Biological Process 边 + 基因注释）
GO_files/goID_2_name.tab      29,736 行
samples/                     合成 PLINK 三分割（train 1000 / val 500 / test 300，293 SNP）
                             + 参考 checkpoint（pt.0/5/10/15/20/best）+ 参考输出 CSV
notebooks                    GO 构建、GO 坍缩、Sankey、本体高亮、上位效应发现/模拟
```

模型骨架（`src/model/`）：`snp2phenotype.py`、`g2pt.py`、`hierarchical_transformer/{attention,geno2pheno,sys2pheno,hierarchical_transformer}.py`，另含基线 `other_models/{elasticnet,xgboost,dnn}.py`——**基线也在仓库里，对照实验不必自己重写**。

**这是一个质量明显高于平均水平的复现起点**：MIT 许可、全链路脚本、readthedocs 文档、合成样例数据，甚至附带了参考 checkpoint 和参考输出 CSV，可以当作环境一致性的 oracle 使用。

---

## 3. 阻塞项清单（按严重度排序）

### B1 — UK Biobank 数据准入 【阻断 S3，不可绕过】[部分已核实]

- 论文主结果建立在 UK Biobank 个体级基因型与表型上（检索摘要称 423,888 名参与者、203,126 个 SNP；**待核实**）。
- 上游仓库 `samples/README.md` 明确声明 [已核实]：*"No real participant data should be distributed through this repository."* 即**上游不会、也不能提供数据**。
- UKB 访问需要正式申请、机构承诺书、审批周期与费用。这是日历时间上的最大不确定项，且不是技术问题。

### B2 — 训练脚本没有随机种子 【决定"精准"能否成立】[已核实]

实测结论：

- `train_snp2p_model.py` 的 argparse 中**不存在 `--seed`**。
- 全仓库 `torch.manual_seed` 只出现在两处推理脚本：`predict_attention.py:638`、`predict_attention_with_selected_phenotypes.py:261`（均为 `seed = 0`）。
- 训练路径上没有任何 `manual_seed` / `cudnn.deterministic` / `use_deterministic_algorithms` 调用。

后果：权重初始化、dropout mask、DataLoader 打乱顺序均不受控。**同一份数据、同一组超参、同一台机器，两次训练的结果不会相同。** 因此论文中的单点数字在原理上就不可精确复现。

处置建议（写进复现协议）：
1. 在复现分支给训练脚本补 `--seed`，并同时固定 `torch`/`numpy`/`random`/`DataLoader` worker 种子；
2. 复现口径改为 **n≥5 种子重复 → 报告 mean ± 95% CI**；
3. 判定标准改为"论文点估计是否落入我们的 CI"，而非"数字是否相等"。

### B3 — 论文使用的本体与 SNP→基因映射未随代码发布 【阻断 S3 的可比性】[已核实]

- 仓库内只有两样东西：全量 `GO_BP_full.txt`，以及**合成**的 `samples/ontology.txt`、`samples/snp2gene.txt`。
- 论文用的是"按 GWAS 结果坍缩后的本体"。坍缩流程在 `Collapse_Gene_Ontology_Based_on_GWAS_results.ipynb` 中，其输入是一份 GWAS 摘要统计——notebook 里示范下载的是 GWAS Catalog 的 `GCST90257283`（公开可取，**这一点是好消息**）。
- 但论文正文的坍缩很可能基于作者自己在 UKB 上跑出的 GWAS。**本体不是常量，而是数据派生物**：GWAS 不同 → 本体不同 → 模型结构不同 → 结果不可比。这条必须在读 Methods 时首先确认：论文的本体究竟来自哪一份 GWAS。
- SNP→基因映射据称用三条证据（eQTL / cS2G / 最近基因）合并（**待核实**），每条都需要定位到具体版本的外部资源。

### B4 — 外部资源缺失与版本漂移 【可绕，但必须先锁版本】[已核实]

- `Collapse_Gene_Ontology_Based_on_GWAS_results.ipynb` 引用 `../G2PT/GO_files/Homo_sapiens.GRCh37.87.gtf`，**该文件不在仓库内**（需从 Ensembl 取；文件名已指明 GRCh37 release 87，版本是明确的）。
- `go_file.ipynb` 从 `http://purl.obolibrary.org/obo/go.obo` 下载本体，**该 URL 永远指向最新版 GO**，无版本锁。所幸仓库已附带生成好的 `GO_BP_full.txt`，应当**直接使用仓库附带文件，不要重跑下载**，否则 GO 版本漂移会静默改变本体。

### B5 — 环境无法按 environment.yml 重建 【可绕，需自建最小依赖集】[已核实]

`environment.yml` 是作者机器的完整 `conda env export`，问题有四：

1. 末行 `prefix: /cellar/users/i5lee/miniconda3/envs/G2PT_github`（作者绝对路径，README 与 samples 文档里的示例命令也硬编码了这个路径）；
2. 依赖 `pytorch-nightly` 频道，含 `torchtriton=3.0.0+dedb7bdf33` 这类 nightly 构建号，**nightly 包会被下架**；
3. Python 3.8 已于 2024-10 EOL；
4. conda 段与 pip 段的 CUDA 版本不一致（conda 侧 `pytorch-cuda=12.1`，pip 侧 `torch==2.4.1+cu124` 与 `nvidia-*-cu12==12.4.*`）。

实际依赖面其实很窄（按 import 实测）：`torch` / `numpy` / `pandas` / `sklearn` / `scipy` / `tqdm` / `sgkit`（读 PLINK） / `networkx` / `prettytable` / `statsmodels`，可视化另需 `pygraphviz` + `matplotlib`。**建议自建一个锁定版本的最小环境，而不是复原 environment.yml。** 但见 B7——这个"最小环境"并不自由。

### B6 — 论文超参未在代码中固化 【需要 PDF 才能解除】[已核实]

仓库里存在两套互相矛盾的参数，且都**不是**论文的 UKB 配置：

| 来源 | hidden | heads | lr | wd | dropout | batch | epochs | patience |
|------|--------|-------|-----|-----|---------|-------|--------|----------|
| `train_model.sh`（合成 demo） | 64 | – | 1e-4 | 1e-4 | 0.2 | 128 | 21 | – |
| 代码默认值（`ModelConfig`/`TrainerConfig`） | 256 | 4 | 1e-3 | 1e-3 | 0.2 | 128 | 300 | 10 |

论文的 UKB 实验用哪一组、以及 `--sys2env/--env2sys/--sys2gene/--gene2pheno/--sys2pheno/--mlm/--use_moe/--z-weight/--cov-effect` 这些开关的取值，**只能从 Methods 里读**。这是现在最需要 PDF 的一条。

### B7 — 训练路径有四个未文档化的硬依赖，且 torch 版本被反向锁死 【实测踩到】[已核实]

按 README，装环境只有 `conda env create -f environment.yml` 一条路，而这条路因 B5 大概率走不通。自建最小环境时按顺序连撞四个**模块级无条件导入**，每个都让脚本在 import 阶段直接崩：

| # | 缺失模块 | 触发位置 |
|---|---------|---------|
| 1 | `xformers` | `hierarchical_transformer.py:4`（`LD_infuser/LDRoBERTa.py:6` 同）|
| 2 | `mlflow` | `src/utils/trainer/snp2p_trainer.py:4` |
| 3 | `obonet` | `src/utils/tree/tree.py:9` |
| 4 | `transformers` | `src/utils/data/dataset/tokenizers.py:5` |

四个都只存在于 `environment.yml` 的 pip 段，README 正文一字未提。AST 扫描得到的完整第三方导入面是 **23 个模块**（清单见 `ENVIRONMENT_NOTES.zh-CN.md`）。

第 1 条的连带后果更重要：安装 `xformers==0.0.28.post1` 会**把 torch 自动降级到 2.4.1** 并装上 `triton==3.0.0`。这反过来**验证了作者环境的 `xformers 0.0.28.post1` + `torch 2.4.1` 是自洽的**（此前对这一点的存疑可以撤销），唯一差别是 CUDA 构建（PyPI 给 `+cu121`，作者记录 `+cu124`）。

**复现协议必须把 torch 与 xformers 当成一对绑定版本冻结**，不能各自独立选版本。

### B8 — 训练脚本无条件要求 CUDA，但这是一行漏写的守卫 【已定位，可一行修复】[已核实]

四个依赖补齐后 import 全通，脚本进入 `main()`，停在：

```
File "train_snp2p_model.py", line 239, in main
    torch.cuda.set_device(args.local_rank)
RuntimeError: Found no NVIDIA driver on your system.
```

关键在于**同一文件的另外两处 CUDA 调用都有守卫**：`:114-115` 是 `if torch.cuda.is_available():`，`:286` 是 `elif args.world_size == 1 and torch.cuda.is_available():`。只有 `:239` 漏了。补上这一行守卫后，**合成链路在纯 CPU 上确实跑起来了**（实测：PLINK 数据正常载入，293 variants / 500 samples，epoch 1 train_loss=0.5498）。

另外 `--cuda` 这个 CLI 参数已经退化——它只触发一条 warning，真正的设备由 `ddp_setup()` + `set_device` 决定；`main()` 里原来的非分布式分支被整块注释掉，脚本已重构为 torchrun-only。

**后果：** 照原样，连 293 个 SNP 的合成 demo 都需要一张 CUDA 卡。排期前必须先确认硬件；复现分支应把这行守卫补上并记为一处**显式偏离**。


---

## 4. 可以直接开工的部分

- **S1 冒烟**：`bash train_model.sh` → `bash predict_model.sh`，全部依赖仓库自带的合成数据，无需任何外部准入。仓库附带的 `samples/output_model.pt.20`、`output_model.PHENOTYPE.head_0.csv`、`*.sys_importance.csv` 可作为**参考输出**，用来判定我们的环境与作者是否一致。
- **S2 上位效应模拟**：`src/utils/analysis/epistasis_simulation.py` 的模拟函数带**固定种子**（`seed=123` / `seed=42`），配合 `evaluate_epistasis_retrieval.py` 可以在完全不碰 UKB 的情况下复现论文中"能检出上位交互"这一类声明。**这是 S3 之前最有价值的靶子。**
- **基线对照**：`src/model/other_models/` 已含 elasticnet / xgboost / dnn；`samples/run_plink_score.sh` 提供 PLINK PRS 全流程（GWAS → score → r² 评估）。论文的对照组基本不用自己实现。

---

## 5. 判定

| 复现目标 | 判定 | 前置条件 |
|---------|------|---------|
| 跑通官方代码（训练链路启动、数据载入、loss 下降） | **已验证** | 见 `ENVIRONMENT_NOTES.zh-CN.md` |
| 产出与作者参考输出**一致**的合成结果 | **未验证** | 需贴近作者版本的环境 + 与 `samples/` 参考输出比对 |
| 复现上位效应检出类声明 | **大概率可行** | 论文 Methods 中的模拟参数 |
| 复现 UKB 主结果的**数值** | **不可行** | 违反 B2，原理上不成立 |
| 复现 UKB 主结果的**结论**（统计等价） | **取决于 B1** | UKB 准入 + 多种子协议 + 本体来源确认 |

---

## 6. 建议路径（三阶段，每阶段带 gate）

**Stage A｜环境与冒烟（不依赖任何审批）**
锁定最小依赖版本 → 跑 `train_model.sh` + `predict_model.sh` → 与 `samples/` 附带参考输出比对。
*Gate A：* 我们的输出与作者参考输出在容差内一致，则环境可信；不一致则先定位环境差异，**不要带着未知环境差异进入 Stage B**。

**Stage B｜模拟层复现（不依赖 UKB）**
读 Methods 取模拟参数 → 跑 `epistasis_simulation` + `evaluate_epistasis_retrieval` → 对齐论文的检出率/富集类指标。
*Gate B：* 模拟层指标落在论文声明的范围内。

**Stage C｜UKB 主结果（依赖 B1 解除）**
本体来源确认 → SNP→基因映射重建 → 多种子训练 → 统计等价判定。
*Gate C：* 论文点估计落入我们的多种子 CI。

---

## 7. 需要你提供的三件事

1. **论文 PDF**（建议 v3，2026-01-15）放到本目录下的 `paper/` 或仓库任意位置并告诉我路径。本会话无法访问 bioRxiv / PMC / Zenodo（出口策略拦截），所有"待核实"项都靠它解除。
2. **UKB 申请状态**：已有 application、正在申请、还是完全没有。这直接决定 Stage C 是排期还是搁置。
3. **目标硬件**：GPU 型号与显存。上游 README 自述稀疏注意力下 32GB 显存可支持 batch 256；你的机器决定 Stage C 要不要改批大小（改了就得在复现协议里记为偏离项）。

---

## 附：与本仓库主线的关系

G2PT 与本仓库的 Grassmann/spectral 主线是**同一问题的两条技术路线**（SNP → 基因/系统 → 表型）。G2PT 走的是显式生物本体 + 层次化注意力；本仓库走的是谱结构 + Grassmann 记忆。复现 G2PT 的直接收益是得到一条**可信的强对照基线**，而不只是复刻一篇论文。建议 Stage A/B 完成后，即把 G2PT 纳入本仓库既有的 gate 体系作为对照臂。
