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

## 1. 本次审计的证据等级（已更新）

**论文原文已读到。** 用户提供了 **v2（2025-04-11）** 全文 PDF，Methods 已逐条抽进 `PAPER_SPEC.zh-CN.md`，原先所有 `[待核实]` 条目**均已解除**。

两点仍需注意：

1. **版本差**：读到的是 v2，最新是 **v3（2026-01-15）**。v2 只做了 **TG/HDL 一个性状**；v3 摘要另称覆盖 LDL、T2D 与跨人群外推。**复现目标若包含后三项，必须再拿 v3。**
2. **Supplementary 未获取**：本会话出口策略拦截 biorxiv.org / PMC / zenodo.org，补充材料没读到。网格搜索的**具体网格点**很可能在那里（见 B6 残留）。

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

### B1 — UK Biobank 数据准入 【阻断 S3，不可绕过】[原文已核实]

- 论文主结果建立在 UK Biobank 个体级基因型与表型上：**423,888 名参与者**、**203,126 个 SNP**、限 **Caucasian ancestry**（已从 v2 原文核实）。作者使用的 UKB 申请号是 **51436 与 26041**——这是他们的，我们必须有自己的。
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

### B3 — ~~论文使用的本体与 SNP→基因映射未随代码发布~~ 【已解除】[原文已核实]

v2 Methods 把本体来源写死了，这条从"阻塞"降为"照做即可"：

- **GO Biological Process, version 2023-07-27**，下载地址原文给定（`release.geneontology.org/2023-07-27/`）；
- 用 **DDOT** 剪枝，只保留与"有显著 SNP 映射的基因"相关的 term；
- 排除 **映射基因数 < 5** 的 system；
- 再排除 **注释基因集合与某个子系统完全相同** 的 system。

SNP→基因映射也已写死：**cS2G（2023 年 11 月版）∪ GTEx v7 eQTL（7 个指定组织）∪ 最近基因（hg19）**。详见 `PAPER_SPEC.zh-CN.md` §3–§4。

**但本体仍然是数据派生物**：剪枝依据是"哪些基因有显著 SNP 映射"，而显著性来自作者自己在 UKB 上跑的 BOLT-LMM。所以本体**必须在我们自己的 GWAS 结果上重建**，不能指望拿到一份现成的。这一点决定了 Stage C 必须从 GWAS 开始，不能从本体开始。

**顺带否定一条歧路：** 上游 `Collapse_Gene_Ontology_Based_on_GWAS_results.ipynb` 下载的是 GWAS Catalog 的 `GCST90257283`——那只是流程示范，**不是论文路径**。

### B4 — 外部资源版本 【已全部锁定；但需更正我先前的一条建议】[原文已核实]

v2 的 Data Availability 把每个外部资源都给了地址与时间：

| 资源 | 论文指定版本 |
|------|------------|
| Gene Ontology | **release 2023-07-27** |
| cS2G | **2023 年 11 月**取自 `alkesgroup.broadinstitute.org/cS2G` |
| eQTL | **GTEx v7**，限 7 个组织 |
| 基因组坐标 | **GRCh37 / hg19** |

**更正：** 本审计早先写过"直接使用仓库附带的 `GO_BP_full.txt`，不要重跑下载"。**这条建议现在要反过来。** 论文锁定的是 GO 2023-07-27，而仓库附带的 `GO_BP_full.txt` 是 `go_file.ipynb` 从无版本锁的 `purl.obolibrary.org/obo/go.obo` 下载后生成的，**版本未知**。

实测该文件的规模指纹，供比对用：

```
55,030  条 term→term 边（default）
134,564 条 term→gene 边（gene）
27,596  个不同 GO term
17,775  个不同基因符号
```

**动作：** 先从 `release.geneontology.org/2023-07-27/` 取 GO 与当期 GOA，按 `go_file.ipynb` 的逻辑重建一份，与上述指纹比对。指纹一致则可用仓库附带文件；不一致则必须用重建的那份，并把差异记入偏离清单。`Homo_sapiens.GRCh37.87.gtf` 仓库内缺失，按文件名从 Ensembl release 87 取即可（版本本身是明确的）。

### B5 — 环境无法按 environment.yml 重建 【可绕，需自建最小依赖集】[已核实]

`environment.yml` 是作者机器的完整 `conda env export`，问题有四：

1. 末行 `prefix: /cellar/users/i5lee/miniconda3/envs/G2PT_github`（作者绝对路径，README 与 samples 文档里的示例命令也硬编码了这个路径）；
2. 依赖 `pytorch-nightly` 频道，含 `torchtriton=3.0.0+dedb7bdf33` 这类 nightly 构建号，**nightly 包会被下架**；
3. Python 3.8 已于 2024-10 EOL；
4. conda 段与 pip 段的 CUDA 版本不一致（conda 侧 `pytorch-cuda=12.1`，pip 侧 `torch==2.4.1+cu124` 与 `nvidia-*-cu12==12.4.*`）。

实际依赖面其实很窄（按 import 实测）：`torch` / `numpy` / `pandas` / `sklearn` / `scipy` / `tqdm` / `sgkit`（读 PLINK） / `networkx` / `prettytable` / `statsmodels`，可视化另需 `pygraphviz` + `matplotlib`。**建议自建一个锁定版本的最小环境，而不是复原 environment.yml。** 但见 B7——这个"最小环境"并不自由。

### B6 — 论文超参 【PARTIALLY CLOSED】[v2 原文已核实 + v3 二手]

v2 Methods 给出的部分：

| 项 | 论文取值 | 对应 CLI / 默认值 |
|----|---------|------------------|
| 嵌入维度 d | **64** | `--hidden-dims 64`（**不是**代码默认的 256；恰好等于 `train_model.sh` 的 demo 值）|
| 传播阶段头数 | **4** | `--n-heads 4`（与代码默认一致）|
| 翻译阶段 | **Differential Attention，1 头** | 上游 README 的 "Future work" 已勾选 Differential Transformer |
| 损失 | **MSE** | `--regression` |
| 优化器 | **AdamW**（decoupled weight decay + L2）| 代码即用 AdamW |
| 划分 | **嵌套交叉验证，train:val:test = 3:1:1** | 需自建，脚本只接受已分好的三份 bfile |
| 网格搜索维度 | **p-value 阈值** × **训练 epoch 数** | 阈值影响 SNP 集，epoch 影响取哪个 checkpoint |

**v3 补上的部分**（二手来源，用户转述自 Review Commons 作者回复；详见 `PAPER_SPEC.zh-CN.md` §10）：

- **learning rate**：grid {1e-3, 1e-4, 1e-5} → 选 **1e-5**
- **weight decay**：grid {0, 1e-5, 1e-3} → 选 **1e-3**
- **embedding dim**：grid {64, 128} → 选 **64**（与 v2 原文一致，可信）
- **CV**：**5-fold** nested，每折 60/20/20（比例与 v2 原文一致）

**仍然缺的**：`batch size`、`dropout`、`训练 epoch / 早停的正式设置`、以及论文描述的 5 个传播步骤对应代码里哪组 `--sys2env/--env2sys/--sys2gene/...` 开关。

**注意 lr=1e-5 这个值的分量**：它比代码默认低两个数量级。配合 B10 的初始化分歧，"用哪套初始化"在这个学习率下会直接决定训练轨迹。

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

### B9 — 不传 `--target-phenotype` 时，验证指标恒为 0，早停与学习率调度双双失效 【实测确认，影响最深】[已核实]

这是 Stage A 跑通后发现的，也是目前最需要论文 Methods 澄清的一条。

`src/utils/trainer/snp2p_trainer.py:321-330`：

```python
target_performance = 0.
for i, pheno in enumerate(phenotypes):
    ...
    if pheno == self.target_phenotype:
        target_performance = performance
return target_performance
```

`--target-phenotype` 的 argparse 默认值是 `None`，而表型列名是 `PHENOTYPE`，**永远匹配不上**，于是 `evaluate()` 恒返回 `0.0`。这个 `0.0` 同时喂给两个地方（`:196-197`）：

- `self.scheduler.step(performance)` → `ReduceLROnPlateau` 永远收到 0.0，**学习率永不下调**；
- `self.early_stopping(performance, model)` → `EarlyStopping` 永远收到 0.0。

再看 `src/utils/trainer/utils.py` 的 `EarlyStopping.__call__`：首次调用走 `best_score is None` 分支，把**当次权重**存为 `best_weights`；此后 `0.0 > 0.0 + 1e-5` 恒为假，一路走 `else` 分支累加 `wait`。后果有三：

1. **`best_weights` 停在第一次验证时的权重，之后再不更新**——`{out}.best` 不是"最好的"模型，而是"最早验证过的"模型。我们这次跑完全程 21 个 epoch，日志里**一次 `New best score` 都没出现**，只有 `No improvement 1/10 … 4/10`，最后照样写出了 `output_model.pt.best`。
2. **`--patience` 变成一个纯粹的步数闸门**：既然"改善"不可能发生，训练必然在第 `patience` 次验证后停止，与验证表现无关。按代码默认值（`epochs=300, val_step=20, patience=10`）算，训练总是停在 epoch 220 附近——**看起来像早停，实则是定时器**。
3. 学习率全程恒定。

**对复现的意义：** `--target-phenotype` 这一个看似无关紧要的开关，实际决定了模型选择策略、训练时长和学习率轨迹。论文的 UKB 实验有没有传这个参数，会导致两个完全不同的训练过程。**这条必须在 Methods 里确认**，它把 B6 从"抄超参"升级成了"抄完整命令行"。

注意上游自带的 `train_model.sh` **没有**传 `--target-phenotype`，因此仓库里附带的 `samples/output_model.pt.best` 大概率就是这条路径产出的。

**读到论文后的修正（重要）：** v2 Methods 写明模型选择是"**对训练 epoch 做网格搜索**"，即按验证集在若干 epoch 检查点里挑，**不依赖代码里的 EarlyStopping**。所以这个失效的早停逻辑**大概率没有污染论文结果**。但它对我们的要求反而更清晰了：

- **必须保留逐 epoch 检查点 `{out}.pt.N`，在验证集上自己选 epoch；**
- **绝不能用 `{out}.best`**，它只是第一次验证的快照；
- 复现协议里要把"epoch 选择"写成一个显式的网格搜索步骤，而不是"开早停让它自己停"。


### B10 — 论文描述的 FFN 与初始化，在公开代码里不存在 【新发现，对精准复现杀伤力等同 B2】[代码已核实]

v3 的修订说明称新增了"更详细的 HiGT initialization"。但把转述的两条规格拿去对代码，**都对不上**：

| 论文（二手转述）| 上游代码实际 | 位置 |
|---------------|------------|------|
| 两层 position-wise FFN，内部维度 4×，**GeLU** | **SwiGLU**：`w12: d→2×inner` 分成门控与值，`F.silu(u) * v`，再 `proj: inner→d`，**无 bias** | `hierarchical_transformer.py:8-20` |
| SNP/gene/system 嵌入用 **uniform Xavier** | **没有任何自定义初始化**。裸 `nn.Embedding(...)` → PyTorch 默认 **N(0,1)**；唯一的 `init_method` 是 `kaiming_uniform_(a=√5)`（即 PyTorch Linear 的默认值），且标注为 Differential Transformer 的辅助函数 | `snp2phenotype.py:82-83`、`g2pt.py:22-23`、`attention.py:8-9` |

**这不是版本差。** 我把仓库历史拉深到 227 个提交后确认：

- SwiGLU 在 **2025-04-28** 引入（`7dd9726`，实现 xFormers 那次），早于 v2 定稿；
- v3 发布时点前最后一个提交 **`af5dd26`（2026-01-13）** 的 `hierarchical_transformer.py`，FFN **已经是 SwiGLU**，与今天的 HEAD 逐字相同；
- `git log -S "xavier"` 显示 xavier 最后一次变动是 **2025-06-27 的 "remove DR and G2P model"**——即它属于**已被删除的旧模型分支**，不在当前 G2PT 路径上。

GeLU 在代码里确实存在，但在**预测头**（`g2pt.py:63`、`sys2pheno.py:25`、`geno2pheno.py:27`），不在 HiGT block 的 FFN。所以要么是论文用标准 Transformer 术语描述了一个实际是 SwiGLU 的模块，要么是二手转述不准。

#### 为什么初始化这条必须当成一级问题

`nn.Embedding` 默认是 **N(0,1)，std = 1.0**；Xavier uniform 对形状 (V, 64) 的张量给出 `bound = √(6/(64+V))`：

| 张量 | 词表 V | xavier std | **默认 / xavier** |
|------|-------|-----------|------------------|
| SNP 嵌入（p=1e-8，V=7,289）| 7,289 | 0.0165 | **61×** |
| SNP 嵌入（p=1e-5，V=16,532）| 16,532 | 0.0110 | **91×** |
| gene 嵌入 | 253 | 0.0794 | 13× |
| system 嵌入 | 20 | 0.1543 | 6× |

（代码里 `snp_embedding` 的词表是 `n_snps*3+2`，因为每个 SNP × 3 种合子性状态。）

在 **lr = 1e-5** 的 AdamW 下，初始尺度差 60–90 倍不是风格差异，是两个不同的优化问题：从 std=1 起步，注意力 logits 初始方差大、softmax 早期接近 one-hot，而 1e-5 的步长几乎不可能把嵌入拉回小尺度；从 std≈0.01 起步则是"从近似均匀注意力长出结构"。**两者在固定 epoch 预算下不会收敛到同一个解。**

#### 处置

复现分支必须把初始化做成**显式开关**（`--init {default,xavier}`），两种都跑，并在协议里记为一处 **待判定的规格分歧**，而不是默认沿用代码行为。FFN 同理：若要严格按论文文字，需要把 SwiGLU 换成 GeLU 两层 FFN，这是一处**实现级偏离**，必须记录。


---

## 4. 可以直接开工的部分

- **S1 冒烟**：`bash train_model.sh` → `bash predict_model.sh`，全部依赖仓库自带的合成数据，无需任何外部准入。仓库附带的 `samples/output_model.pt.20`、`output_model.PHENOTYPE.head_0.csv`、`*.sys_importance.csv` 可作为**参考输出**，用来判定我们的环境与作者是否一致。
- **S2 上位效应模拟**：`src/utils/analysis/epistasis_simulation.py` 的模拟函数带**固定种子**（`seed=123` / `seed=42`），配合 `evaluate_epistasis_retrieval.py` 可以在完全不碰 UKB 的情况下复现论文中"能检出上位交互"这一类声明。**这是 S3 之前最有价值的靶子。**
- **基线对照**：`src/model/other_models/` 已含 elasticnet / xgboost / dnn；`samples/run_plink_score.sh` 提供 PLINK PRS 全流程（GWAS → score → r² 评估）。论文的对照组基本不用自己实现。

---

## 5. 判定

| 复现目标 | 判定 | 前置条件 |
|---------|------|---------|
| 跑通官方代码（21 epoch 全程走完，exit 0） | **已验证** | 见 `ENVIRONMENT_NOTES.zh-CN.md` |
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

## 7. 待确认事项（已更新）

论文 PDF 已收到（v2），原先三问中的第一问解除。剩下的：

### 7.1 UKB 服务器上到底有什么 —— 请逐条对照

复现 Stage C 需要下列输入。**这个清单是按 v2 Methods 反推的**，请按"有 / 没有 / 不确定"逐条回答：

| # | 需要的东西 | 论文要求 |
|---|-----------|---------|
| 1 | 个体级**基因型**（SNP array，PLINK bed/bim/fam 或可转换格式）| 覆盖 203,126 个 SNP，坐标 **hg19/GRCh37** |
| 2 | **血清甘油三酯 TG**（mmol/L）| 用于 log₂(TG/HDL) |
| 3 | **血清 HDL 胆固醇**（mmol/L）| 同上 |
| 4 | **sex、age** | 协变量 |
| 5 | **遗传主成分 PC1–PC10** | 协变量；UKB 有现成字段，也可自算（自算即偏离）|
| 6 | **祖先/人群标签** | 论文限 Caucasian ancestry 子集 |
| 7 | 有效样本量 | 论文 423,888；我们的子集有多少？ |
| 8 | **BOLT-LMM** 是否可用（或可安装）| SNP 筛选必须用它；这一步是 Stage C 的起点 |
| 9 | 服务器算力 | GPU 型号与显存、是否多卡、是否支持 torchrun |

外部资源（与 UKB 无关，可独立准备）：cS2G、GTEx v7 eQTL（7 组织）、GO release 2023-07-27、Ensembl GRCh37.87 GTF。

> 说明：你提到"之前的会话窗口已经问过一遍"——**那次对话的内容不在本会话里，我看不到**。上面这份清单是按论文重新推的，请直接在这份上回答。

### 7.2 v3 找不到 —— 不阻塞，按 v2 开工

结论：**v2 足够支撑 TG/HDL 这条主线的复现，不要为了等 v3 停工。** v3（2026-01-15，确实存在）相对 v2 主要是把性状范围从 TG/HDL 扩到 LDL、T2D 与跨人群外推；而我们 Stage A/B/C 的规格全部来自 v2 的 Methods，已经齐了。

若仍要取，直链见 `FETCH_LIST.zh-CN.md`。

### 7.3 真正该优先取的是 v2 的 Supplementary

**它的优先级高于 v3 全文。** 现在唯一还卡着 B6 的是网格搜索的具体取值（lr、wd、batch size、dropout、epoch 网格范围、CV 折数），这类内容通常在补充材料里。而且检索显示 v2 Supplementary 还包含**消融实验**（完整 G2PT vs 去掉 system / 去掉 system+gene / 去掉层次连接 / 去掉 SNP→gene 映射 的 R² 对照）——这本身就是一组很好的复现靶子，比"把主结果数字对上"更容易验证实现是否正确。

直链在 `FETCH_LIST.zh-CN.md`，本会话出口被拦，取不到。

## 附：与本仓库主线的关系

G2PT 与本仓库的 Grassmann/spectral 主线是**同一问题的两条技术路线**（SNP → 基因/系统 → 表型）。G2PT 走的是显式生物本体 + 层次化注意力；本仓库走的是谱结构 + Grassmann 记忆。复现 G2PT 的直接收益是得到一条**可信的强对照基线**，而不只是复刻一篇论文。建议 Stage A/B 完成后，即把 G2PT 纳入本仓库既有的 gate 体系作为对照臂。
