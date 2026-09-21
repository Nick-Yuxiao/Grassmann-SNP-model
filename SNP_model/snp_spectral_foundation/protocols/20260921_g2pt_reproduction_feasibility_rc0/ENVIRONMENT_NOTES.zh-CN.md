# Stage A 环境重建实测记录

**目的：** 验证"不用 `environment.yml`、改用自建最小环境"这条路能否跑通上游 `train_model.sh` 的合成数据链路。
**平台：** Linux x86_64，**无 GPU**（纯 CPU），Python 3.11.15，venv（非 conda）。
**上游：** `idekerlab/G2PT` @ `fc13abd`，只读浅克隆，不纳入本仓库版本管理。

> ⚠️ **本记录不能替代作者环境的一致性验证。** 这里只回答"代码能不能跑起来、缺什么"，**不回答"结果是否与作者一致"**。后者必须靠与 `samples/` 附带的参考输出比对，且应在贴近作者版本的环境里做。

## 依赖链：连撞四次未文档化的硬导入

README 的安装说明只有 `conda env create -f environment.yml` 一条路（因 B5 大概率走不通）。自建最小环境时，按顺序撞到四个**模块级无条件导入**，每一个都会让训练脚本在 import 阶段直接崩：

| # | 缺失模块 | 触发位置 | 是否在 README 正文出现 |
|---|---------|---------|---------------------|
| 1 | `xformers` | `src/model/hierarchical_transformer/hierarchical_transformer.py:4`（`LD_infuser/LDRoBERTa.py:6` 同）| 否，只在 `environment.yml` 的 pip 段 |
| 2 | `mlflow` | `src/utils/trainer/snp2p_trainer.py:4` | 否 |
| 3 | `obonet` | `src/utils/tree/tree.py:9` | 否 |
| 4 | `transformers` | `src/utils/data/dataset/tokenizers.py:5` | 否 |

AST 扫描 `src/` + 四个顶层脚本得到的**完整第三方导入面共 23 个**：

```
dash google_genai matplotlib mlflow networkx numpy obonet openai pandas plotly
prettytable pygraphviz scipy seaborn sgkit sklearn statsmodels torch tqdm
transformers xformers xgboost yaml
```

（`openai` / `google_genai` 属 LLM 解释模块，`dash` / `plotly` / `seaborn` / `pygraphviz` 属可视化，`xgboost` 属基线——训练主路径不需要，但复现全部图表需要。）

## torch 与 xformers 的版本互锁（已验证）

装 `xformers==0.0.28.post1`（作者 pip 段里的版本）时，pip **自动把 torch 从 2.14.0 降级到 2.4.1**，并装上 `triton==3.0.0`：

```
Successfully installed ... torch-2.4.1 triton-3.0.0 xformers-0.0.28.post1
>>> torch 2.4.1+cu121  xformers 0.0.28.post1
```

**结论：作者环境里 `xformers==0.0.28.post1` + `torch==2.4.1` 这一对是自洽的**（此前审计里对这一点的存疑可以撤销）。唯一差别是 CUDA 构建：PyPI 默认给 `+cu121`，作者记录的是 `+cu124`。

**对复现协议的要求：torch 不是自由变量，必须与 xformers 作为一对绑定版本冻结。**

## 阻断点：训练脚本无条件要求 CUDA（→ 审计 B8）

四个依赖补齐后，import 全部通过，脚本进入 `main()`，然后停在：

```
File "train_snp2p_model.py", line 239, in main
    torch.cuda.set_device(args.local_rank)
RuntimeError: Found no NVIDIA driver on your system.
```

关键细节：**同一文件里另外两处 CUDA 调用都做了保护**——`ddp_setup()` 内的 `train_snp2p_model.py:114-115` 是 `if torch.cuda.is_available():`，`:286` 是 `elif args.world_size == 1 and torch.cuda.is_available():`。只有 `:239` 这一处漏了守卫。

所以这不是"设计上只支持 GPU"，而是**一处漏写的 `is_available()` 守卫**挡住了代码其余部分已经支持的 CPU 路径。另外 `--cuda` 这个 CLI 参数实际已经退化：它只触发一条 warning，真正的设备由 `ddp_setup()` + `set_device` 决定（`main()` 里非分布式分支的原始代码被整块注释掉了，脚本已重构为 torchrun-only）。

**后果：连 293 个 SNP 的合成 demo 都需要一张 CUDA 卡**，除非打这一行补丁。已在只读克隆里做过探针式验证（**未提交到任何仓库**，原文件已备份）。

## 给复现协议的结论

1. **torch 与 xformers 作为一对绑定版本冻结**，不能各自独立选版本。
2. 训练最小依赖集（README 未给出，此处为实测）：
   `torch` + `xformers`（版本互锁）、`mlflow`、`obonet`、`transformers`、`numpy`、`pandas`、`scikit-learn`、`scipy`、`tqdm`、`sgkit`、`prettytable`、`networkx`；
   完整复现图表另需 `matplotlib`、`seaborn`、`plotly`、`dash`、`pygraphviz`、`xgboost`。
3. 复现分支应自带锁定版本清单，并在 README 补上这四个漏掉的依赖。
4. **Stage A 需要 GPU**（或给 `:239` 打一行守卫）。这一点在排期时必须先确认硬件。

## 附：参考 checkpoint 的可读性

`samples/output_model.pt.*` 不是纯 `state_dict`，而是整体 pickle，内含 `argparse.Namespace`：

```
WeightsUnpickler error: Unsupported global: GLOBAL argparse.Namespace
```

两个含义：
- **好的一面**：checkpoint 里记录了完整的训练参数命名空间。若作者日后放出 UKB 训练好的权重，超参可以直接从 checkpoint 读出，不必依赖论文正文。
- **坏的一面**：加载它必须 `weights_only=False`，即执行 pickle 代码。复现协议里应写明这一条（来源可信仍应记录为一次显式的信任决定）。
