# Stage A 环境重建实测记录

**目的：** 验证"不用 `environment.yml`、改用自建最小环境"这条路能否跑通上游 `train_model.sh` 的合成数据链路。
**平台：** Linux x86_64，**无 GPU**（纯 CPU），Python 3.11.15，venv（非 conda）。
**上游：** `idekerlab/G2PT` @ `fc13abd`，只读浅克隆，不纳入本仓库版本管理。

> ⚠️ **本记录不能替代作者环境的一致性验证。** 这里装的是各依赖的当代版本，与作者环境（Python 3.8 / torch 2.4.1+cu124 / numpy 1.24 / pandas 2.0.3 / sgkit 0.7.0）差距很大。本记录只回答"代码在现代依赖下能不能跑起来"，**不回答"结果是否与作者一致"**。后者必须靠与 `samples/` 附带的参考输出比对，且应在贴近作者版本的环境里做。

## 尝试 1 —— 失败（xformers）

装了 `torch numpy pandas scikit-learn scipy tqdm networkx prettytable statsmodels matplotlib sgkit bed-reader` 后直接跑训练，立即失败：

```
File "src/model/hierarchical_transformer/hierarchical_transformer.py", line 4, in <module>
    import xformers.ops as xops
ModuleNotFoundError: No module named 'xformers'
```

**结论（→ 审计 B7）：** `xformers` 是训练路径的硬依赖，且 README 的安装说明里完全没有提到它——README 只给了 `conda env create -f environment.yml` 一条路，而这条路因 B5 大概率走不通。`xformers` 只出现在 `environment.yml` 的 pip 段（`xformers==0.0.28.post1`）。

连带后果：**xformers 每个版本硬绑定特定 torch 版本**，所以 torch 版本不是自由变量。首次安装拿到的 `torch 2.14.0+cu130` 必须让位给 xformers 指定的版本。

## 尝试 2 —— 进行中

改装作者 pip 段里的 `xformers==0.0.28.post1`（由它反向决定 torch 版本），再补 `mlflow`（`src/utils/trainer/snp2p_trainer.py:4` 的硬导入，同样不在 README 里），然后跑上游 `train_model.sh` 的等价命令（去掉 `--cuda 0`，本机无 GPU）。

结果待回填。

## 给复现协议的结论（当前）

1. **torch 与 xformers 必须作为一对绑定版本冻结**，不能各自独立选版本。
2. 训练最小依赖集（README 未完整给出，此处为实测）：
   `torch` + `xformers`（版本互锁）、`mlflow`、`numpy`、`pandas`、`scikit-learn`、`scipy`、`tqdm`、`sgkit`、`prettytable`、`networkx`；
   可视化另需 `matplotlib` + `pygraphviz`。
3. 复现分支应当自带一份锁定版本清单，并在 README 里写清 `xformers`/`mlflow` 这两个 README 漏掉的依赖。
