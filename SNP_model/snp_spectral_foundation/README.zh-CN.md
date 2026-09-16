# SNP Spectral Foundation：可执行候选骨架

> ⚠️ **当前主线（2026-09-15）：** 本目录的模型代码保留为候选工程资产，但当前研究不再先做新的 foundation-model attribution audit，也不启动新架构训练。权威路线是 0526 Stage 2 phenotype bridge，见 `protocols/20260915_0526_stage2_bridge_rc1/PROTOCOL.zh-CN.md`；该协议当前已冻结但因真实 phenotype 资产未绑定而禁止正式运行。

## 结论先行

本目录把目标流水线实现成了可训练的 PyTorch 原型，但把 Grassmann 当作**候选分支**而非既定真理。原因是本项目既有真实 HapMap3 gate 已表明：当前 Grassmann memory 在相同 16-byte 预算下被 aligned PCA 在 8/8 个域稳定支配；主要损失是有符号、位置对齐的系数。因此默认模型同时保存 spectral structure 与独立 residual sidecar，并保留 `aligned_control` 做强对照。

## 对应关系

```text
raw dosage {-1,0,1,2}
  → genotype + global block ID + within-block position + train-only AF/MAF embedding
  → shared local Transformer + masked-genotype head
  → per-person/per-block randomized low-rank spectral decomposition
  → Grassmann projector probes + eigenvalues + MAF summary
  ⊕ signed position-aligned dosage residual sidecar
  → block memory
  → latent-slot global mixer（对 block 数近似线性，而非 N² dense attention）
  → trait-ID-conditioned block reweighting
  → trait-effect representation
  → learned spectral/residual fusion
  → phenotype prediction
```

这里的 “phenotype-conditioned” 严格实现为 **trait identity conditioning**，不读取受试者的真实结局。若 effect representation 用于解释或关联分析，仍必须在数据流程中 cross-fit。

## 代码入口

- `src/snp_spectral_foundation/model.py`：完整模型、谱分解、全局 mixer 与融合。
- `src/snp_spectral_foundation/training.py`：mask、预训练 loss、marginal baseline gate。
- `scripts/smoke_train.py`：合成 LD 数据上的端到端冒烟运行。
- `tests/test_model.py`：shape、反向传播与 `U → UQ` Grassmann 不变性。
- `EXPERIMENT_CONTRACT.zh-CN.md`：逐层 gate、干预和消融合同。

## 运行

在本目录执行：

```powershell
$env:PYTHONPATH = "$PWD\src"
python -m unittest discover -s tests -v
python scripts\smoke_train.py
```

冒烟数据只证明代码链路能运行，不构成模型优越性或生物学证据。真实数据接入时必须从训练个体计算定向 ALT allele frequency、MAF 与 LD block，并先过 contextual pretraining gate。注意 dosage residual 的期望值是 `2 × ALT AF`，不能误用丢失等位基因方向的 MAF。

`block_ids` 是固定 SNP panel 中的全局 block 编号。分块/流式读取时必须传入真实编号，不能让每个 chunk 都从 0 重新编号；不同 genome build 或 panel 也不能直接复用同一编号表。

## 当前没有伪装完成的部分

- 没有在本目录复制或下载 SNPBag/GenoBERT checkpoint；公开模型审计应作为独立 benchmark adapter 接入。
- 没有把 dosage 当 phased haplotype；phase attribution 需要新的输入契约。
- 没有宣称全基因组可训练；当前实现的 latent-slot global mixer 提供线性扩展接口，但仍需真实吞吐与显存 profiling。
- 谱分解使用固定 sketch 的 randomized subspace iteration，避免为每个 person-block 构造完整 `d_model × d_model` covariance；它是可微近似，精度/吞吐仍需按 rank 与 power iteration 做 profile。
- 没有用当前已打开的 chr22 test 调参；正式比较应转向未使用染色体或外部 cohort。
