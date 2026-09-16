# E0 local encoder 实现审计

_2026-09-16 · 适用范围：E0-X 未见 genomic region 迁移_

---

## 📋 判定

现有历史 `LocalGenotypeEncoder` 的默认结构不能原样作为 E0 primary：它同时使用 block 内 position lookup 和 learned global block embedding。后者会给训练过的 block 一个专属参数，并使未见 block 没有等价表示，因此 E0-X 的结果无法区分“共享 local rule”与“坐标/ID 记忆”。

代码已增加向后兼容的两种模式：

- `legacy_lookup`：保留旧 checkpoint 所需的 position/block lookup，只用于复现旧结果。
- `relative_continuous`：以 block/window 内相对位置、相邻 marker bp gap 和 window span 编码位置；E0 同时设置 `use_local_block_embedding=false`。

E0 模式若收到 `block_ids` 会显式报错。自动测试还确认：给所有 marker 坐标加同一个常数后，hidden state 与 logits 保持一致。这两个检查防止 E0 primary 静默退回坐标 lookup。

## ⚙️ 三个冻结候选配方

| Recipe | 结构 | 回答的问题 |
| --- | --- | --- |
| `E1_PRIMARY_LOCAL_TRANSFORMER` | 256 SNP，4 层，`d=128`，8 heads，无 ancestry conditioning | 默认跨区段共享 local encoder |
| `E2_ANCESTRY_CONDITIONED` | E1 + 仅由 train families 拟合的 PC-FiLM | 显式 population conditioning 是否改善跨 ancestry transfer |
| `E3_CAPACITY_SENSITIVITY` | 256 SNP，6 层，`d=192`，8 heads | E1 失败是否只是明显容量不足 |

三者均从随机初始化训练，保留 train-only ALT AF、MAF/LD/haplotype-related context 等 genotype structure；均禁止 global block ID、variant lookup 和 block-specific parameters。`E2` 的 PC 只用于 conditioning，不将 ancestry 从基本表示中“抹掉”。

同一配方最多 3 个 training seeds。只用 validation family × block macro CE 选一个 recipe/checkpoint，并用 one-standard-error rule 优先更小模型；test 不参与选择。三配方均未达到 validation 门槛后，本轮结束，不能继续追加 recipe 救援。

## 🖥️ 当前资源约束

只读盘点显示主机可见 20 个逻辑 CPU 和一个约 6 GB 显存的 RTX 3060 Laptop GPU 候选，但 Windows 与现有 WSL `epinformer_repro` 环境中的 PyTorch 均为 CPU-only，CUDA runtime 尚不可用。C 盘工作文件系统只余约 0.46 GB，而 F 盘约余 189 GB。

因此：

- protocol、代码和小日志保留在仓库；VCF、tensor cache、checkpoint 与 run outputs 放 F 盘。
- 当前只允许 CPU 上的 unit tests、manifest、M0 小规模 dry run。
- 正式训练前必须在实际执行环境验证 CUDA runtime，并以 development-only 窗口冻结 micro-batch、gradient accumulation 和 precision。
- 不通过显存压力测试推测 batch size；有效每 update masked-token 数固定，micro-batch 只负责适配显存。

## 🚫 尚未实现的部分

`E2` 的 train-PC FiLM adapter、正式 streaming VCF/PGEN loader、窗口中心 loss mask、M1a 与 M1b adapter 尚未实现。它们的缺失阻止 `RUN_AUTHORIZED`，但不改变已经冻结的 estimand、split 和 E1 输入契约。

E1 的 train-only real-VCF smoke 已完成：真实 VCF → train-only ALT AF → deterministic mask → integer relative positions → masked CE → backward/optimizer update 全链通过。该运行只含 4 个 train samples、256 SNP 和 2 updates，不构成性能结果。
