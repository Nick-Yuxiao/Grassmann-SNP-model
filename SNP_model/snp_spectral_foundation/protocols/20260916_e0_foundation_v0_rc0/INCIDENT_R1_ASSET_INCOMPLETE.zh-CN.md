# R1 事件报告：首个 UKB 目录只含 imputed 资产

_2026-09-16 · 里程碑 `R1`（未完成）· 状态 `ASSET_INCOMPLETE` · 事件报告，不是定时状态报告_

---

## 📋 结论

服务器 `c20-5090` 上第一个被授权只读扫描的 UKB 目录已完成 inventory。该目录**不能支撑 E0 正式运行**：它只含 imputation v3 的 hard-call，且同一挂载点下没有 kinship 与 genetic map。R1 因此停在 `ASSET_INCOMPLETE`，`R2` 仍然不可启动。

这不是 E0 协议的失败，也不改变任何 estimand；它是一条资产可用性事实。

## 🧬 已扫描资产

| 项目 | 值 |
| --- | --- |
| 主机 | `c20-5090.cbicrcluster.tianmouc.com` · Linux 5.15.0 · glibc 2.35 |
| 扫描根 | `/data3/ukb_all/genotype/新建文件夹` |
| 目录总量 | 3.4 TB |
| Fileset 命名 | `ukb_imp_chr{11..22}_v3.{bed,bim,fam,log}`（首屏 50 行所见；完整染色体覆盖待确认） |
| 单染色体 `.bed` | 15.8 GB（chr21）至 564 GB（chr11） |
| `.bim` 量级 | chr22 为 39 MB，约 1.26×10⁶ variants，符合 imputed 密度 |
| `.fam` | 所有染色体均为 12,671,752 bytes，指向同一套约 4.87×10⁵ 样本 |
| 分类结果 | `imputed_candidate`，`e0_primary_target_eligible = false` |
| 输出盘 | `/home/tyuxiao` 可用 70,943 GiB；`/data3` 仅余 3.0 TB |

## 🚫 四条阻断

| 阻断 | 性质 | 影响 |
| --- | --- | --- |
| 无 array-scale 直接分型 fileset | **科学性阻断** | E0 primary target 必须是观测分型 |
| 无可解析 kinship/relatedness 文件 | **硬阻断** | family connected components 无法冻结，`R2` family 轴停摆 |
| 无 genetic map | 工程阻断 | `1 cM` guard 规则无法执行 |
| 未探测到 CUDA-enabled PyTorch | 疑似误报 | inventory 由系统 `python3` 3.10.12 运行，非 miniforge 环境；需在含 torch 的环境复测 |

## ⚠️ 为什么 imputed 不能做 E0 primary target

E0 的 estimand 是

\[
p(G_{masked}\mid G_{observed\ context}, position, AF),
\]

其中 `G` 必须是观测到的 diploid 分型。imputation v3 的 hard-call 是 haplotype-HMM 后验的众数，不是观测量。若对其 masking：

1. **循环比较**：`M1b` 冻结为 Beagle 5.5 haplotype-copying/HMM。用 HMM 的输出作为 target，E 与 `B*` 实际比较的是谁更能复现 IMPUTE4 的后验，而不是谁更能建模真实 genotype 分布。`δ_S`、`δ_X` 因此不再度量 neural headroom。
2. **信息泄漏**：imputation 在全 panel 上完成，masked 位点的取值已由其邻域 LD 通过 reference panel 决定；这与 E0-X「未见区段迁移」的可识别性直接冲突。
3. **guard 失效**：held-out run 的 guard 只能隔断 panel 内的 LD，无法隔断 imputation 阶段已经跨过边界的信息传递。

因此 inventory 将此类 fileset 标为不合格，而不是降级使用。

## ▶️ 解除阻断的两条路径

### 路径 A（首选）：取得 array 直接分型资产

| 用途 | UKB 标准文件 |
| --- | --- |
| E0 primary panel | `ukb_cal_chrN_v2.bed`、`ukb_snp_chrN_v2.bim`、`ukb*_cal_*.fam` |
| family 轴（必需） | `ukb_rel_a<app>_s<n>.dat` |
| strata 平衡 | `ukb_sqc_v2.txt` |
| `M1b` reference | `ukb_hap_chrN_v2.bgen` |

### 路径 B（退路）：从 imputed 容器中提取 array-genotyped 位点

UKB imputation v3 包含直接分型位点，并在 `ukb_mfi_chrN_v3.txt` 中以 `INFO = 1` 标记。可只保留该子集构成 panel。采用此路径必须：

- 在 `CONTRACT.json` 中把 estimand 范围显式限制为「array-genotyped sites extracted from an imputed container」
- 记录这是退路而非首选，并说明 hard-call 与原始 array call 可能存在的不一致
- 单独审计 REF/ALT 方向与 strand

路径 B 不能解除 kinship 缺口。**在 kinship 到位前，`R2` 无论走哪条路径都不可启动。**

### genetic map

genetic map 不是 UKB 受限资产。优先复用开发阶段已冻结哈希的镜像 `geneticMap-GRCh37`（commit `42e69e59f1d3a31725378e61d8f8e82a51e88d63`），以保持与 `DEVELOPMENT_REPORT.zh-CN.md` 的哈希链一致。

## 📅 下一步

1. 确认扫描根的完整染色体覆盖与是否存在非 `imp` 文件
2. 在 `/data3` 与用户主目录范围内搜索 `*_cal_*`、`ukb_snp_*`、`*_hap_*`、`*rel*.dat`、`*.kin0`、`*sqc*`、`*mfi*`
3. 在含 torch 的 miniforge 环境复测 CUDA，清除第四条阻断
4. 取得 kinship 后重跑 R1，并在 `PASS` 后才进入 R2 family 轴

在上述第 1–2 步有结果前，不扩大扫描范围，不生成任何 manifest，不签发 `RUN_AUTHORIZED`。

## 🔗 相关

- [E0 协议](PROTOCOL.zh-CN.md)
- [R0.5 状态报告](STATUS_20260916.md)
- [R1 执行包](r1_ukb_binding/README.zh-CN.md)
