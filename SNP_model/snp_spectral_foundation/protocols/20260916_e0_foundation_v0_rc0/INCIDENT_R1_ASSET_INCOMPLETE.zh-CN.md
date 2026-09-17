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
| Fileset 命名 | `ukb_imp_chr{11..22}_v3.{bed,bim,fam,log}` |
| 染色体覆盖 | 仅 `chr11`–`chr22`；48 个文件 = 12 染色体 × 4 文件。**`chr1`–`chr10` 不在此目录** |
| Fileset 数 | 12，全部判为 `imputed_candidate` |
| 样本数 | 487,409（各染色体一致） |
| Variant 数 | `chr11` 4,628,348 至 `chr22` 1,255,683 |
| Build hint | 有判别力的染色体一致指向 `GRCh37`，无冲突 |
| 单染色体 `.bed` | 15.8 GB（chr21）至 564 GB（chr11） |
| `e0_primary_target_eligible` | 12 个 fileset 全为 `false` |
| 输出盘 | `/home/tyuxiao` 可用 70,943 GiB；`/data3` 仅余 3.0 TB |

## 🗂️ `genotype/` 全树（第二轮扫描）

`/data3/ukb_all/genotype/` 共七个实体子目录，合计约 11.1 TB。全部为 imputation v3 的 PLINK 转换产物，**没有任何 array、haplotype、kinship、sample-QC 或 mfi 文件**。

| 子目录 | 体积 | 染色体 |
| --- | --- | --- |
| `新建文件夹` | 3.4 TB | `chr11`–`chr22` |
| `data2_external_drive_copies` | 2.4 TB | `chr2`（`gsq_root`）、`chr3`、`chr4`（`onetouch_root`） |
| `data1_ukb_crosschr_fixed` | 1.6 TB | `chr3`、`chr5`、`chr11` |
| `external_new_chr8_9_10` | 1.6 TB | `chr8`、`chr9`（`gsq`）、`chr10`（`onetouch`） |
| `mnt_onetouch` | 1.3 TB | `chr6`、`chr7` |
| `mnt_gsq` | 841 GB | `chr1` |
| `$RECYCLE.BIN`、`System Volume Information` | 20 KB | Windows 残留，非数据 |

结论：`chr1`–`chr22` 全部齐备（`chr3`、`chr11` 等存在跨目录重复副本），所有 `.fam` 均为 12,671,752 bytes，即同一套 487,409 人。**覆盖范围不再是阻断；数据类型仍然是阻断。**

目录来源可追溯：`mnt_gsq`、`mnt_onetouch`、`data2_external_drive_copies` 指向外接硬盘，迁移记录在 `/data3/ukb_all/logs/`（`rsync_all_ukb_20260503.log`、`source_manifest_all_ukb_20260503.tsv`）。原盘上是否还有 array 或 mfi 资产未随迁，需查该 manifest。

## 🧾 `phenotype/` 中的原始交付物

`/data3/ukb_all/phenotype/mnt_gsq_xinjianwenjianjia/` 保存了 UKB 主数据集的 SAS 原件，合计约 380 GB，其中 `ukb_alldata.sas7bdat` 为 123 GB。`phenotype/derived/` 下是已按 `chr1` fam 顺序对齐的派生表，包括 `ukb_covariates_chr1fam_aligned.tsv`（212 MB）。

这批资产对 E0 有两处潜在价值，均待验证：

1. **分层变量**：若 covariates 或 SAS 主表含 genetic PCs、genotyping array、batch、sex，可直接作为 `R2` 的 `--strata-file`。
2. **亲缘信息的替代来源**：UKB field `22021`（genetic kinship to other participants）是每个个体的汇总标记，而非配对边。它不能重建 connected components，但可支持一条替代设计——只保留 `22021 == 0` 的无亲缘子集，使 family-disjoint 由构造保证。该替代方案需作为协议修正单独登记，并明确它把 `δ_S` 的解释从「新家系」收窄为「新的无亲缘个体」。

`ukb_event_labels_*`、`ukb5_lab_traits_*` 与 SAS 主表中的 trait 字段属于 outcome，E0 全程不读；上述检查只涉及列名，不读取任何数据行。

## 🧩 迁移脚本解释了缺口，而不是证明资产不存在

`/data3/ukb_all/logs/run_all_ukb_migration_20260503.sh` 的拷贝范围是写死的两个 glob：

```text
ukb_imp_chr*_v3.{bed,bim,fam,log}
*.sas7bdat
```

因此 array calls、haplotype、`ukb_rel_*.dat`、`ukb_sqc_*`、`ukb_mfi_*` **在构造上不可能被迁移过来**，与它们在原盘是否存在无关。`rsync` 日志中这些模式命中 0 次，只说明脚本没尝试过，不构成原盘没有的证据。

原始来源路径（来自 `source_manifest_all_ukb_20260503.tsv`）：

| 来源 | 内容 |
| --- | --- |
| `/data1/jingyixi/ukb_crosschr_fixed/` | `chr3`、`chr5`、`chr11` |
| `/data2/external_drive_copies/{gsq_root,onetouch_root}/` | `chr2`、`chr3`、`chr4` |
| `/mnt/gsq/` | `chr1`；`/mnt/gsq/新建文件夹/` 为 SAS 原件 |
| `/mnt/onetouch/` | `chr6`、`chr7` |

这些位置尚未被检查，是下一步搜索 array/kinship/mfi 的首选目标。

## 🔑 授权与转换链（来自 PLINK 日志）

`mnt_gsq/ukb_imp_chr1_v3.log` 给出完整的上游链条：

| 项目 | 值 |
| --- | --- |
| Application ID | **44430**（由 `--sample ukb44430_imp_chr1_v3_s487298.sample` 推得） |
| 上游输入 | `ukb_imp_chr1_v3.bgen`（BGEN v1.2）+ 官方 `.sample` |
| 转换工具 | PLINK v2.00a2.1 32-bit（2020-01-14，主机 `DESKTOP-UUS2K98`，`D:\UKB`） |
| chr1 variants | 7,402,791 |
| 样本 | 487,409（264,303 女 / 222,995 男 / 111 ambiguous） |

**必须记录的转换风险**：日志含

```text
Warning: No --bgen REF/ALT mode specified ('ref-first', 'ref-last', or 'ref-unknown').
```

UKB imputed BGEN 为 ref-first。未显式指定时，`.bim` 的 A1/A2 与官方 REF/ALT 的对应关系不确定。对 E0 内部任务而言，全文件统一的方向翻转只是把 genotype class `0` 与 `2` 一致地互换，masked-genotype CE 不受影响；但它使以下操作在未审计前不可进行：

- 与外部 reference panel（含 `M1b` 的 Beagle reference）比对
- 与官方 `mfi` AF 或任何外部 AF 交叉验证
- 任何跨 cohort 的 allele 方向断言

因此 `CONTRACT.json` 必须显式记录 allele 方向为「PLINK 转换自定义，未对齐官方 REF/ALT」，并在需要外部比对时先做方向审计。

## 📊 分层变量已齐备

`phenotype/derived/ukb_covariates_chr1fam_aligned.tsv` 共 53 列，按 `chr1` fam 顺序对齐：

| 列 | 字段 | 用途 |
| --- | --- | --- |
| `fid`、`iid`、`eid`、`n_eid` | 标识 | 与 genotype fam 对齐 |
| `n_31_0_0`、`n_22001_0_0` | 自报 sex、genetic sex | strata |
| `n_22000_0_0` | genotype measurement batch | strata（可区分 UKBiLEVE 与 Axiom array） |
| `n_22006_0_0` | genetic ethnic grouping | strata |
| `n_54_0_0` | assessment centre | strata |
| `n_22009_0_1` – `n_22009_0_40` | **genetic PC 1–40** | strata 与 `M1a` 协变量 |
| `n_34_0_0`、`n_52_0_0`、`n_21003_0_0`、`n_21022_0_0` | 出生年月、年龄 | 非 E0 所需 |

`R2` 的 `--strata-file` 输入因此已经具备。**该表不含 field `22021`**，亲缘信息仍需另找来源。

## 🔍 全盘搜索结果

在 `/data3` 与 `/home/tyuxiao` 的 `maxdepth 4` 范围内搜索 `*_cal_*`、`ukb_snp_*`、`*_hap_*`、`*rel*.dat`、`*.kin0`、`*sqc*`、`*mfi*`，**返回为空**。该深度已覆盖 `genotype` 各子目录中的文件。

`/data3/ukb_all/` 下存在四个同级目录：`genotype`、`phenotype`、`logs`、`runs`；`genotype` 本身还有多个尚未列出的子目录，`chr1`–`chr10` 与其余资产可能在其中。全部资产属主为 `jingyixi`（组 `tpc`），即本集群的数据管理方。

## 🚫 四条阻断

| 阻断 | 性质 | 影响 |
| --- | --- | --- |
| 全部 22 条染色体均只有 imputed fileset | **科学性阻断** | E0 primary target 必须是观测分型 |
| 无可解析 kinship/relatedness 文件 | **硬阻断** | family connected components 无法冻结，`R2` family 轴停摆。应用号 44430 已知，官方 `ukbgene rel` 可直接下载该文件 |
| 无 genetic map | 工程阻断 | `1 cM` guard 规则无法执行 |
| ~~未探测到 CUDA-enabled PyTorch~~ | **已解除（2026-09-17）** | 系误报：inventory 由系统 `python3` 3.10.12 运行。在 `/home/tyuxiao/miniforge3/envs/grassbench`（Python 3.11.15）中实测 `torch 2.11.0+cu128`、`torch.cuda.is_available() == True`。micro-batch 与 precision 的绑定不再被环境阻断，仍需实测吞吐与显存 |

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

1. 列出 `/data3/ukb_all/genotype` 其余子目录与 `phenotype`、`runs`，定位 `chr1`–`chr10` 及可能的 array/kinship 资产
2. 向数据管理方 `jingyixi` 申请 `ukb_cal_*`、`ukb_rel_*.dat`、`ukb_sqc_v2.txt`，以及退路所需的 `ukb_mfi_chr*_v3.txt`
3. ~~在含 torch 的 miniforge 环境复测 CUDA~~ —— 已完成，见上表
4. 传入已冻结哈希的 `geneticMap-GRCh37` 镜像，清除 genetic map 阻断
5. 取得 kinship 后重跑 R1，`PASS` 后才进入 R2 family 轴

在第 1–2 步有结果前，不生成任何 manifest，不签发 `RUN_AUTHORIZED`。

即便走路径 B，当前可见资产也只覆盖 `chr11`–`chr22`；primary claim 的染色体范围需相应限定或补齐。

## 🔗 相关

- [E0 协议](PROTOCOL.zh-CN.md)
- [R0.5 状态报告](STATUS_20260916.md)
- [R1 执行包](r1_ukb_binding/README.zh-CN.md)
