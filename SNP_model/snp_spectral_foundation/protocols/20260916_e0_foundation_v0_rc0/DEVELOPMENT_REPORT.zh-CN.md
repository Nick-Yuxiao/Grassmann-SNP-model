# E0 双轴 manifest 开发验证报告

_2026-09-16 · `R0.5` · `DEVELOPMENT_ONLY_NON_EVIDENCE`_

---

## 📋 结论

E0 的 family × genomic-block 双轴 manifest 生成与独立泄漏验证链已经在公开数据上跑通，但正式 UKB 运行仍未获授权。当前结果只证明工程结构可执行，不证明 encoder 有效，也不构成 UKB 上的 family-disjoint 或 unseen-locus 证据。

本轮还发现并修正了一个实质性设计问题：随机散布 held-out blocks 后逐块增加 `1 cM` guard 会使 buffer 膨胀。正式协议改为按 strata 抽取多个连续 held-out runs，仅在 run 外边界设置 guard；run 内 LD blocks 继续作为评价和 two-way bootstrap 的 genomic 聚类单位。

## 🧬 开发数据与边界

| 项目 | 冻结值 |
| --- | --- |
| 数据 | 1000 Genomes Phase 3 chr22 VCF |
| Build / contig | GRCh37 / `22` |
| 样本数 | 2,504 |
| VCF SHA-256 | `5d52145d6e3be132dc19ff525af2ae2e792511e60c49c148cb8eee669f114c2` |
| Panel SHA-256 | `b4023dc6ee2d62ee89c8d4d347db4d348e65518d66d346574cdae7a4bbd76858` |
| Relationship file SHA-256 | `359d59d1263c208fe8876817cc1ddc98aa3096c2c2212d3e13fdfb92b25ece3d` |
| Genetic map SHA-256 | `057f35f9e53bcaca83b04fb9a4c76948f19b0f5cef0d88b548314b9670245d50` |
| Coordinate contract | VCF `POS` 1-based；manifest `start0/end0` 0-based half-open；`start1/end1` 1-based inclusive |

开发 VCF 的前 100 个 records 共检查 250,400 个 callable GT，phased fraction 为 `1.0`。该抽查只验证开发输入能够支持 haplotype-aware 代码路径，不代表全文件已经完成 phase、missingness 或 allele QC。

公开 related-individual 文件中的个体均不在当前 2,504-sample VCF panel，因此实际 relationship edge 为 `0`，所有 family components 均为 singleton。由此，公开数据没有验证多成员 connected component 的真实分布；这部分只由合成单元测试覆盖，必须在 UKB kinship 上重新运行。

## 🧪 生成与验证结果

| 检查 | 结果 |
| --- | --- |
| Family split | train 1,753；validation 376；test 375 |
| Family components | 2,504；跨 split component 0；最大 component size 1 |
| Development blocks | train 93；validation 22；test 22 |
| Buffers | validation 5；test 5 |
| Block definition | 0.5 cM engineering blocks；1 cM guards |
| Coordinate conversion errors | 0 |
| Cross-split core/guard overlaps | 0 |
| Independent validator | `PASS` |
| Protocol automated tests | 11/11 passed |

Family manifest SHA-256 为 `1731a4aad4abab33b3ddd3b56f30b7bd10f96647a09e47fe3f70064af542d556`；block manifest SHA-256 为 `1c5d8755569c03d0a0d91704fcbfc6b95180864302b9aff90d71d65417348231`。

验证器不输出样本标识，只输出 aggregate counts、错误数与 manifest hashes。故意构造的 family crossing、guard overlap、非法坐标与 block overlap 均会使测试失败。

### M0 validation-only dry run

在不读取 development test 的前提下，使用 1,753 个 train samples 估计 M0，在 376 个 validation samples × 2 个 validation blocks 上生成一个 mask seed：

| 项目 | 结果 |
| --- | --- |
| 合格 biallelic SNPs | 2,829 |
| Masked targets | 212,589 |
| Family × block × seed cells | 752 |
| M0a empirical CE | 0.114902 |
| M0b HWE/AF CE | 0.126097 |
| M0a accuracy | 0.956756 |
| M0b accuracy | 0.950225 |
| Independent validation | `PASS` |

这些值只验证 evaluator、train-only probability fitting 和 mask identity；它们不是 encoder 性能，也不用于正式 E0 判决。独立验证器确认输出哈希、split labels、mask key 唯一性、core 坐标、概率归一化和 cell target 总数一致。大文件位于 `F:\Yuxiao Tan\Document\Grassmann_model_external\e0_runs\m0_validation_2blocks_rc0`。

M0 artifact hashes：mask `2ce63f4aca1dd785a4b74a604d5c52554189335f312646586bbbd14d121d5f25`；variant stats `5f69a19da9c5437f892d29e2a6762c9020d85ba8cff652536a9241034a903304`；cell metrics `01d8fd56ad9a150d8090067ba862edb92735f1cffb1dacc0390b5849fef15152`。

### E1 train-only real-VCF smoke

E1 primary 的真实数据最短链已通过：从公开 VCF 读取一个 train block，以全部 1,753 个 train samples 估计 ALT AF，取 4 个 train samples × 256 SNP，完成两次 masked-CE forward/backward/update。模型为冻结 E1 结构，828,291 parameters，`relative_continuous` 位置编码且无 local block embedding；loss 与 gradient 均有限。运行未读取 validation/test、未保存 checkpoint，输出位于 `F:\Yuxiao Tan\Document\Grassmann_model_external\e0_runs\e1_realdata_smoke_rc0`。

两步 loss 从 `1.1045` 到 `0.3030` 只说明微型 batch 可以反向传播和过拟合，不能解释为泛化、收敛或 neural headroom。

## 🚫 为什么尚不能启动正式 E0

- 尚无本机可访问的 UKB genotype、kinship、授权和正式 genetic-map 绑定。
- 开发 blocks 是固定 0.5 cM 工程单元，不是由 train families 或冻结 reference 构造的 LD blocks。
- 开发 guard 尚未执行正式的 `max(receptive field, adjacent LD block, 1 cM)` 与跨边界 `r²` 审计。
- 公开 panel 没有多成员 family component，不能替代 UKB 亲缘泄漏验证。
- M1b haplotype-copying/HMM 的软件、版本、reference 和 benchmark subset 尚未冻结。
- 计算资源与正式窗口吞吐尚未测算。

所以当前允许推进 manifest/loader、baseline adapter 和小规模 dry run；不允许读取 test outcome、训练正式 encoder 或声称 Foundation-v0。

## ▶️ UKB 到位后的最短路径

1. 只读 inventory：绑定 build、phasing、样本、variant、kinship、授权和文件哈希计划。
2. 以 kinship connected component 冻结 family manifest，并确认跨 split component 为 0。
3. 只用 train families 或冻结 reference 构造 LD blocks；按 chromosome、位置、block size、MAF、LD density 分层抽取连续 held-out runs。
4. 依据 encoder receptive field、相邻 LD block 和 `1 cM` 最大值构造 guards，再用 train families 完成跨边界 LD 审计。
5. 独立验证全部 manifests，冻结 M0/M1/HMM 与最多三个 encoder recipes，之后才签发 `RUN_AUTHORIZED`。

## 🔗 资产

- Development run：[`development/1000g_chr22_rc3_cm`](development/1000g_chr22_rc3_cm/)
- Manifest generator：[`scripts/build_dev_manifests.py`](scripts/build_dev_manifests.py)
- Independent validator：[`scripts/validate_manifests.py`](scripts/validate_manifests.py)
- Tests：[`tests`](tests/)
- Encoder implementation audit：[`IMPLEMENTATION_AUDIT.zh-CN.md`](IMPLEMENTATION_AUDIT.zh-CN.md)
- Frozen recipe template：[`MODEL_RECIPES.template.json`](MODEL_RECIPES.template.json)
- Frozen baseline template：[`BASELINES.template.json`](BASELINES.template.json)
- 5W1H execution sheet：[`E0_5W1H.zh-CN.md`](E0_5W1H.zh-CN.md)
- SNPBag 官方代码镜像：`F:\Yuxiao Tan\Document\Grassmann_model_external\SNPBag`，commit `984b8e2af7198005142eeea6f94e848ea2be0d28`
- GRCh37 genetic map 镜像：`F:\Yuxiao Tan\Document\Grassmann_model_external\geneticMap-GRCh37`，commit `42e69e59f1d3a31725378e61d8f8e82a51e88d63`
- Beagle 5.5：`F:\Yuxiao Tan\Document\Grassmann_model_external\beagle_5.5_27Feb25.75f`；jar SHA-256 `7319f4af9638be05c18dcc1bfb8fb41a58a09293507ebf0d54617d0e40df5a70`

SNPBag 官方仓库当前提供代码、demo 与示例 embedding，但没有可直接下载的 full-genome pretrained weights；因此它不能替代本项目从头训练的 E0，也暂不能执行公开 checkpoint 的 `C vs B` phenotype benchmark。[SNPBag repository](https://github.com/augix/SNPBag)；[公开权重请求 issue](https://github.com/augix/SNPBag/issues/3)。
