# 功能注释来源（F arm）

## 0. 服务器现状

之前在 `/data3` 和 home 下的 `find` **全空**：服务器上目前没有任何功能注释文件。
F 必须先取公共数据。这些文件不是 UKB 数据，没有 UKB 数据协议问题，
可以在有外网的机器上下载再拷进去。

## 1. 三条硬约束

1. **基因组版本必须是 GRCh37 / hg19。** panel 是 GRCh37，本包**不做 lift-over**。
   GRCh38 的注释会 join 到 0 个位点，`c2` 会直接停下来。这是预期行为。
2. **必须 phenotype-free。** 不能用由性状关联导出的任何东西
   （GWAS Catalog、同性状 eQTL、PRS/PGS 权重、TWAS、fine-mapping 后验）。
   `c0` 对名字或路径里含 `gwas / catalog / prs / pgs / eqtl / sqtl / pqtl / twas /
   finemap / heritab / sumstat / ldsc / magma` 的 `--local` track 直接拒绝。
3. **只用位置 join。** A1/A2 来自一次没有显式 REF/ALT 模式的 BGEN 转换，
   allele-aware join 不可靠。

## 2. 内置 catalogue

`c0_prepare_annotations.py` 里写死了下面几条，每条都带"为什么它是 phenotype-free"
和"它的局限"。默认只取前三条（`k562_chromhmm_broad`、`phastcons100way`、`tss_proximity`）。

| track | 内容 | 红系相关性 | 局限 |
|-------|------|-----------|------|
| `k562_chromhmm_broad` | Broad ChromHMM 染色质状态，K562 | **是**（K562 是红白血病系，红系染色质的标准替身） | 是癌系，不是原代红系 |
| `roadmap_e123_k562_chromhmm` | Roadmap 15-state core marks，E123 = K562 | **是** | 同上；与上一条重叠很大 |
| `dnase_clusters` | ENCODE 跨细胞系 DNase 可及性 cluster | 否 | 泛组织，作为"通用调控"的对照 |
| `phastcons100way` | 100 物种保守元件 | 否 | 强但通用，算是功能轴上的近似 null |
| `tss_proximity` | 到最近 TSS 的距离，`exp(-d/10kb)` | 否 | 只是基因结构；用来看 F 里有多少其实只是基因密度 |

**为什么 ChromHMM 要按状态拆开**：`split_by_name=True` 会把 segmentation 按状态名
拆成一条条二值 track（strong enhancer、active promoter、heterochromatin……）。
合成一条会把增强子和异染色质加在一起，互相抵消。拆开之后，F 的每一列有明确含义，
`F_columns.tsv` 里能读出是哪个状态、哪条染色体。

## 3. URL 与下载

`c0` 用 `curl --fail` 下载。**这些 URL 我无法在本环境验证**（出站代理挡掉了
`hgdownload.soe.ucsc.edu` 和 `egg2.wustl.edu`），所以下载前先自己确认一次：

```bash
curl -sS -I "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/wgEncodeBroadHmmK562HMM.txt.gz" | head -1
curl -sS -I "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/phastConsElements100way.txt.gz" | head -1
curl -sS -I "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz" | head -1
```

拿到 `200` 再跑 `c0`。如果某条 404 了，去 UCSC 的
`https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/` 目录页找当前文件名，
然后用 `--local <track>=<path>` 传进去。`c0` 不会因为某条 track 失败而把整批作废，
但会把每条的状态写进 `ANNOTATION_SUMMARY.json`。

**UCSC `database/*.txt.gz` 的格式**：这类表第一列是 `bin`，不是染色体。`c0` 的
`ucsc_bed_bin` 格式会跳过它。如果你手工下了 `.bed` 而不是 `.txt.gz`，用
`--local-format bed`。

## 4. 服务器上的执行顺序

有外网的机器上：

```bash
python3 scripts/c0_prepare_annotations.py \
  --panel-dir "$PANEL" \
  --out-dir "$ANN" \
  --track k562_chromhmm_broad \
  --track phastcons100way \
  --track tss_proximity
```

没有外网的机器上：先把文件拷进去，然后

```bash
python3 scripts/c0_prepare_annotations.py \
  --panel-dir "$PANEL" --out-dir "$ANN" --offline \
  --local k562_chromhmm_broad=/path/wgEncodeBroadHmmK562HMM.txt.gz \
  --local phastcons100way=/path/phastConsElements100way.txt.gz \
  --local tss_proximity=/path/refGene.txt.gz
```

输出：`$ANN/annotation_manifest.tsv`（直接喂给 `c2`）、每条 track 的 `.bed` /
`.points.tsv`、`$ANN/ANNOTATION_SUMMARY.json`。

## 5. 先看 coverage，再往下跑

`ANNOTATION_SUMMARY.json` 里每条 track 有 `coverage` = 被注释到的 panel variant 占比。
panel 是 5 万 bp 间隔 thin 过的 2.5 万个位点，所以：

- ChromHMM 的**开放/活跃**状态覆盖通常是个位数百分比 —— 正常；
- ChromHMM 的**异染色质/低信号**状态覆盖会很高 —— 也正常；
- 某条 track coverage **= 0** —— 几乎一定是基因组版本搞错了，先查这个；
- 所有 track 都 = 0 —— `c0` 会直接报错并指出版本问题。

覆盖低于 `--min-coverage`（默认 0.2%）的 track 会被从 manifest 里剔除，
但仍然出现在 `reports` 里，避免"悄悄少了一条"。

## 6. 之后可以加、但现在不加的

- BLUEPRINT 的原代 erythroblast ChIP-seq（比 K562 更贴近真实红系，但要按样本处理，
  且需要确认 hg19）；
- K562 的 GATA1 / TAL1 / KLF1 ChIP-seq peak（ENCODE 文件名不稳定，建议手工下载后
  用 `--local` 传入）；
- CADD / phyloP 这类连续打分（bigWig，需要 `bigWigAverageOverBed`，比现在这套重）。

加 track 的代价是 F 变宽、`delta P | F` 的解释变复杂。**先用默认三条把管线跑通，
拿到量级，再决定要不要加。**
