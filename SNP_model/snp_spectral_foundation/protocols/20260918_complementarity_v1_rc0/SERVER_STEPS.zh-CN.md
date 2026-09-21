# Concept Gate 服务器执行步骤

_路径按 `c20-5090` 写。全程**不打开** `tq_test` 与 `bridge_holdout`，**不导出**任何个体级数据。_

---

## 0. 解包与环境

```bash
cd /home/tyuxiao/Grassmann_model/SNP_model/incoming
tar -xzf complementarity_v1.tar.gz && cd 20260918_complementarity_v1_rc0
sha256sum -c MANIFEST.sha256

export OUT=/home/tyuxiao/Grassmann_model/e0_runs
export PANEL=$OUT/taskb_panel_rc0
export COV=/data3/ukb_all/phenotype/derived/ukb_covariates_chr1fam_aligned.tsv
export TRAITS=$OUT/taskb_extract/traits.tsv
export SPLIT=$OUT/taskb_split_rc0/tq_split_manifest.tsv
export ANN=$OUT/annotations_hg19
export GATE=$OUT/concept_gate_rc0

python3 -m unittest discover -s tests      # 20 个全过
```

只需要 numpy，和跑 Task B 时同一个解释器。**确认用的是装了 numpy 的那个**：

```bash
python3 -c "import numpy, sys; print(sys.executable, numpy.__version__)"
```

---

## 1. 取功能注释（一次性）

先确认 URL 还活着：

```bash
for u in \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/wgEncodeBroadHmmK562HMM.txt.gz \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/phastConsElements100way.txt.gz \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz ; do
  printf '%s -> ' "$u"; curl -sS -I --max-time 30 "$u" | head -1
done
```

都是 `200 OK` 就直接跑：

```bash
python3 scripts/c0_prepare_annotations.py \
  --panel-dir "$PANEL" --out-dir "$ANN" \
  --track k562_chromhmm_broad --track phastcons100way --track tss_proximity
```

服务器没有外网 → 在别处下好拷进来，再跑 `--offline` 那一套（见
`ANNOTATION_SOURCES.zh-CN.md` 第 4 节）。

**跑完先看 coverage**，全 0 就是基因组版本错了，不要往下走：

```bash
python3 -c "
import json; s=json.load(open('$ANN/ANNOTATION_SUMMARY.json'))
print(s['panel_variants'],'panel variants')
[print(f\"{r['track']:45s} {r.get('coverage')}\") for r in s['reports']]"
```

---

## 2. 构造 P（人群结构载体）

`--covariates` 在这里**只用来量 P 与 PC1–40 的重叠**，不参与拟合。

```bash
PCS=""; for i in $(seq 1 40); do PCS="$PCS --covariate-column n_22009_0_$i"; done

python3 scripts/c1_build_population_features.py \
  --panel-dir "$PANEL" \
  --split-manifest "$SPLIT" \
  --covariates "$COV" --covariate-id-column eid \
  --global-pc-prefix n_22009_0_ --global-pc-count 40 \
  --block-size 100 --pcs-per-block 3 \
  --out-dir "$OUT/features_P"
```

2.5 万位点、block 100、每块 3 个 → 约 750 列。

**必须看这一行**：

```bash
python3 -c "
import json; s=json.load(open('$OUT/features_P/P_SUMMARY.json'))
print(json.dumps(s['overlap_with_global_pcs'], indent=2, ensure_ascii=False))"
```

`top_canonical_correlations` 接近 1 → P 基本是 A 已有的 PC 的重复，`delta P | F`
没有解释力，把这件事写进结果，不要假装它是个发现。

---

## 3. 构造 F（功能载体）

```bash
python3 scripts/c2_build_functional_features.py \
  --panel-dir "$PANEL" \
  --split-manifest "$SPLIT" \
  --annotation-manifest "$ANN/annotation_manifest.tsv" \
  --group-by chromosome \
  --out-dir "$OUT/features_F"

python3 -c "
import json; s=json.load(open('$OUT/features_F/F_SUMMARY.json'))
print(s['shape']); print(s['tracks_with_zero_coverage'])"
```

---

## 4. 跑 Concept Gate

**先用 MCV 冒烟**（MCV 在本步是 development trait，不是 confirmatory primary）：

```bash
python3 scripts/c3_complementarity_gate.py \
  --panel-dir "$PANEL" \
  --covariates "$COV" --covariate-id-column eid \
  $PCS --covariate-column n_21003_0_0 --square-column n_21003_0_0 \
  --categorical-column n_22001_0_0 --categorical-column n_54_0_0 \
  --phenotype "$TRAITS" --phenotype-id-column n_eid \
  --trait-column n_30040_0_0 \
  --population-features "$OUT/features_P/P_population.npy" \
  --functional-features "$OUT/features_F/F_functional.npy" \
  --tq-results "$OUT/taskb_tq_mcv/TQ_RESULTS.json" \
  --split-manifest "$SPLIT" \
  --eval-split validation \
  --inner-folds 5 --bootstrap 2000 \
  --out-dir "$GATE/n_30040_0_0"
```

阈值从 `TQ_RESULTS.json` 的 `selected_p_threshold` 读入，**不重新挑**。
`--tq-results` 的 trait 与 `--trait-column` 不一致会直接报错。

跑通之后，对其余 trait 换 `--trait-column` / `--tq-results` / `--out-dir` 重复。
串行跑，PF 的设计矩阵不小。

内存不够会看到：

```
The PF design would need about NN.N GB (... rows x ... columns).
```

→ 回到第 2 步，把 `--block-size` 从 100 调到 200 或 400，P 的列数按比例减半、减到四分之一。

---

## 5. 回传（只回传这些）

```bash
cat "$GATE"/*/arms.csv
cat "$GATE"/*/comparisons.csv
python3 -c "
import json,glob
for p in sorted(glob.glob('$GATE/*/GATE_RESULTS.json')):
    r=json.load(open(p))
    print(r['trait_column'], r['verdict'], r['status'])
    print('  frozen p:', r['frozen_threshold']['value'],
          'variants:', r['frozen_threshold']['variants_selected'])
    print('  participants:', r['participants']['train'], 'train /',
          r['participants']['eval'], 'eval')"
```

再加上第 1–3 步的三个 summary：

```bash
cat "$ANN/ANNOTATION_SUMMARY.json"
cat "$OUT/features_P/P_SUMMARY.json"
cat "$OUT/features_F/F_SUMMARY.json"
```

这些全部是汇总量，不含参与者标识，也不含任何个体级预测或表型值。

**不要贴**：`tq_split_manifest.tsv`、`traits.tsv`、`P_samples.tsv`、
`P_population.npy`、`F_functional.npy`、`panel_samples.tsv`。
前两个含 UKB eid，后四个是个体级数据。

---

## 6. 读结果的时候记住

`GATE_RESULTS.json` 的 `status` 永远是 `DEVELOPMENTAL_NOT_CONFIRMATORY`。
validation 到这一步已经被用了两次（选 B 的阈值 + 本次评估），所以 PASS 的意思是
"值得写一个 confirmatory 协议并在封存 split 上跑一次"，不是"已经证实"。

`tq_test` 那一次打开，仍然被 withdrawal 名单这一条冻结条款挡着——那条不因为
Concept Gate 跑完而解除。
