# Panel 染色体覆盖范围

_2026-09-17 · 依据全副本一致性校验 · 状态 `20_OF_22_AUTOSOMES`_

---

## 🔍 校验方法

对 `/data3/ukb_all/genotype` 下的每一个 `ukb_imp_chr*_v3.bed` 副本，检验

\[
\text{size} \overset{?}{=} 3 + \left\lceil \frac{n_{\text{samples}}}{4} \right\rceil \times n_{\text{variants}},
\]

其中 `n_variants` 数自该 fileset 自己的 `.bim`，`n_samples` 数自其 `.fam`。PLINK 1 binary 的大小由这两个数完全确定，因此不等即为不完整。

## ❌ 三个残缺副本

| 副本 | 实际字节 | 应为字节 | 完整度 | 有无替代 |
| --- | --- | --- | --- | --- |
| `data1_ukb_crosschr_fixed/ukb_imp_chr11_v3.bed` | 111,621,963,776 | 563,978,088,847 | 19.8% | ✅ `新建文件夹` 那份完好 |
| `新建文件夹/ukb_imp_chr18_v3.bed` | 125,128,671,232 | 316,766,499,890 | 39.5% | ❌ 无第二份 |
| `新建文件夹/ukb_imp_chr21_v3.bed` | 15,789,457,408 | 153,675,885,777 | 10.3% | ❌ 无第二份 |

`data1_ukb_crosschr_fixed` 中的 chr3 与 chr5 副本均完好，因此该目录并非整体损坏，chr11 是单个文件的传输问题。

## 📐 对 panel 的影响

Panel 覆盖 **20 条常染色体**：`chr1`–`chr17`、`chr19`、`chr20`、`chr22`。

排除的 `chr18`（78 Mb）与 `chr21`（48 Mb）合计 126 Mb，占常染色体总长 2,879 Mb 的 **4.4%**。

对本轮资格检验而言这个损失不实质：问题是「raw genotype 相对协变量是否增加样本外信息」，少 4.4% 的基因组会让捕获的方差略微下降，不改变该判断的方向。但它必须写进 `CONFIG.json`，且结论的适用范围相应限定为这 20 条染色体。

## 🚫 为什么不使用残缺文件的前半部分

技术上可以只读 `index < ⌊(size − 3) / bytes_per_variant⌋` 的位点，从而利用 chr18 的前 39.5%。本轮不这么做：那会在 panel 中引入「只取染色体起始段」的位置偏倚，而换来的覆盖增益不到 2%。对一个资格检验来说，干净的排除优于有偏的补全。

## 📋 写入 `CONFIG.json`

```json
"panel": {
  "chromosomes_included": [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,19,20,22],
  "chromosomes_excluded": [18, 21],
  "exclusion_reason": "incomplete .bed transfers with no intact duplicate on this server",
  "excluded_fraction_of_autosome": 0.044
}
```

## 📨 需要告知数据管理方

`新建文件夹` 下的 `ukb_imp_chr18_v3.bed` 与 `ukb_imp_chr21_v3.bed` 是不完整的传输，需要重新拷贝；`data1_ukb_crosschr_fixed` 下的 `ukb_imp_chr11_v3.bed` 同样不完整（该条已有可用替代）。三者的 `.bim` 与 `.fam` 都完整，只有 `.bed` 被截断，符合传输中断而非源数据损坏的特征。
