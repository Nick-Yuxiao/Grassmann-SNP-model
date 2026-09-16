# TQ1 / E0 parallel readiness rc1 结果

## 冻结裁决

最终状态：**NOT_READY_TASK**。

- **E0 PASS**：当前 8 个 local genotype encoders 在固定 locus panels 上明确学到了可泛化的局部 genotype context。
- **TQ1 TASK-INELIGIBLE**：固定 GEUVADIS chr18 expression panel 上，dosage 没有稳定超过 covariates。
- **新 Bridge 未授权**；Grassmann、C/D/E 均未评价。
- 原 230 人 pool 已被 TQ1 消耗，禁止再用于 C/D/E Bridge。

## H:TQ1

230 个原 `bridge_test` individuals 上：

| 指标 | 结果 |
|---|---:|
| macro R²(A) | 0.213539 |
| macro R²(B) | 0.196092 |
| B−A | −0.017447 |
| paired individual-bootstrap 95% CI | [−0.042849, 0.008168] |
| B>A traits | 2/8 |

预先指定的 family-cluster sensitivity 为 −0.017447，95% CI [−0.040912, 0.005596]，与 primary 判定一致。因此旧 50-person Task Gate 的不确定阳性方向没有在更大的 230-person pool 中复现；不能把该 panel 称为合格 Bridge task。

## E0

280 名 development 外 individuals、230 个 families、8 个固定 panels、每 panel 5 个固定 20% masks：

| 指标 | 结果 |
|---|---:|
| empirical marginal CE | 0.847176 |
| real-context model CE | 0.241862 |
| population-matched donor-context CE | 1.726043 |
| empirical marginal accuracy | 0.597285 |
| real-context model accuracy | 0.910832 |
| donor-context accuracy | 0.559998 |
| CE(empirical)−CE(real) | 0.605314; 95% CI [0.593348, 0.617644] |
| CE(donor)−CE(real) | 1.484181; 95% CI [1.426118, 1.543223] |

两个 co-primary contrast 均通过，且 8/8 panels 的方向一致。结果排除了“这批 local encoders 只会输出训练集位点频率、没有使用同一个体上下文”这一解释。

## 科学边界

E0 PASS 只适用于当前 8 个小型、per-locus、40-epoch encoders；它们不是已经完成的 genome-wide UKB foundation model。E0 不证明 hidden states 对 phenotype 有增量价值，也不支持 Grassmann。TQ1 FAIL 只否定当前固定 GEUVADIS task panel 的资格；它不否定 genotype pretraining、LD/MAF/haplotype 学习或所有 phenotype tasks。

下一次只有在更大、预先有稳定 genotype signal 的独立 phenotype 数据上重新完成 A/B qualification 后，才能把通过 E0 的 encoder family 或其规模化版本带入一个全新的 Bridge。

