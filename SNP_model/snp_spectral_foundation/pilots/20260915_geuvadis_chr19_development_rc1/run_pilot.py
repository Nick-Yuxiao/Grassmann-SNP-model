from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import platform
import random
import sys
from dataclasses import asdict
from pathlib import Path
from typing import Iterable

import numpy as np
import sklearn
import torch
from sklearn.decomposition import PCA
from sklearn.model_selection import GroupShuffleSplit


PILOT_DIR = Path(__file__).resolve().parent
ROOT = PILOT_DIR.parents[2]
PACKAGE_ROOT = ROOT / "snp_spectral_foundation"
SRC = PACKAGE_ROOT / "src"
sys.path.insert(0, str(SRC))

from snp_spectral_foundation.config import ModelConfig  # noqa: E402
from snp_spectral_foundation.model import LocalGenotypeEncoder, encode_genotypes  # noqa: E402
from snp_spectral_foundation.training import (  # noqa: E402
    contextual_gate_metrics,
    masked_genotype_loss,
    sample_ssl_mask,
)


SEED = 5_262_026
BOOTSTRAP_SEED = 5_263_026
ALPHAS = np.asarray(
    [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100, 1000, 10000, 100000, 1000000],
    dtype=np.float64,
)
BLOCK_SIZE = 32
MAX_SNPS = 1280
MIN_SNPS = 256
RANK = 4

ASSET_DIR = PACKAGE_ROOT / "pilot_assets" / "geuvadis"
EXTRACTED = ASSET_DIR / "GeuvadisTranscriptExpr" / "inst" / "extdata"
ARCHIVE = ASSET_DIR / "GeuvadisTranscriptExpr_1.40.0.tar.gz"
GENOTYPE_TSV = EXTRACTED / "genotypes_CEU_chr19.tsv"
EXPRESSION_TSV = EXTRACTED / "TrQuantCount_CEU_chr19.tsv"
SNP_SUBSET = EXTRACTED / "snp_id_subset.txt"
GENE_SUBSET = EXTRACTED / "gene_id_subset.txt"
RELATIONSHIPS = ROOT / "data" / "raw" / "relationships_w_pops_041510.txt"
PROTOCOL_DIR = PACKAGE_ROOT / "protocols" / "20260915_0526_stage2_bridge_rc2"
RESULT_DIR = PILOT_DIR / "results"


def sha256_file(path: Path, chunk_size: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def json_dump(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def read_nonempty_lines(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def expression_samples(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        header = next(csv.reader(handle, delimiter="\t"))
    if header[:2] != ["TargetID", "Gene_Symbol"]:
        raise ValueError("unexpected expression header")
    if len(set(header[2:])) != len(header[2:]):
        raise ValueError("duplicate expression sample IDs")
    return header[2:]


def load_relationships(path: Path) -> dict[str, dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return {row["IID"]: row for row in csv.DictReader(handle, delimiter="\t")}


def make_group_split(
    sample_ids: list[str], relationships: dict[str, dict[str, str]], seed: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    groups = [
        relationships[sample]["FID"]
        if sample in relationships and relationships[sample].get("FID")
        else f"singleton:{sample}"
        for sample in sample_ids
    ]
    all_indices = np.arange(len(sample_ids))
    first = GroupShuffleSplit(n_splits=1, test_size=0.15, random_state=seed)
    train_val_rel, test_rel = next(first.split(all_indices, groups=np.asarray(groups)))
    train_val = all_indices[train_val_rel]
    test = all_indices[test_rel]
    train_val_groups = np.asarray(groups)[train_val]
    second = GroupShuffleSplit(n_splits=1, test_size=0.15 / 0.85, random_state=seed + 1)
    train_rel, val_rel = next(second.split(train_val, groups=train_val_groups))
    train = train_val[train_rel]
    val = train_val[val_rel]
    if min(len(train), len(val), len(test)) == 0:
        raise ValueError("empty split")
    group_sets = [set(np.asarray(groups)[idx]) for idx in (train, val, test)]
    if group_sets[0] & group_sets[1] or group_sets[0] & group_sets[2] or group_sets[1] & group_sets[2]:
        raise AssertionError("family leakage across splits")
    return np.sort(train), np.sort(val), np.sort(test), groups


def load_and_select_trait(
    path: Path,
    expected_samples: list[str],
    candidate_genes: set[str],
    train_idx: np.ndarray,
) -> tuple[np.ndarray, dict[str, object]]:
    best: tuple[float, str, str, np.ndarray, float] | None = None
    n_candidates = 0
    n_qc = 0
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if header[2:] != expected_samples:
            raise ValueError("expression sample order changed")
        for row in reader:
            if len(row) != len(header) or row[1] not in candidate_genes:
                continue
            n_candidates += 1
            raw = np.asarray(row[2:], dtype=np.float64)
            train_raw = raw[train_idx]
            if not np.isfinite(train_raw).all() or float(np.mean(train_raw > 0)) < 0.80:
                continue
            transformed = np.log1p(raw)
            variance = float(np.var(transformed[train_idx], ddof=0))
            n_qc += 1
            candidate = (variance, row[0], row[1], transformed, float(np.mean(train_raw > 0)))
            if best is None or variance > best[0] or (variance == best[0] and row[0] < best[1]):
                best = candidate
    if best is None:
        raise ValueError("no transcript passed frozen train-only trait QC")
    return best[3], {
        "target_id": best[1],
        "gene_id": best[2],
        "transform": "log1p_expected_transcript_count",
        "train_variance": best[0],
        "train_nonzero_fraction": best[4],
        "candidate_transcripts": n_candidates,
        "qc_passing_transcripts": n_qc,
        "selection_data_scope": "train_only",
    }


def load_candidate_genotypes(
    path: Path, expected_samples: list[str], desired_ids: set[str]
) -> tuple[np.ndarray, list[dict[str, object]]]:
    rows: list[np.ndarray] = []
    meta: list[dict[str, object]] = []
    last_position = -1
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if header[:4] != ["chr", "start", "end", "snpId"] or header[4:] != expected_samples:
            raise ValueError("unexpected genotype header or sample order")
        for row in reader:
            if len(row) != len(header):
                raise ValueError("malformed genotype row")
            chrom, start_text, end_text, snp_id = row[:4]
            position = int(start_text)
            if chrom != "19" or position != int(end_text):
                raise ValueError("coordinate contract violation")
            if position < last_position:
                raise ValueError("genotype source is not coordinate sorted")
            last_position = position
            if snp_id not in desired_ids:
                continue
            dosage = np.asarray(row[4:], dtype=np.int8)
            if np.any((dosage < -1) | (dosage > 2)):
                raise ValueError(f"invalid dosage at {snp_id}")
            rows.append(dosage)
            meta.append({"chrom": chrom, "position_1based": position, "snp_id": snp_id})
    found = {str(item["snp_id"]) for item in meta}
    missing = desired_ids - found
    if missing:
        raise ValueError(f"{len(missing)} frozen subset SNP IDs were absent")
    return np.stack(rows, axis=1), meta


def choose_panel_indices(n_eligible: int, max_snps: int = MAX_SNPS, block_size: int = BLOCK_SIZE) -> np.ndarray:
    n_use = min(max_snps, (n_eligible // block_size) * block_size)
    if n_use < MIN_SNPS:
        raise ValueError(f"only {n_use} eligible SNPs; frozen minimum is {MIN_SNPS}")
    if n_use == n_eligible:
        return np.arange(n_use, dtype=np.int64)
    selected = np.linspace(0, n_eligible - 1, num=n_use, dtype=np.int64)
    if len(np.unique(selected)) != n_use:
        raise AssertionError("evenly spaced panel selection produced duplicates")
    return selected


def qc_and_select_panel(
    genotype: np.ndarray, meta: list[dict[str, object]], train_idx: np.ndarray
) -> tuple[np.ndarray, list[dict[str, object]], np.ndarray, dict[str, object]]:
    train = genotype[train_idx]
    observed = train >= 0
    counts = observed.sum(axis=0)
    missingness = 1.0 - counts / len(train_idx)
    alt_frequency = np.divide(
        np.where(observed, train, 0).sum(axis=0),
        2.0 * counts,
        out=np.full(genotype.shape[1], np.nan, dtype=np.float64),
        where=counts > 0,
    )
    maf = np.minimum(alt_frequency, 1.0 - alt_frequency)
    eligible = np.flatnonzero((missingness <= 0.10) & (maf >= 0.05) & np.isfinite(maf))
    take = choose_panel_indices(len(eligible))
    selected = eligible[take]
    panel = genotype[:, selected]
    panel_meta = [dict(meta[i]) for i in selected]
    panel_af = alt_frequency[selected]
    for item, miss, af_value in zip(panel_meta, missingness[selected], panel_af, strict=True):
        item["train_missingness"] = float(miss)
        item["train_alt_frequency"] = float(af_value)
        item["train_maf"] = float(min(af_value, 1.0 - af_value))
    return panel, panel_meta, panel_af, {
        "candidate_snps": int(genotype.shape[1]),
        "eligible_snps": int(len(eligible)),
        "selected_snps": int(len(selected)),
        "blocks": int(len(selected) // BLOCK_SIZE),
        "filter_scope": "train_only",
        "missingness_max": 0.10,
        "maf_min": 0.05,
    }


def impute_dosage(panel: np.ndarray, panel_af: np.ndarray) -> np.ndarray:
    result = panel.astype(np.float32)
    missing = result < 0
    imputed = np.broadcast_to((2.0 * panel_af).astype(np.float32), result.shape)
    result[missing] = imputed[missing]
    return result


def fixed_covariates(
    sample_ids: list[str], relationships: dict[str, dict[str, str]], dosage: np.ndarray, train_idx: np.ndarray
) -> tuple[np.ndarray, dict[str, object]]:
    sex = np.full(len(sample_ids), np.nan, dtype=np.float64)
    for i, sample in enumerate(sample_ids):
        value = relationships.get(sample, {}).get("sex")
        if value in {"1", "2"}:
            sex[i] = 1.0 if value == "2" else 0.0
    sex_missing = ~np.isfinite(sex)
    train_sex_mean = float(np.nanmean(sex[train_idx]))
    sex[sex_missing] = train_sex_mean

    geno_mean = dosage[train_idx].mean(axis=0, dtype=np.float64)
    geno_sd = dosage[train_idx].std(axis=0, dtype=np.float64)
    geno_sd[geno_sd < 1e-8] = 1.0
    standardized = (dosage.astype(np.float64) - geno_mean) / geno_sd
    n_components = min(10, len(train_idx) - 1, dosage.shape[1])
    pca = PCA(n_components=n_components, svd_solver="full")
    pca.fit(standardized[train_idx])
    pcs = pca.transform(standardized)
    covariates = np.column_stack((pcs, sex, sex_missing.astype(np.float64))).astype(np.float32)
    return covariates, {
        "genotype_pc_count": int(n_components),
        "train_pc_explained_variance_ratio": [float(x) for x in pca.explained_variance_ratio_],
        "sex_coding": "female=1_male=0",
        "sex_observed": int((~sex_missing).sum()),
        "sex_missing": int(sex_missing.sum()),
        "sex_train_mean_imputation": train_sex_mean,
    }


def set_reproducibility(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.set_num_threads(min(4, max(1, os.cpu_count() or 1)))
    torch.use_deterministic_algorithms(True)


def train_encoder(
    panel: np.ndarray,
    panel_af: np.ndarray,
    train_idx: np.ndarray,
    val_idx: np.ndarray,
) -> tuple[LocalGenotypeEncoder, ModelConfig, list[dict[str, float]], dict[str, float]]:
    n_blocks = panel.shape[1] // BLOCK_SIZE
    cfg = ModelConfig(
        max_blocks=n_blocks,
        snps_per_block=BLOCK_SIZE,
        d_model=16,
        n_heads=4,
        local_layers=2,
        local_ff_dim=64,
        spectral_rank=RANK,
        spectral_oversample=4,
        dropout=0.0,
    )
    model = LocalGenotypeEncoder(cfg)
    optimizer = torch.optim.AdamW(model.parameters(), lr=0.002, weight_decay=0.0001)
    geno = torch.from_numpy(panel.reshape(len(panel), n_blocks, BLOCK_SIZE).astype(np.int64))
    af = torch.from_numpy(panel_af.reshape(n_blocks, BLOCK_SIZE).astype(np.float32))
    mask_generator = torch.Generator(device="cpu").manual_seed(SEED + 100)
    order_rng = np.random.default_rng(SEED + 101)
    val_geno = geno[val_idx]
    val_mask_generator = torch.Generator(device="cpu").manual_seed(SEED + 102)
    val_mask = sample_ssl_mask(val_geno, 0.20, val_mask_generator)
    history: list[dict[str, float]] = []
    for epoch in range(1, 41):
        model.train()
        losses: list[float] = []
        for batch_idx in np.array_split(order_rng.permutation(train_idx), math.ceil(len(train_idx) / 16)):
            batch = geno[batch_idx]
            mask = sample_ssl_mask(batch, 0.20, mask_generator)
            hidden, logits = model(encode_genotypes(batch, mask), af)
            del hidden
            loss = masked_genotype_loss(logits, batch, mask)
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            optimizer.step()
            losses.append(float(loss.detach()))
        model.eval()
        with torch.no_grad():
            _, val_logits = model(encode_genotypes(val_geno, val_mask), af)
            val_loss = float(masked_genotype_loss(val_logits, val_geno, val_mask))
        history.append({"epoch": float(epoch), "train_loss": float(np.mean(losses)), "val_loss": val_loss})
        if epoch == 1 or epoch % 5 == 0:
            print(f"pretrain epoch {epoch:02d}/40 train_loss={np.mean(losses):.6f} val_loss={val_loss:.6f}", flush=True)
    model.eval()
    with torch.no_grad():
        _, val_logits = model(encode_genotypes(val_geno, val_mask), af)
        metrics = contextual_gate_metrics(val_logits, val_geno, val_mask, af)
        metrics["masked_cross_entropy"] = float(masked_genotype_loss(val_logits, val_geno, val_mask))
    return model, cfg, history, metrics


def extract_hidden(
    model: LocalGenotypeEncoder, panel: np.ndarray, panel_af: np.ndarray, batch_size: int = 16
) -> np.ndarray:
    n_blocks = panel.shape[1] // BLOCK_SIZE
    geno = torch.from_numpy(panel.reshape(len(panel), n_blocks, BLOCK_SIZE).astype(np.int64))
    af = torch.from_numpy(panel_af.reshape(n_blocks, BLOCK_SIZE).astype(np.float32))
    chunks: list[np.ndarray] = []
    model.eval()
    with torch.no_grad():
        for start in range(0, len(panel), batch_size):
            hidden, _ = model(encode_genotypes(geno[start : start + batch_size]), af)
            chunks.append(hidden.cpu().numpy().astype(np.float32))
    return np.concatenate(chunks, axis=0)


def aligned_hidden_pca(hidden: np.ndarray, train_idx: np.ndarray, rank: int) -> tuple[np.ndarray, dict[str, object]]:
    train_tokens = hidden[train_idx].reshape(-1, hidden.shape[-1]).astype(np.float64)
    mean = train_tokens.mean(axis=0)
    covariance = np.cov(train_tokens - mean, rowvar=False, ddof=1)
    values, vectors = np.linalg.eigh(covariance)
    order = np.argsort(values)[::-1][:rank]
    basis = vectors[:, order]
    scores = np.einsum("nbsd,dr->nbsr", hidden.astype(np.float64) - mean, basis).astype(np.float32)
    return scores.reshape(len(hidden), -1), {
        "rank": rank,
        "train_hidden_channel_eigenvalues": [float(values[i]) for i in order],
        "fit_scope": "train_only",
    }


def grassmann_features(hidden: np.ndarray, rank: int) -> np.ndarray:
    centered = hidden.astype(np.float64) - hidden.astype(np.float64).mean(axis=2, keepdims=True)
    covariance = np.einsum("nbsd,nbse->nbde", centered, centered) / float(hidden.shape[2] - 1)
    values, vectors = np.linalg.eigh(covariance)
    values = values[..., -rank:][..., ::-1].copy()
    vectors = vectors[..., -rank:][..., ::-1].copy()
    projector = np.einsum("nbdr,nber->nbde", vectors, vectors)
    triangle = np.triu_indices(hidden.shape[-1])
    projector_features = projector[..., triangle[0], triangle[1]]
    spectrum = np.log1p(np.clip(values, 0.0, None))
    spectrum /= np.maximum(spectrum.sum(axis=-1, keepdims=True), 1e-12)
    return np.concatenate((projector_features, spectrum), axis=-1).reshape(len(hidden), -1).astype(np.float32)


def fit_ridge_arm(
    x: np.ndarray,
    y: np.ndarray,
    train_idx: np.ndarray,
    val_idx: np.ndarray,
    test_idx: np.ndarray,
) -> tuple[np.ndarray, dict[str, object]]:
    x_train = x[train_idx].astype(np.float64)
    mean = x_train.mean(axis=0)
    sd = x_train.std(axis=0)
    keep = sd >= 1e-8
    if not np.any(keep):
        raise ValueError("arm has no nonconstant features")
    z_train = (x_train[:, keep] - mean[keep]) / sd[keep]
    z_val = (x[val_idx, :][:, keep].astype(np.float64) - mean[keep]) / sd[keep]
    z_test = (x[test_idx, :][:, keep].astype(np.float64) - mean[keep]) / sd[keep]
    y_mean = float(y[train_idx].mean())
    y_sd = float(y[train_idx].std(ddof=0))
    if y_sd < 1e-12:
        raise ValueError("selected trait is constant in train")
    y_train_std = (y[train_idx] - y_mean) / y_sd
    gram = z_train @ z_train.T
    val_cross = z_val @ z_train.T
    test_cross = z_test @ z_train.T
    identity = np.eye(len(train_idx), dtype=np.float64)
    best: tuple[float, float, np.ndarray] | None = None
    for alpha in ALPHAS:
        dual = np.linalg.solve(gram + float(alpha) * identity, y_train_std)
        val_prediction = y_mean + y_sd * (val_cross @ dual)
        mse = float(np.mean((y[val_idx] - val_prediction) ** 2))
        if best is None or mse < best[0] - 1e-12 or (abs(mse - best[0]) <= 1e-12 and alpha > best[1]):
            best = (mse, float(alpha), dual)
    assert best is not None
    test_prediction = y_mean + y_sd * (test_cross @ best[2])
    return test_prediction, {
        "input_features": int(x.shape[1]),
        "nonconstant_train_features": int(keep.sum()),
        "selected_alpha": best[1],
        "validation_mse": best[0],
        "fit_scope": "train_only_alpha_on_validation",
    }


def predictive_r2(y_true: np.ndarray, prediction: np.ndarray, train_mean: float) -> float:
    denominator = float(np.sum((y_true - train_mean) ** 2))
    if denominator <= 0:
        raise ValueError("nonpositive predictive-R2 denominator")
    return 1.0 - float(np.sum((y_true - prediction) ** 2)) / denominator


def bootstrap_statistics(
    y_test: np.ndarray,
    predictions: dict[str, np.ndarray],
    train_mean: float,
    replicates: int = 5000,
    seed: int = BOOTSTRAP_SEED,
) -> tuple[dict[str, dict[str, float]], dict[str, dict[str, float]]]:
    rng = np.random.default_rng(seed)
    sample_indices = rng.integers(0, len(y_test), size=(replicates, len(y_test)))
    y_boot = y_test[sample_indices]
    denominator = np.sum((y_boot - train_mean) ** 2, axis=1)
    if np.any(denominator <= 0):
        raise ValueError("degenerate bootstrap replicate")
    boot_r2: dict[str, np.ndarray] = {}
    arm_summary: dict[str, dict[str, float]] = {}
    for arm, prediction in predictions.items():
        residual = y_boot - prediction[sample_indices]
        values = 1.0 - np.sum(residual**2, axis=1) / denominator
        boot_r2[arm] = values
        arm_summary[arm] = {
            "r2": predictive_r2(y_test, prediction, train_mean),
            "ci95_low": float(np.quantile(values, 0.025)),
            "ci95_high": float(np.quantile(values, 0.975)),
        }
    contrast_pairs = {"C-B": ("C", "B"), "B-A": ("B", "A"), "D-B": ("D", "B"), "E-B": ("E", "B"), "E-D": ("E", "D")}
    contrast_summary: dict[str, dict[str, float]] = {}
    for label, (left, right) in contrast_pairs.items():
        values = boot_r2[left] - boot_r2[right]
        contrast_summary[label] = {
            "delta_r2": arm_summary[left]["r2"] - arm_summary[right]["r2"],
            "ci95_low": float(np.quantile(values, 0.025)),
            "ci95_high": float(np.quantile(values, 0.975)),
        }
    return arm_summary, contrast_summary


def write_panel(path: Path, panel_meta: Iterable[dict[str, object]]) -> None:
    fields = ["chrom", "position_1based", "snp_id", "train_missingness", "train_alt_frequency", "train_maf"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(panel_meta)


def write_predictions(
    path: Path, sample_ids: list[str], test_idx: np.ndarray, y: np.ndarray, predictions: dict[str, np.ndarray]
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        fields = ["sample_id", "observed_log1p_count", "pred_A", "pred_B", "pred_C", "pred_D", "pred_E"]
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        for offset, index in enumerate(test_idx):
            writer.writerow([sample_ids[index], f"{y[index]:.12g}"] + [f"{predictions[a][offset]:.12g}" for a in "ABCDE"])


def report_text(results: dict[str, object]) -> str:
    primary = results["contrasts"]["C-B"]
    disposition = "支持" if results["primary_support"] else "不支持"
    rows = []
    for arm in "ABCDE":
        metric = results["arms"][arm]
        rows.append(f"| {arm} | {metric['r2']:.6f} | [{metric['ci95_low']:.6f}, {metric['ci95_high']:.6f}] | {results['decoders'][arm]['selected_alpha']:g} |")
    return f"""# GEUVADIS chr19 development Pilot rc1 结果

## 冻结主结论

本次结果**{disposition} 0526 Stage 2 bridge**。冻结主比较 `C-B` 的
ΔR² = **{primary['delta_r2']:.6f}**，paired 95% bootstrap CI =
[{primary['ci95_low']:.6f}, {primary['ci95_high']:.6f}]。判定严格使用运行前规则：
点估计大于 0 且区间下界大于 0；没有追加 0.005 门槛、P1 前置门或多重校正门。

## 数据与模型

- 数据：Bioconductor `GeuvadisTranscriptExpr 1.40.0`，CEU chr19，真实个体级 genotype–transcript expression。
- 样本：train/validation/test = {results['split']['n_train']}/{results['split']['n_validation']}/{results['split']['n_test']}，家系隔离。
- 性状：`{results['trait']['target_id']}` / `{results['trait']['gene_id']}`，仅由 train 选择。
- SNP：{results['panel']['selected_snps']} 个，{results['panel']['blocks']} blocks；缺失填补、AF、MAF 和 PCA 只由 train 拟合。
- Encoder：40 epochs masked-genotype pretraining；checkpoint 在提取 H 和 phenotype decoder 前保存并冻结。

## 各臂 held-out predictive R²

| Arm | R² | paired-bootstrap marginal 95% CI | alpha |
| --- | ---: | ---: | ---: |
{chr(10).join(rows)}

## 解释边界

这是小样本真实分子表型 development Pilot，不是外部复制、全基因组有效性、临床预测或因果证据。
`E` 是 Grassmann 候选机制臂；无论它输赢，都不能改写唯一主问题 `C-B`。后续是否进入新数据或 biology-prior
阶段，只按本文件已经给出的冻结结果决定，不事后修改本次标准。
"""


def run() -> None:
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    if (RESULT_DIR / "RESULTS.json").exists():
        raise RuntimeError("RESULTS.json already exists; frozen Pilot cannot be rerun")
    required = [ARCHIVE, GENOTYPE_TSV, EXPRESSION_TSV, SNP_SUBSET, GENE_SUBSET, RELATIONSHIPS]
    for path in required:
        if not path.is_file():
            raise FileNotFoundError(path)
    archive_hash = sha256_file(ARCHIVE)
    expected_hash = "619bf74ac7d9f226ffd8ccde218bb40bd32889f66336294ce91ac5fb9b845d6d"
    if archive_hash != expected_hash:
        raise ValueError("Bioconductor source archive hash mismatch")

    set_reproducibility(SEED)
    samples = expression_samples(EXPRESSION_TSV)
    relationships = load_relationships(RELATIONSHIPS)
    train_idx, val_idx, test_idx, groups = make_group_split(samples, relationships, SEED)
    print(f"split train/validation/test={len(train_idx)}/{len(val_idx)}/{len(test_idx)}", flush=True)

    candidate_genes = set(read_nonempty_lines(GENE_SUBSET))
    y, trait_info = load_and_select_trait(EXPRESSION_TSV, samples, candidate_genes, train_idx)
    print(f"train-only selected trait={trait_info['target_id']} gene={trait_info['gene_id']}", flush=True)

    desired_snps = set(read_nonempty_lines(SNP_SUBSET))
    genotype, genotype_meta = load_candidate_genotypes(GENOTYPE_TSV, samples, desired_snps)
    panel, panel_meta, panel_af, panel_info = qc_and_select_panel(genotype, genotype_meta, train_idx)
    print(f"panel candidates/eligible/selected={panel_info['candidate_snps']}/{panel_info['eligible_snps']}/{panel_info['selected_snps']}", flush=True)
    write_panel(RESULT_DIR / "PANEL.tsv", panel_meta)

    dosage = impute_dosage(panel, panel_af)
    covariates, covariate_info = fixed_covariates(samples, relationships, dosage, train_idx)
    model, config, pretrain_history, pretrain_metrics = train_encoder(panel, panel_af, train_idx, val_idx)

    checkpoint_path = RESULT_DIR / "ENCODER_FROZEN.pt"
    torch.save(
        {
            "state_dict": model.state_dict(),
            "config": asdict(config),
            "seed": SEED,
            "panel_snp_ids": [item["snp_id"] for item in panel_meta],
            "train_sample_ids": [samples[i] for i in train_idx],
            "source_archive_sha256": archive_hash,
            "phenotype_seen_by_encoder": False,
        },
        checkpoint_path,
    )
    checkpoint_hash = sha256_file(checkpoint_path)
    for parameter in model.parameters():
        parameter.requires_grad_(False)
    print(f"encoder frozen checkpoint_sha256={checkpoint_hash}", flush=True)

    hidden = extract_hidden(model, panel, panel_af)
    aligned, aligned_info = aligned_hidden_pca(hidden, train_idx, RANK)
    grassmann = grassmann_features(hidden, RANK)
    features = {
        "A": covariates,
        "B": np.concatenate((covariates, dosage), axis=1),
        "C": np.concatenate((covariates, hidden.reshape(len(samples), -1)), axis=1),
        "D": np.concatenate((covariates, aligned), axis=1),
        "E": np.concatenate((covariates, grassmann), axis=1),
    }
    predictions: dict[str, np.ndarray] = {}
    decoder_info: dict[str, dict[str, object]] = {}
    for arm in "ABCDE":
        predictions[arm], decoder_info[arm] = fit_ridge_arm(features[arm], y, train_idx, val_idx, test_idx)
        print(f"arm {arm}: alpha={decoder_info[arm]['selected_alpha']:g} val_mse={decoder_info[arm]['validation_mse']:.6g}", flush=True)

    train_mean = float(y[train_idx].mean())
    arm_summary, contrasts = bootstrap_statistics(y[test_idx], predictions, train_mean)
    primary_support = bool(contrasts["C-B"]["delta_r2"] > 0 and contrasts["C-B"]["ci95_low"] > 0)
    split_info = {
        "seed": SEED,
        "unit": "pedigree_family_or_singleton",
        "n_total": len(samples),
        "n_train": len(train_idx),
        "n_validation": len(val_idx),
        "n_test": len(test_idx),
        "train_sample_ids": [samples[i] for i in train_idx],
        "validation_sample_ids": [samples[i] for i in val_idx],
        "test_sample_ids": [samples[i] for i in test_idx],
        "train_family_count": len({groups[i] for i in train_idx}),
        "validation_family_count": len({groups[i] for i in val_idx}),
        "test_family_count": len({groups[i] for i in test_idx}),
    }
    results: dict[str, object] = {
        "analysis_id": "0526-stage2-geuvadis-development-pilot-rc2",
        "pilot_id": "20260915_geuvadis_chr19_development_rc1",
        "status": "VALID_COMPLETE",
        "claim_class": "real_molecular_phenotype_development_pilot_not_external_replication",
        "split": split_info,
        "trait": trait_info,
        "panel": panel_info,
        "covariates": covariate_info,
        "encoder": {
            "config": asdict(config),
            "epochs": 40,
            "mask_probability": 0.20,
            "batch_size": 16,
            "optimizer": "AdamW",
            "learning_rate": 0.002,
            "weight_decay": 0.0001,
            "validation_metrics": pretrain_metrics,
            "checkpoint_sha256": checkpoint_hash,
            "phenotype_seen_during_pretraining": False,
        },
        "D_aligned_control": aligned_info,
        "E_grassmann": {
            "rank": RANK,
            "feature": "projector_upper_triangle_plus_normalized_log1p_eigenvalues",
            "basis_invariant": True,
        },
        "decoders": decoder_info,
        "arms": arm_summary,
        "contrasts": contrasts,
        "primary_contrast": "C-B",
        "primary_support_rule": "delta_r2_gt_0_and_paired_ci95_low_gt_0",
        "primary_support": primary_support,
        "bootstrap_replicates": 5000,
        "bootstrap_seed": BOOTSTRAP_SEED,
        "test_open_count": 1,
        "interpretation_limits": [
            "not_external_replication",
            "not_clinical_prediction",
            "not_genome_wide_validation",
            "not_causal_evidence",
            "E_does_not_gate_primary_C_minus_B",
        ],
    }
    write_predictions(RESULT_DIR / "TEST_PREDICTIONS.tsv", samples, test_idx, y, predictions)
    json_dump(RESULT_DIR / "PRETRAIN_HISTORY.json", pretrain_history)
    json_dump(RESULT_DIR / "RESULTS.json", results)
    (RESULT_DIR / "REPORT.zh-CN.md").write_text(report_text(results), encoding="utf-8")

    binding = {
        "source_archive": str(ARCHIVE.relative_to(ROOT)),
        "source_archive_sha256": archive_hash,
        "extracted_input_sha256": {
            str(path.relative_to(ROOT)): sha256_file(path)
            for path in [GENOTYPE_TSV, EXPRESSION_TSV, SNP_SUBSET, GENE_SUBSET, RELATIONSHIPS]
        },
        "protocol_sha256": {
            path.name: sha256_file(path)
            for path in sorted(PROTOCOL_DIR.iterdir())
            if path.is_file()
        },
        "runtime": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scikit_learn": sklearn.__version__,
            "torch": torch.__version__,
            "cuda_available": torch.cuda.is_available(),
            "torch_threads": torch.get_num_threads(),
        },
    }
    json_dump(RESULT_DIR / "RUN_BINDING.json", binding)
    manifest_paths = [
        Path(__file__),
        PILOT_DIR / "README.md",
        *sorted(PROTOCOL_DIR.glob("*")),
        checkpoint_path,
        RESULT_DIR / "PANEL.tsv",
        RESULT_DIR / "TEST_PREDICTIONS.tsv",
        RESULT_DIR / "PRETRAIN_HISTORY.json",
        RESULT_DIR / "RESULTS.json",
        RESULT_DIR / "REPORT.zh-CN.md",
        RESULT_DIR / "RUN_BINDING.json",
    ]
    manifest_lines = [f"{sha256_file(path)}  {path.relative_to(ROOT).as_posix()}" for path in manifest_paths]
    (RESULT_DIR / "RUN_MANIFEST.sha256").write_text("\n".join(manifest_lines) + "\n", encoding="utf-8")
    print(
        f"primary C-B delta_R2={contrasts['C-B']['delta_r2']:.6f} "
        f"CI=[{contrasts['C-B']['ci95_low']:.6f}, {contrasts['C-B']['ci95_high']:.6f}] "
        f"support={primary_support}",
        flush=True,
    )


if __name__ == "__main__":
    run()
