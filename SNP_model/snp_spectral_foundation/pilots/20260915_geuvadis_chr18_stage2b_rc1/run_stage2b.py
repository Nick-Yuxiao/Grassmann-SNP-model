from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import os
import platform
import random
import sys
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pgenlib
import sklearn
import torch
from sklearn.decomposition import PCA


PILOT_DIR = Path(__file__).resolve().parent
ROOT = PILOT_DIR.parents[2]
PACKAGE = ROOT / "snp_spectral_foundation"
SRC = PACKAGE / "src"
sys.path.insert(0, str(SRC))

from snp_spectral_foundation.config import ModelConfig  # noqa: E402
from snp_spectral_foundation.model import LocalGenotypeEncoder, encode_genotypes  # noqa: E402
from snp_spectral_foundation.training import (  # noqa: E402
    contextual_gate_metrics,
    masked_genotype_loss,
    sample_ssl_mask,
)


ASSET = PACKAGE / "pilot_assets" / "geuvadis_tensorqtl_chr18"
STEM = "GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.chr18"
PGEN = ASSET / f"{STEM}.pgen"
PVAR = ASSET / f"{STEM}.pvar"
PSAM = ASSET / f"{STEM}.psam"
EXPRESSION = ASSET / "GEUVADIS.445_samples.expression.bed.gz"
SAMPLE_PANEL = ASSET / "integrated_call_samples_v3.20130502.ALL.panel"
PEDIGREE = ASSET / "integrated_call_samples_v3.20250704.ALL.ped"
MANIFEST = PILOT_DIR / "SAMPLE_MANIFEST.tsv"
CONFIG_PATH = PILOT_DIR / "CONFIG_FROZEN.json"
PROTOCOL_PATH = PILOT_DIR / "FROZEN_PROTOCOL.zh-CN.md"
DESIGN_PATH = PILOT_DIR / "DESIGN_FREEZE.json"
RESULT_DIR = PILOT_DIR / "results"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def dump_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def set_reproducibility(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.set_num_threads(min(4, max(1, os.cpu_count() or 1)))
    torch.use_deterministic_algorithms(True)


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_pvar() -> tuple[np.ndarray, list[str], list[str], list[str]]:
    positions: list[int] = []
    ids: list[str] = []
    refs: list[str] = []
    alts: list[str] = []
    header = None
    with PVAR.open("r", encoding="utf-8", newline="") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line.rstrip("\n").lstrip("#").split("\t")
                if header[:5] != ["CHROM", "POS", "ID", "REF", "ALT"]:
                    raise ValueError("unexpected PVAR header")
                continue
            row = line.rstrip("\n").split("\t")
            if header is None or row[0] != "chr18":
                raise ValueError("PVAR coordinate/contig contract violation")
            pos = int(row[1])
            if positions and pos < positions[-1]:
                raise ValueError("PVAR is not coordinate sorted")
            if "," in row[4]:
                raise ValueError("multiallelic PVAR record")
            positions.append(pos)
            ids.append(row[2])
            refs.append(row[3])
            alts.append(row[4])
    if len(positions) != 367_759 or len(set(ids)) != len(ids):
        raise ValueError("unexpected PVAR count or duplicate IDs")
    return np.asarray(positions, dtype=np.int64), ids, refs, alts


def read_pgen(n_samples: int, n_variants: int) -> np.ndarray:
    genotype = np.empty((n_samples, n_variants), dtype=np.int8)
    with pgenlib.PgenReader(str(PGEN).encode()) as reader:
        if reader.get_raw_sample_ct() != n_samples or reader.get_variant_ct() != n_variants:
            raise ValueError("PGEN/PVAR/PSAM dimension mismatch")
        reader.read_range(0, n_variants, genotype, sample_maj=1)
    if np.any((genotype < -9) | (genotype > 2)):
        raise ValueError("invalid hardcall domain")
    return genotype


def expression_header() -> list[str]:
    with gzip.open(EXPRESSION, "rt", encoding="utf-8", newline="") as handle:
        header = next(csv.reader(handle, delimiter="\t"))
    if header[:4] != ["#chr", "start", "end", "gene_id"]:
        raise ValueError("unexpected expression BED header")
    return header[4:]


def load_expression_subset(sample_ids: list[str], selected_genes: set[str] | None = None) -> list[dict[str, object]]:
    all_ids = expression_header()
    lookup = {sample: i + 4 for i, sample in enumerate(all_ids)}
    columns = [lookup[sample] for sample in sample_ids]
    records: list[dict[str, object]] = []
    with gzip.open(EXPRESSION, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader)
        for row in reader:
            if row[0] != "chr18":
                continue
            gene = row[3]
            if selected_genes is not None and gene not in selected_genes:
                continue
            start = int(row[1])
            end = int(row[2])
            if end != start + 1:
                raise ValueError(f"expression row is not a one-base TSS interval: {gene}")
            values = np.asarray([float(row[i]) for i in columns], dtype=np.float64)
            records.append({"gene_id": gene, "tss_1based": start + 1, "values": values})
    if selected_genes is not None and {str(x["gene_id"]) for x in records} != selected_genes:
        raise ValueError("selected expression traits were not recovered exactly")
    return records


def af_and_dosage(genotype: np.ndarray, train_idx: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    train = genotype[train_idx]
    observed = train >= 0
    count = observed.sum(axis=0)
    af = np.divide(np.where(observed, train, 0).sum(axis=0), 2.0 * count, out=np.full(genotype.shape[1], np.nan), where=count > 0)
    missingness = 1.0 - count / len(train_idx)
    dosage = genotype.astype(np.float32)
    missing = dosage < 0
    dosage[missing] = np.broadcast_to((2.0 * af).astype(np.float32), dosage.shape)[missing]
    return af, missingness, dosage


def correlations(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    xc = x.astype(np.float64) - x.mean(axis=0, dtype=np.float64)
    yc = y.astype(np.float64) - y.mean()
    denom = np.sqrt(np.sum(xc * xc, axis=0) * np.sum(yc * yc))
    return np.divide(xc.T @ yc, denom, out=np.zeros(x.shape[1], dtype=np.float64), where=denom > 0)


def base_demographics(rows: list[dict[str, str]], indices: np.ndarray, fit_indices: np.ndarray) -> np.ndarray:
    populations = ["CEU", "FIN", "GBR", "TSI", "YRI"]
    if set(row["pop"] for row in rows if row["primary_role"] != "EXCLUDED_RC2_TEST") != set(populations):
        raise ValueError("unexpected population set")
    sex = np.asarray([1.0 if row["sex"] == "2" else 0.0 for row in rows], dtype=np.float64)
    pop = np.asarray([[float(row["pop"] == p) for p in populations[1:]] for row in rows], dtype=np.float64)
    design = np.column_stack([np.ones(len(rows)), sex, pop])
    # The arguments document that no held-out rows are used for fitting; demographics themselves need no fitted transform.
    if not set(fit_indices).issubset(set(indices)) and len(indices) != len(rows):
        raise ValueError("invalid demographic fit scope")
    return design[:, 1:].astype(np.float32)


def residualize_demographics(y: np.ndarray, demographics: np.ndarray, dev_indices: np.ndarray) -> np.ndarray:
    design = np.column_stack([np.ones(len(dev_indices)), demographics[dev_indices].astype(np.float64)])
    beta = np.linalg.lstsq(design, y, rcond=None)[0]
    return y - design @ beta


def select_traits_and_panels(
    records: list[dict[str, object]],
    genotype: np.ndarray,
    positions: np.ndarray,
    variant_ids: list[str],
    dev_indices: np.ndarray,
    train_indices: np.ndarray,
    demographics: np.ndarray,
    config: dict[str, object],
) -> list[dict[str, object]]:
    af_all, missing_all, dosage_all = af_and_dosage(genotype, train_indices)
    maf_all = np.minimum(af_all, 1.0 - af_all)
    candidates: list[dict[str, object]] = []
    for record in records:
        y = np.asarray(record["values"], dtype=np.float64)
        if len(y) != len(dev_indices) or not np.isfinite(y).all() or float(np.var(y)) <= 0:
            continue
        tss = int(record["tss_1based"])
        cis = np.flatnonzero(
            (positions >= tss - int(config["cis_radius_bp"]))
            & (positions <= tss + int(config["cis_radius_bp"]))
            & (maf_all >= float(config["maf_min"]))
            & (missing_all <= float(config["missingness_max"]))
        )
        if len(cis) < int(config["snps_per_trait"]):
            continue
        r = correlations(dosage_all[np.ix_(dev_indices, cis)], y)
        order = sorted(range(len(cis)), key=lambda j: (-abs(float(r[j])), variant_ids[int(cis[j])]))
        chosen = cis[np.asarray(order[: int(config["snps_per_trait"])], dtype=np.int64)]
        chosen = chosen[np.argsort(positions[chosen], kind="stable")]
        candidates.append(
            {
                "gene_id": str(record["gene_id"]),
                "tss_1based": tss,
                "y_development": y,
                "residual": residualize_demographics(y, demographics, dev_indices),
                "screen_score_max_abs_r": float(np.max(np.abs(r))),
                "cis_qc_variant_count": int(len(cis)),
                "variant_indices": chosen,
            }
        )
    candidates.sort(key=lambda x: (-float(x["screen_score_max_abs_r"]), str(x["gene_id"])))
    selected: list[dict[str, object]] = []
    for candidate in candidates:
        if any(abs(int(candidate["tss_1based"]) - int(other["tss_1based"])) < int(config["min_tss_separation_bp"]) for other in selected):
            continue
        if any(abs(float(np.corrcoef(candidate["residual"], other["residual"])[0, 1])) > float(config["max_abs_development_residual_trait_correlation"]) for other in selected):
            continue
        selected.append(candidate)
        if len(selected) == int(config["trait_count"]):
            break
    if len(selected) != int(config["trait_count"]):
        raise RuntimeError(f"DESIGN-INELIGIBLE: only {len(selected)} traits satisfy the frozen diversity rules")
    return selected


def design_sensitivity(n: int, k: int, rho: float, config: dict[str, object]) -> dict[str, object]:
    rng = np.random.default_rng(int(config["seed"]) + 700)
    simulations = int(config["design_simulations"])
    covariance = (1.0 - rho) * np.eye(k) + rho * np.ones((k, k))
    chol = np.linalg.cholesky(covariance)
    values: list[np.ndarray] = []
    batch = 200
    for _ in range(math.ceil(simulations / batch)):
        m = min(batch, simulations - sum(len(x) for x in values))
        if m <= 0:
            break
        draw = lambda: rng.normal(size=(m, n, k)) @ chol.T
        baseline = math.sqrt(float(config["design_assumed_B_r2"])) * draw()
        increment = math.sqrt(float(config["minimum_scientific_bridge_delta_r2"])) * draw()
        noise_var = 1.0 - float(config["design_assumed_B_r2"]) - float(config["minimum_scientific_bridge_delta_r2"])
        outcome = baseline + increment + math.sqrt(noise_var) * draw()
        denom = np.sum((outcome - outcome.mean(axis=1, keepdims=True)) ** 2, axis=1)
        r2_b = 1.0 - np.sum((outcome - baseline) ** 2, axis=1) / denom
        r2_c = 1.0 - np.sum((outcome - baseline - increment) ** 2, axis=1) / denom
        values.append(np.mean(r2_c - r2_b, axis=1))
    delta = np.concatenate(values)
    sd = float(np.std(delta, ddof=1))
    critical = 1.959963984540054 * sd
    power = float(np.mean(delta > critical))
    return {
        "method": "frozen_gaussian_monte_carlo_normal_critical_value",
        "n_bridge": n,
        "traits": k,
        "exchangeable_trait_correlation": rho,
        "simulations": simulations,
        "assumed_B_r2": float(config["design_assumed_B_r2"]),
        "minimum_scientific_delta_r2": float(config["minimum_scientific_bridge_delta_r2"]),
        "mean_simulated_delta_r2": float(np.mean(delta)),
        "sd_simulated_estimator": sd,
        "one_sided_critical_delta": critical,
        "estimated_power": power,
        "required_power": float(config["design_power_min"]),
        "eligible": bool(power >= float(config["design_power_min"])),
        "limitation": "model-based planning sensitivity, not an empirical guarantee",
    }


def fixed_covariates(
    rows: list[dict[str, str]], dosage: np.ndarray, train_idx: np.ndarray, config: dict[str, object]
) -> tuple[np.ndarray, dict[str, object]]:
    demographics = base_demographics(rows, np.arange(len(rows)), train_idx)
    mean = dosage[train_idx].mean(axis=0, dtype=np.float64)
    sd = dosage[train_idx].std(axis=0, dtype=np.float64)
    sd[sd < 1e-8] = 1.0
    standardized = (dosage.astype(np.float64) - mean) / sd
    n_components = min(int(config["genotype_pc_count"]), len(train_idx) - 1, dosage.shape[1])
    pca = PCA(n_components=n_components, svd_solver="full")
    pca.fit(standardized[train_idx])
    pcs = pca.transform(standardized)
    return np.column_stack([demographics, pcs]).astype(np.float32), {
        "sex_columns": 1,
        "population_columns": ["FIN", "GBR", "TSI", "YRI"],
        "genotype_pc_count": n_components,
        "pc_fit_scope": "dev_train_only",
        "pc_explained_variance_ratio": [float(x) for x in pca.explained_variance_ratio_],
    }


def train_encoder(panel: np.ndarray, af: np.ndarray, train_idx: np.ndarray, val_idx: np.ndarray, config: dict[str, object], trait_offset: int):
    cfg = ModelConfig(
        max_blocks=panel.shape[1] // int(config["block_size"]),
        snps_per_block=int(config["block_size"]),
        d_model=int(config["d_model"]),
        n_heads=int(config["n_heads"]),
        local_layers=int(config["local_layers"]),
        local_ff_dim=int(config["local_ff_dim"]),
        spectral_rank=int(config["spectral_rank"]),
        spectral_oversample=4,
        dropout=float(config["dropout"]),
    )
    model = LocalGenotypeEncoder(cfg)
    optimizer = torch.optim.AdamW(model.parameters(), lr=float(config["learning_rate"]), weight_decay=float(config["weight_decay"]))
    n_blocks = cfg.max_blocks
    geno = torch.from_numpy(panel.reshape(len(panel), n_blocks, cfg.snps_per_block).astype(np.int64))
    af_tensor = torch.from_numpy(af.reshape(n_blocks, cfg.snps_per_block).astype(np.float32))
    mask_rng = torch.Generator(device="cpu").manual_seed(int(config["seed"]) + 1000 + trait_offset)
    order_rng = np.random.default_rng(int(config["seed"]) + 2000 + trait_offset)
    val_rng = torch.Generator(device="cpu").manual_seed(int(config["seed"]) + 3000 + trait_offset)
    val_mask = sample_ssl_mask(geno[val_idx], float(config["mask_probability"]), val_rng)
    history = []
    for epoch in range(1, int(config["pretrain_epochs"]) + 1):
        model.train()
        losses = []
        order = order_rng.permutation(train_idx)
        for batch_indices in np.array_split(order, math.ceil(len(order) / int(config["batch_size"]))):
            batch = geno[batch_indices]
            mask = sample_ssl_mask(batch, float(config["mask_probability"]), mask_rng)
            _, logits = model(encode_genotypes(batch, mask), af_tensor)
            loss = masked_genotype_loss(logits, batch, mask)
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            optimizer.step()
            losses.append(float(loss.detach()))
        model.eval()
        with torch.no_grad():
            _, logits = model(encode_genotypes(geno[val_idx], val_mask), af_tensor)
            val_loss = float(masked_genotype_loss(logits, geno[val_idx], val_mask))
        history.append({"epoch": epoch, "train_loss": float(np.mean(losses)), "validation_loss": val_loss})
    model.eval()
    with torch.no_grad():
        _, logits = model(encode_genotypes(geno[val_idx], val_mask), af_tensor)
        validation_metrics = contextual_gate_metrics(logits, geno[val_idx], val_mask, af_tensor)
        validation_metrics["masked_cross_entropy"] = float(masked_genotype_loss(logits, geno[val_idx], val_mask))
    return model, cfg, history, validation_metrics


def extract_hidden(model: LocalGenotypeEncoder, panel: np.ndarray, af: np.ndarray, block_size: int) -> np.ndarray:
    n_blocks = panel.shape[1] // block_size
    geno = torch.from_numpy(panel.reshape(len(panel), n_blocks, block_size).astype(np.int64))
    af_tensor = torch.from_numpy(af.reshape(n_blocks, block_size).astype(np.float32))
    chunks = []
    with torch.no_grad():
        for start in range(0, len(panel), 32):
            hidden, _ = model(encode_genotypes(geno[start : start + 32]), af_tensor)
            chunks.append(hidden.cpu().numpy().astype(np.float32))
    return np.concatenate(chunks, axis=0)


def aligned_features(hidden: np.ndarray, train_idx: np.ndarray, rank: int) -> tuple[np.ndarray, list[float]]:
    tokens = hidden[train_idx].reshape(-1, hidden.shape[-1]).astype(np.float64)
    mean = tokens.mean(axis=0)
    covariance = np.cov(tokens - mean, rowvar=False, ddof=1)
    values, vectors = np.linalg.eigh(covariance)
    order = np.argsort(values)[::-1][:rank]
    scores = np.einsum("nbsd,dr->nbsr", hidden.astype(np.float64) - mean, vectors[:, order])
    return scores.reshape(len(hidden), -1).astype(np.float32), [float(values[i]) for i in order]


def grassmann_features(hidden: np.ndarray, rank: int) -> np.ndarray:
    centered = hidden.astype(np.float64) - hidden.astype(np.float64).mean(axis=2, keepdims=True)
    covariance = np.einsum("nbsd,nbse->nbde", centered, centered) / float(hidden.shape[2] - 1)
    values, vectors = np.linalg.eigh(covariance)
    values = values[..., -rank:][..., ::-1].copy()
    vectors = vectors[..., -rank:][..., ::-1].copy()
    projector = np.einsum("nbdr,nber->nbde", vectors, vectors)
    triangle = np.triu_indices(hidden.shape[-1])
    p = projector[..., triangle[0], triangle[1]]
    spectrum = np.log1p(np.clip(values, 0.0, None))
    spectrum /= np.maximum(spectrum.sum(axis=-1, keepdims=True), 1e-12)
    return np.concatenate([p, spectrum], axis=-1).reshape(len(hidden), -1).astype(np.float32)


def fit_ridge_predictions(x: np.ndarray, y_dev: np.ndarray, train_idx: np.ndarray, val_idx: np.ndarray, predict_idx: np.ndarray, alphas: list[float]):
    x_train = x[train_idx].astype(np.float64)
    mean = x_train.mean(axis=0)
    sd = x_train.std(axis=0)
    keep = sd >= 1e-8
    z_train = (x_train[:, keep] - mean[keep]) / sd[keep]
    z_val = (x[val_idx][:, keep].astype(np.float64) - mean[keep]) / sd[keep]
    z_predict = (x[predict_idx][:, keep].astype(np.float64) - mean[keep]) / sd[keep]
    local_train = np.searchsorted(np.sort(np.concatenate([train_idx, val_idx])), train_idx)
    del local_train
    y_train = y_dev[: len(train_idx)]
    y_val = y_dev[len(train_idx) :]
    y_mean = float(y_train.mean())
    y_sd = float(y_train.std(ddof=0))
    y_train_std = (y_train - y_mean) / y_sd
    gram = z_train @ z_train.T
    val_cross = z_val @ z_train.T
    predict_cross = z_predict @ z_train.T
    identity = np.eye(len(train_idx))
    best = None
    for alpha in alphas:
        dual = np.linalg.solve(gram + float(alpha) * identity, y_train_std)
        val_prediction = y_mean + y_sd * (val_cross @ dual)
        mse = float(np.mean((y_val - val_prediction) ** 2))
        if best is None or mse < best[0] - 1e-12 or (abs(mse - best[0]) <= 1e-12 and alpha > best[1]):
            best = (mse, float(alpha), dual)
    assert best is not None
    return y_mean + y_sd * (predict_cross @ best[2]), {"alpha": best[1], "validation_mse": best[0], "input_features": int(x.shape[1]), "nonconstant_features": int(keep.sum()), "train_mean": y_mean}


def summarize_gate(y: list[np.ndarray], predictions: dict[str, list[np.ndarray]], train_means: list[float], contrast: tuple[str, str], seed: int, replicates: int):
    n = len(y[0])
    if any(len(v) != n for v in y):
        raise ValueError("trait holdout sizes differ")
    arms = list(predictions)
    trait_metrics: list[dict[str, object]] = []
    point = {arm: [] for arm in arms}
    for k, truth in enumerate(y):
        denominator = float(np.sum((truth - train_means[k]) ** 2))
        if denominator <= 0:
            raise ValueError("nonpositive R2 denominator")
        row = {"trait_index": k}
        for arm in arms:
            value = 1.0 - float(np.sum((truth - predictions[arm][k]) ** 2)) / denominator
            point[arm].append(value)
            row[f"r2_{arm}"] = value
        trait_metrics.append(row)
    macro = {arm: float(np.mean(point[arm])) for arm in arms}
    rng = np.random.default_rng(seed)
    boot = {arm: np.zeros(replicates) for arm in arms}
    for start in range(0, replicates, 250):
        stop = min(replicates, start + 250)
        sample = rng.integers(0, n, size=(stop - start, n))
        for arm in arms:
            across_traits = []
            for k, truth in enumerate(y):
                yb = truth[sample]
                pb = predictions[arm][k][sample]
                denom = np.sum((yb - train_means[k]) ** 2, axis=1)
                across_traits.append(1.0 - np.sum((yb - pb) ** 2, axis=1) / denom)
            boot[arm][start:stop] = np.mean(np.stack(across_traits, axis=1), axis=1)
    left, right = contrast
    delta_boot = boot[left] - boot[right]
    delta = macro[left] - macro[right]
    return {
        "n_individuals": n,
        "macro_r2": macro,
        "trait_metrics": trait_metrics,
        "contrast": f"{left}-{right}",
        "macro_delta_r2": delta,
        "paired_bootstrap_ci95": [float(np.quantile(delta_boot, 0.025)), float(np.quantile(delta_boot, 0.975))],
        "pass": bool(delta > 0 and np.quantile(delta_boot, 0.025) > 0),
        "bootstrap_replicates": replicates,
        "bootstrap_seed": seed,
    }


def main() -> None:
    if (RESULT_DIR / "FINAL_STATUS.json").exists():
        raise RuntimeError("frozen Stage 2B already has a final status and cannot be rerun")
    config = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
    set_reproducibility(int(config["seed"]))
    rows = load_tsv(MANIFEST)
    psam_ids = [row["#IID"] for row in load_tsv(PSAM)]
    if [row["sample_id"] for row in rows] != psam_ids or expression_header() != psam_ids:
        raise ValueError("PSAM/expression/manifest sample order mismatch")
    role_indices = {
        role: np.asarray([i for i, row in enumerate(rows) if row["primary_role"] == role], dtype=np.int64)
        for role in ["development", "task_gate", "bridge_test", "EXCLUDED_RC2_TEST"]
    }
    train_idx = np.asarray([i for i, row in enumerate(rows) if row["development_role"] == "dev_train"], dtype=np.int64)
    val_idx = np.asarray([i for i, row in enumerate(rows) if row["development_role"] == "dev_validation"], dtype=np.int64)
    dev_idx = np.concatenate([train_idx, val_idx])
    if not np.array_equal(np.sort(dev_idx), np.sort(role_indices["development"])):
        raise ValueError("development role mismatch")
    # Ordering y_development as train then validation makes decoder scope explicit.
    outcome_dev_order = np.concatenate([train_idx, val_idx])
    demographics = base_demographics(rows, np.arange(len(rows)), train_idx)
    positions, variant_ids, refs, alts = read_pvar()
    genotype = read_pgen(len(rows), len(positions))
    print(f"asset audit passed: samples={len(rows)} variants={len(positions)}", flush=True)

    # Development phenotype is the only outcome opened before design eligibility and model freeze.
    development_records = load_expression_subset([psam_ids[i] for i in outcome_dev_order])
    selected = select_traits_and_panels(development_records, genotype, positions, variant_ids, outcome_dev_order, train_idx, demographics, config)
    observed_corr = max(abs(float(np.corrcoef(a["residual"], b["residual"])[0, 1])) for i, a in enumerate(selected) for b in selected[i + 1 :])
    planning_rho = max(float(config["design_trait_correlation_floor"]), observed_corr)
    sensitivity = design_sensitivity(len(role_indices["bridge_test"]), len(selected), planning_rho, config)
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    dump_json(RESULT_DIR / "DESIGN_SENSITIVITY.json", sensitivity)

    trait_manifest = []
    for item in selected:
        idx = np.asarray(item["variant_indices"], dtype=np.int64)
        trait_manifest.append({
            "gene_id": item["gene_id"],
            "tss_1based": item["tss_1based"],
            "screen_score_max_abs_r": item["screen_score_max_abs_r"],
            "cis_qc_variant_count": item["cis_qc_variant_count"],
            "variant_indices": [int(x) for x in idx],
            "variant_ids": [variant_ids[int(x)] for x in idx],
            "positions_1based": [int(positions[int(x)]) for x in idx],
            "refs": [refs[int(x)] for x in idx],
            "alts": [alts[int(x)] for x in idx],
        })
    dump_json(RESULT_DIR / "TRAIT_PANEL.json", trait_manifest)
    if not sensitivity["eligible"]:
        final = {"analysis_id": config["analysis_id"], "status": "DESIGN-INELIGIBLE", "design_sensitivity": sensitivity, "task_gate_opened": False, "bridge_gate_opened": False}
        dump_json(RESULT_DIR / "FINAL_STATUS.json", final)
        print(json.dumps(final, indent=2), flush=True)
        return

    task_idx = role_indices["task_gate"]
    bridge_idx = role_indices["bridge_test"]
    predictions_task = {arm: [] for arm in "ABCDE"}
    predictions_bridge = {arm: [] for arm in "ABCDE"}
    train_means: list[float] = []
    model_records = []
    checkpoint_dir = RESULT_DIR / "checkpoints"
    checkpoint_dir.mkdir(exist_ok=True)
    for trait_number, item in enumerate(selected):
        gene = str(item["gene_id"])
        snp_idx = np.asarray(item["variant_indices"], dtype=np.int64)
        raw_panel = genotype[:, snp_idx]
        af, _, dosage = af_and_dosage(raw_panel, train_idx)
        covariates, covariate_info = fixed_covariates(rows, dosage, train_idx, config)
        model, model_cfg, history, pretrain_metrics = train_encoder(raw_panel, af, train_idx, val_idx, config, trait_number)
        checkpoint = checkpoint_dir / f"{gene}.pt"
        torch.save({"state_dict": model.state_dict(), "config": asdict(model_cfg), "gene_id": gene, "variant_ids": trait_manifest[trait_number]["variant_ids"], "phenotype_seen_by_encoder": False}, checkpoint)
        for parameter in model.parameters():
            parameter.requires_grad_(False)
        hidden = extract_hidden(model, raw_panel, af, int(config["block_size"]))
        aligned, eigenvalues = aligned_features(hidden, train_idx, int(config["spectral_rank"]))
        grassmann = grassmann_features(hidden, int(config["spectral_rank"]))
        features = {
            "A": covariates,
            "B": np.concatenate([covariates, dosage], axis=1),
            "C": np.concatenate([covariates, hidden.reshape(len(rows), -1)], axis=1),
            "D": np.concatenate([covariates, aligned], axis=1),
            "E": np.concatenate([covariates, grassmann], axis=1),
        }
        y_dev = np.asarray(item["y_development"], dtype=np.float64)
        decoders = {}
        for arm in "ABCDE":
            pred_all, decoder = fit_ridge_predictions(features[arm], y_dev, train_idx, val_idx, np.concatenate([task_idx, bridge_idx]), list(config["ridge_alphas"]))
            predictions_task[arm].append(pred_all[: len(task_idx)])
            predictions_bridge[arm].append(pred_all[len(task_idx) :])
            decoders[arm] = decoder
        train_means.append(float(y_dev[: len(train_idx)].mean()))
        model_records.append({
            "gene_id": gene,
            "checkpoint_sha256": sha256_file(checkpoint),
            "pretrain_validation": pretrain_metrics,
            "pretrain_history": history,
            "covariates": covariate_info,
            "aligned_eigenvalues": eigenvalues,
            "feature_dimensions": {arm: int(features[arm].shape[1]) for arm in "ABCDE"},
            "decoders": decoders,
        })
        print(f"frozen model {trait_number + 1}/8 gene={gene} contextual_lift={pretrain_metrics['contextual_lift']:.4f}", flush=True)
    dump_json(RESULT_DIR / "FROZEN_MODELS.json", model_records)

    # Only now open Task Gate outcomes; Bridge outcomes are still not read.
    selected_genes = {str(x["gene_id"]) for x in selected}
    task_records = {str(x["gene_id"]): x for x in load_expression_subset([psam_ids[i] for i in task_idx], selected_genes)}
    y_task = [np.asarray(task_records[str(item["gene_id"])]["values"], dtype=np.float64) for item in selected]
    task = summarize_gate(y_task, {arm: predictions_task[arm] for arm in "AB"}, train_means, ("B", "A"), int(config["bootstrap_seed"]), int(config["bootstrap_replicates"]))
    task["status"] = "PASS" if task["pass"] else "TASK-INELIGIBLE"
    dump_json(RESULT_DIR / "TASK_GATE_RESULTS.json", task)
    print(f"TASK {task['status']} macro_B-A={task['macro_delta_r2']:.6f} CI={task['paired_bootstrap_ci95']}", flush=True)
    if not task["pass"]:
        final = {"analysis_id": config["analysis_id"], "status": "TASK-INELIGIBLE", "task_gate": task, "bridge_gate_opened": False, "test_open_count": 0}
        dump_json(RESULT_DIR / "FINAL_STATUS.json", final)
        write_report(final, selected, sensitivity)
        write_manifest()
        return

    # Authorized single opening of Bridge outcomes.
    bridge_records = {str(x["gene_id"]): x for x in load_expression_subset([psam_ids[i] for i in bridge_idx], selected_genes)}
    y_bridge = [np.asarray(bridge_records[str(item["gene_id"])]["values"], dtype=np.float64) for item in selected]
    bridge = summarize_gate(y_bridge, predictions_bridge, train_means, ("C", "B"), int(config["bootstrap_seed"]) + 1, int(config["bootstrap_replicates"]))
    bridge["secondary_contrasts"] = {}
    for left, right in [("D", "B"), ("E", "B"), ("E", "D")]:
        secondary = summarize_gate(y_bridge, {left: predictions_bridge[left], right: predictions_bridge[right]}, train_means, (left, right), int(config["bootstrap_seed"]) + 1, int(config["bootstrap_replicates"]))
        bridge["secondary_contrasts"][f"{left}-{right}"] = {"macro_delta_r2": secondary["macro_delta_r2"], "paired_bootstrap_ci95": secondary["paired_bootstrap_ci95"]}
    bridge["status"] = "PASS" if bridge["pass"] else "NO-GO"
    dump_json(RESULT_DIR / "BRIDGE_GATE_RESULTS.json", bridge)
    final = {
        "analysis_id": config["analysis_id"],
        "status": "BRIDGE-PASS" if bridge["pass"] else "BRIDGE-NO-GO",
        "task_gate": task,
        "bridge_gate": bridge,
        "bridge_gate_opened": True,
        "test_open_count": 1,
        "claim_stopped_if_no_go": "generic genotype-pretrained full-H exceeds same-panel dosage in this high-signal GEUVADIS chr18 replication setting",
        "claims_not_stopped": ["SNP foundation modelling", "learning LD/MAF/haplotype", "independently justified phenotype models", "all possible non-H representations"],
    }
    dump_json(RESULT_DIR / "FINAL_STATUS.json", final)
    write_report(final, selected, sensitivity)
    write_manifest()
    print(f"BRIDGE {bridge['status']} macro_C-B={bridge['macro_delta_r2']:.6f} CI={bridge['paired_bootstrap_ci95']}", flush=True)


def write_report(final: dict[str, object], selected: list[dict[str, object]], sensitivity: dict[str, object]) -> None:
    task = final.get("task_gate", {})
    bridge = final.get("bridge_gate")
    genes = ", ".join(str(x["gene_id"]) for x in selected)
    bridge_text = "Bridge phenotype 未打开。"
    if isinstance(bridge, dict):
        bridge_text = f"Bridge `C-B` = **{bridge['macro_delta_r2']:.6f}**, paired 95% CI = **[{bridge['paired_bootstrap_ci95'][0]:.6f}, {bridge['paired_bootstrap_ci95'][1]:.6f}]**；判定 **{bridge['status']}**。"
    text = f"""# 0526 Stage 2B GEUVADIS chr18 rc1 结果

## 冻结结论

最终状态：**{final['status']}**。

Task `B-A` = **{task.get('macro_delta_r2', float('nan')):.6f}**，paired 95% CI = **{task.get('paired_bootstrap_ci95')}**；判定 **{task.get('status')}**。

{bridge_text}

## 审计摘要

- Development / Task / Bridge / RC2 excluded = 150 / 50 / 230 / 15；家庭零交叉。
- traits：{genes}。
- 设计灵敏度：预设 ΔR²=0.02，Monte Carlo power={sensitivity['estimated_power']:.4f}，要求 ≥0.80。
- A-E 使用每个 trait 完全相同的 256 个 cis SNP、样本、covariates、decoder family 和一次性 test opening。
- RC2 原结论保持：`0526 Stage 2 bridge, RC2: NO-GO in this setting.`

## 解释边界

本结果只裁决本冻结 setting 的 full-H bridge。它不把 MAF/LD/haplotype 学习判为失败，也不裁决整个 SNP foundation-model programme。发布的 normalized expression 在分组前已由上游处理，因此这不是 raw-expression preprocessing-sealed replication。
"""
    (RESULT_DIR / "REPORT.zh-CN.md").write_text(text, encoding="utf-8")


def write_manifest() -> None:
    binding = {
        "assets": {path.name: {"bytes": path.stat().st_size, "sha256": sha256_file(path)} for path in [PGEN, PVAR, PSAM, EXPRESSION, SAMPLE_PANEL, PEDIGREE]},
        "runtime": {"python": platform.python_version(), "platform": platform.platform(), "numpy": np.__version__, "sklearn": sklearn.__version__, "torch": torch.__version__, "pgenlib": getattr(pgenlib, "__version__", "0.94.1-installed"), "cuda": torch.cuda.is_available()},
        "protocol_sha256": sha256_file(PROTOCOL_PATH),
        "config_sha256": sha256_file(CONFIG_PATH),
        "sample_manifest_sha256": sha256_file(MANIFEST),
    }
    dump_json(RESULT_DIR / "RUN_BINDING.json", binding)
    paths = [Path(__file__), PILOT_DIR / "freeze_design.py", PROTOCOL_PATH, CONFIG_PATH, DESIGN_PATH, MANIFEST]
    paths += sorted(path for path in RESULT_DIR.rglob("*") if path.is_file() and path.name != "RUN_MANIFEST.sha256")
    lines = [f"{sha256_file(path)}  {path.relative_to(ROOT).as_posix()}" for path in paths]
    (RESULT_DIR / "RUN_MANIFEST.sha256").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
