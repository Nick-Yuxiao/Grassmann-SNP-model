"""TQ-G1: annotation-free gene-layer readout Pilot on GEUVADIS chr18.

Primary question: does pooling a frozen shared genotype encoder's hidden state
into a *gene* unit add held-out predictive R2 over the same-window additive
dosage ridge, when the gene is the unit of replication?

This runner never opens a sealed role, never uses functional annotation, and
never selects genes or cis SNPs on phenotype signal. Run it once:

    python run_tqg1.py

Artifacts land in ``results/``. A completed FINAL_STATUS.json blocks a rerun.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import platform
import random
import sys
from pathlib import Path

import numpy as np
import sklearn
import torch

PILOT_DIR = Path(__file__).resolve().parent
PACKAGE = PILOT_DIR.parents[1]
sys.path.insert(0, str(PACKAGE / "src"))
sys.path.insert(0, str(PILOT_DIR))

from snp_spectral_foundation.config import ModelConfig  # noqa: E402
from snp_spectral_foundation.model import LocalGenotypeEncoder, encode_genotypes  # noqa: E402
from snp_spectral_foundation.training import (  # noqa: E402
    contextual_gate_metrics,
    masked_genotype_loss,
    sample_ssl_mask,
)
from tqg1_core import (  # noqa: E402
    af_missingness_dosage,
    decile_strata,
    deterministic_gene_order,
    fit_ridge,
    gene_pool_features,
    matched_null_spearman,
    paired_gene_bootstrap_difference,
    r2_against_train_mean,
    spearman,
    two_way_clustered_bootstrap,
    uniform_thin,
)

CONFIG_PATH = PILOT_DIR / "CONFIG_FROZEN.json"
BINDING_PATH = PILOT_DIR / "PRE_RUN_BINDING.json"
PROTOCOL_PATH = PILOT_DIR / "FROZEN_PROTOCOL.zh-CN.md"
RESULT_DIR = PILOT_DIR / "results"
AUTHORIZATION_PATH = PILOT_DIR / "OPEN_AUTHORIZATION.json"

MASK_EVERYTHING = True


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
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


# ---------------------------------------------------------------------------
# asset binding
# ---------------------------------------------------------------------------


def resolve_assets(binding: dict[str, object]) -> dict[str, Path]:
    raw = binding.get("asset_dir")
    if not raw:
        raise SystemExit(
            "PRE_RUN_BINDING.json: asset_dir is unset. Point it at the directory holding the "
            "six chr18 GEUVADIS assets, then rerun."
        )
    asset_dir = Path(str(raw)).expanduser()
    if not asset_dir.is_absolute():
        asset_dir = (PILOT_DIR / asset_dir).resolve()
    if not asset_dir.is_dir():
        raise SystemExit(f"asset_dir does not exist: {asset_dir}")
    names = binding["assets"]
    assert isinstance(names, dict)
    paths = {key: asset_dir / str(name) for key, name in names.items()}
    missing = [str(p) for p in paths.values() if not p.is_file()]
    if missing:
        raise SystemExit("missing assets:\n  " + "\n  ".join(missing))
    return paths


def audit_assets(paths: dict[str, Path], binding: dict[str, object]) -> dict[str, dict[str, object]]:
    expected = binding.get("expected_sha256") or {}
    assert isinstance(expected, dict)
    audit: dict[str, dict[str, object]] = {}
    for path in paths.values():
        digest = sha256_file(path)
        want = expected.get(path.name)
        if want and want != digest:
            raise SystemExit(
                f"asset hash mismatch for {path.name}\n  expected {want}\n  observed {digest}"
            )
        audit[path.name] = {"bytes": path.stat().st_size, "sha256": digest, "hash_enforced": bool(want)}
    return audit


# ---------------------------------------------------------------------------
# readers
# ---------------------------------------------------------------------------


def read_pvar(path: Path, chromosome: str, expected_count: int) -> tuple[np.ndarray, list[str]]:
    positions: list[int] = []
    ids: list[str] = []
    header: list[str] | None = None
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line.rstrip("\n").lstrip("#").split("\t")
                if header[:5] != ["CHROM", "POS", "ID", "REF", "ALT"]:
                    raise ValueError("unexpected PVAR header")
                continue
            row = line.rstrip("\n").split("\t")
            if header is None or row[0] != chromosome:
                raise ValueError("PVAR coordinate/contig contract violation")
            pos = int(row[1])
            if positions and pos < positions[-1]:
                raise ValueError("PVAR is not coordinate sorted")
            if "," in row[4]:
                raise ValueError("multiallelic PVAR record")
            positions.append(pos)
            ids.append(row[2])
    if len(positions) != expected_count or len(set(ids)) != len(ids):
        raise ValueError(f"unexpected PVAR count {len(positions)} or duplicate IDs")
    return np.asarray(positions, dtype=np.int64), ids


def read_pgen(path: Path, n_samples: int, n_variants: int) -> np.ndarray:
    import pgenlib

    genotype = np.empty((n_samples, n_variants), dtype=np.int8)
    with pgenlib.PgenReader(str(path).encode()) as reader:
        if reader.get_raw_sample_ct() != n_samples or reader.get_variant_ct() != n_variants:
            raise ValueError("PGEN/PVAR/PSAM dimension mismatch")
        reader.read_range(0, n_variants, genotype, sample_maj=1)
    if np.any((genotype < -9) | (genotype > 2)):
        raise ValueError("invalid hardcall domain")
    # The encoder's token contract is {-1, 0, 1, 2}; PGEN emits -9 for missing.
    genotype[genotype < 0] = -1
    return genotype


def expression_header(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        header = next(csv.reader(handle, delimiter="\t"))
    if header[:4] != ["#chr", "start", "end", "gene_id"]:
        raise ValueError("unexpected expression BED header")
    return header[4:]


def load_expression(path: Path, chromosome: str, sample_ids: list[str]) -> list[dict[str, object]]:
    """Read expression for the named samples only.

    Sealed roles are never passed in, so their outcomes are never in memory.
    """
    all_ids = expression_header(path)
    lookup = {sample: i + 4 for i, sample in enumerate(all_ids)}
    columns = [lookup[sample] for sample in sample_ids]
    records: list[dict[str, object]] = []
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader)
        for row in reader:
            if row[0] != chromosome:
                continue
            start, end = int(row[1]), int(row[2])
            if end != start + 1:
                raise ValueError(f"expression row is not a one-base TSS interval: {row[3]}")
            records.append(
                {
                    "gene_id": row[3],
                    "tss_1based": start + 1,
                    "values": np.asarray([float(row[i]) for i in columns], dtype=np.float64),
                }
            )
    return records


# ---------------------------------------------------------------------------
# covariates
# ---------------------------------------------------------------------------


POPULATIONS = ["CEU", "FIN", "GBR", "TSI", "YRI"]


def demographic_covariates(rows: list[dict[str, str]]) -> np.ndarray:
    sex = np.asarray([1.0 if row["sex"] == "2" else 0.0 for row in rows], dtype=np.float64)
    pop = np.asarray(
        [[float(row["pop"] == p) for p in POPULATIONS[1:]] for row in rows], dtype=np.float64
    )
    return np.column_stack([sex, pop]).astype(np.float32)


def genotype_pcs(
    dosage_panel: np.ndarray, train_idx: np.ndarray, n_components: int
) -> tuple[np.ndarray, list[float]]:
    from sklearn.decomposition import PCA

    mean = dosage_panel[train_idx].mean(axis=0, dtype=np.float64)
    sd = dosage_panel[train_idx].std(axis=0, dtype=np.float64)
    sd[sd < 1e-8] = 1.0
    standardized = (dosage_panel.astype(np.float64) - mean) / sd
    k = min(n_components, len(train_idx) - 1, dosage_panel.shape[1])
    pca = PCA(n_components=k, svd_solver="full")
    pca.fit(standardized[train_idx])
    return pca.transform(standardized).astype(np.float32), [
        float(x) for x in pca.explained_variance_ratio_
    ]


# ---------------------------------------------------------------------------
# shared encoder
# ---------------------------------------------------------------------------


def build_config(config: dict[str, object], n_blocks: int) -> ModelConfig:
    return ModelConfig(
        max_blocks=max(1, n_blocks),
        snps_per_block=int(config["block_size"]),
        d_model=int(config["d_model"]),
        n_heads=int(config["n_heads"]),
        local_layers=int(config["local_layers"]),
        local_ff_dim=int(config["local_ff_dim"]),
        dropout=float(config["dropout"]),
        local_position_mode=str(config["local_position_mode"]),
        use_local_block_embedding=bool(config["use_local_block_embedding"]),
    )


def pretrain_shared_encoder(
    genotype_panel: np.ndarray,
    af: np.ndarray,
    positions: np.ndarray,
    train_idx: np.ndarray,
    val_idx: np.ndarray,
    config: dict[str, object],
) -> tuple[LocalGenotypeEncoder, ModelConfig, list[dict[str, float]], dict[str, float]]:
    block_size = int(config["block_size"])
    n_blocks = genotype_panel.shape[1] // block_size
    cfg = build_config(config, n_blocks)
    model = LocalGenotypeEncoder(cfg)
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=float(config["learning_rate"]),
        weight_decay=float(config["weight_decay"]),
    )
    geno = torch.from_numpy(genotype_panel.reshape(-1, n_blocks, block_size).astype(np.int64))
    af_tensor = torch.from_numpy(af.reshape(n_blocks, block_size).astype(np.float32))
    pos_tensor = torch.from_numpy(positions.reshape(n_blocks, block_size).astype(np.int64))

    mask_rng = torch.Generator(device="cpu").manual_seed(int(config["seed"]) + 1000)
    order_rng = np.random.default_rng(int(config["seed"]) + 2000)
    val_rng = torch.Generator(device="cpu").manual_seed(int(config["seed"]) + 3000)
    val_mask = sample_ssl_mask(geno[val_idx], float(config["mask_probability"]), val_rng)

    history: list[dict[str, float]] = []
    batch_size = int(config["batch_size"])
    for epoch in range(1, int(config["pretrain_epochs"]) + 1):
        model.train()
        losses: list[float] = []
        order = order_rng.permutation(train_idx)
        for chunk in np.array_split(order, math.ceil(len(order) / batch_size)):
            batch = geno[chunk]
            mask = sample_ssl_mask(batch, float(config["mask_probability"]), mask_rng)
            _, logits = model(
                encode_genotypes(batch, mask), af_tensor, genomic_positions=pos_tensor
            )
            loss = masked_genotype_loss(logits, batch, mask)
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            optimizer.step()
            losses.append(float(loss.detach()))
        model.eval()
        with torch.no_grad():
            _, logits = model(
                encode_genotypes(geno[val_idx], val_mask), af_tensor, genomic_positions=pos_tensor
            )
            val_loss = float(masked_genotype_loss(logits, geno[val_idx], val_mask))
        history.append({"epoch": epoch, "train_loss": float(np.mean(losses)), "validation_loss": val_loss})
        print(f"  epoch {epoch:3d}  train {np.mean(losses):.4f}  validation {val_loss:.4f}", flush=True)

    model.eval()
    with torch.no_grad():
        _, logits = model(
            encode_genotypes(geno[val_idx], val_mask), af_tensor, genomic_positions=pos_tensor
        )
        metrics = contextual_gate_metrics(logits, geno[val_idx], val_mask, af_tensor)
        metrics["masked_cross_entropy"] = float(masked_genotype_loss(logits, geno[val_idx], val_mask))
    return model, cfg, history, metrics


@torch.no_grad()
def encode_window(
    model: LocalGenotypeEncoder,
    panel: np.ndarray,
    af: np.ndarray,
    positions: np.ndarray,
    block_size: int,
    mask_all: bool,
    chunk: int = 64,
) -> np.ndarray:
    """Hidden states for one cis window; ``mask_all`` is the gene perturbation."""
    n_blocks = panel.shape[1] // block_size
    geno = torch.from_numpy(panel.reshape(-1, n_blocks, block_size).astype(np.int64))
    af_tensor = torch.from_numpy(af.reshape(n_blocks, block_size).astype(np.float32))
    pos_tensor = torch.from_numpy(positions.reshape(n_blocks, block_size).astype(np.int64))
    out: list[np.ndarray] = []
    for start in range(0, len(geno), chunk):
        batch = geno[start : start + chunk]
        masked = torch.ones_like(batch, dtype=torch.bool) if mask_all else None
        hidden, _ = model(
            encode_genotypes(batch, masked), af_tensor, genomic_positions=pos_tensor
        )
        out.append(hidden.cpu().numpy().astype(np.float32))
    return np.concatenate(out, axis=0)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def resolve_roles(
    rows: list[dict[str, str]], evaluation_role: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Bind dev_train, dev_validation and the evaluation role, checking leakage."""
    train_idx = np.asarray(
        [i for i, r in enumerate(rows) if r["development_role"] == "dev_train"], dtype=np.int64
    )
    val_idx = np.asarray(
        [i for i, r in enumerate(rows) if r["development_role"] == "dev_validation"], dtype=np.int64
    )
    eval_idx = np.asarray(
        [i for i, r in enumerate(rows) if r["primary_role"] == evaluation_role], dtype=np.int64
    )
    if len(train_idx) == 0 or len(val_idx) == 0 or len(eval_idx) == 0:
        raise SystemExit("empty dev_train, dev_validation or evaluation role")
    overlap = set(eval_idx) & (set(train_idx) | set(val_idx))
    if overlap:
        raise SystemExit("evaluation individuals leak into development")
    dev_families = {rows[i]["family"] for i in np.concatenate([train_idx, val_idx])}
    eval_families = {rows[i]["family"] for i in eval_idx}
    if dev_families & eval_families:
        raise SystemExit("family leakage between development and evaluation")
    print(
        f"roles bound: dev_train={len(train_idx)} dev_validation={len(val_idx)} "
        f"{evaluation_role}={len(eval_idx)}",
        flush=True,
    )
    return train_idx, val_idx, eval_idx


def run_analysis(
    config: dict[str, object],
    rows: list[dict[str, str]],
    positions: np.ndarray,
    genotype: np.ndarray,
    load_records,
    train_idx: np.ndarray,
    val_idx: np.ndarray,
    eval_idx: np.ndarray,
    result_dir: Path,
) -> dict[str, object]:
    """Everything after asset binding. ``load_records`` returns expression rows
    for exactly the sample indices it is handed, so sealed outcomes stay unread.
    """
    af_all, missing_all, dosage_all = af_missingness_dosage(genotype, train_idx)
    maf_all = np.minimum(af_all, 1.0 - af_all)
    qc = np.flatnonzero(
        (maf_all >= float(config["maf_min"]))
        & (missing_all <= float(config["missingness_max"]))
        & np.isfinite(af_all)
    )
    print(f"QC variants on dev_train: {len(qc)} of {len(positions)}", flush=True)

    block_size = int(config["block_size"])

    # ---- covariates -------------------------------------------------------
    pc_panel = uniform_thin(qc, min(int(config["pc_panel_snps"]), len(qc)))
    pcs, pc_ratio = genotype_pcs(
        dosage_all[:, pc_panel], train_idx, int(config["genotype_pc_count"])
    )
    covariates = np.column_stack([demographic_covariates(rows), pcs]).astype(np.float32)

    # ---- shared encoder ---------------------------------------------------
    pretrain_snps = int(config["pretrain_blocks"]) * block_size
    if len(qc) < pretrain_snps:
        raise SystemExit("not enough QC variants for the frozen pretraining panel")
    pretrain_idx = uniform_thin(qc, pretrain_snps)
    print(f"pretraining shared encoder on {pretrain_snps} SNPs in {len(pretrain_idx)//block_size} blocks", flush=True)
    encoder, model_cfg, history, pretrain_metrics = pretrain_shared_encoder(
        genotype[:, pretrain_idx],
        af_all[pretrain_idx],
        positions[pretrain_idx],
        train_idx,
        val_idx,
        config,
    )
    checkpoint_path = result_dir / "ENCODER_FROZEN.pt"
    torch.save({"state_dict": encoder.state_dict(), "config": model_cfg.__dict__}, checkpoint_path)
    for parameter in encoder.parameters():
        parameter.requires_grad_(False)

    # ---- genes ------------------------------------------------------------
    dev_idx = np.concatenate([train_idx, val_idx])
    outcome_idx = np.concatenate([dev_idx, eval_idx])
    records = load_records(outcome_idx)
    dev_slice = slice(0, len(dev_idx))
    eval_slice = slice(len(dev_idx), len(outcome_idx))

    snps_per_gene = int(config["snps_per_gene"])
    if snps_per_gene % block_size:
        raise SystemExit("snps_per_gene must be a multiple of block_size")
    radius = int(config["cis_radius_bp"])

    eligible: dict[str, dict[str, object]] = {}
    for record in records:
        dev_values = np.asarray(record["values"])[dev_slice]
        if not np.isfinite(np.asarray(record["values"])).all() or float(np.var(dev_values)) <= 0:
            continue
        tss = int(record["tss_1based"])
        window = qc[(positions[qc] >= tss - radius) & (positions[qc] <= tss + radius)]
        if len(window) < snps_per_gene:
            continue
        chosen = uniform_thin(window, snps_per_gene)
        eligible[str(record["gene_id"])] = {
            "tss_1based": tss,
            "variant_indices": chosen,
            "cis_qc_variant_count": int(len(window)),
            "values": np.asarray(record["values"]),
        }
    ordered = deterministic_gene_order(sorted(eligible), str(config["gene_selection_salt"]))
    selected = ordered[: int(config["max_genes"])]
    if len(selected) < int(config["min_genes_required"]):
        raise SystemExit(
            f"DESIGN-INELIGIBLE: only {len(selected)} eligible genes, "
            f"{config['min_genes_required']} required"
        )
    selected.sort()
    print(f"genes: {len(eligible)} eligible, {len(selected)} analysed", flush=True)

    # ---- per-gene arms ----------------------------------------------------
    arms = ["A", "B", "C_gene", "C_full"]
    truth = np.zeros((len(selected), len(eval_idx)))
    predictions = {arm: np.zeros((len(selected), len(eval_idx))) for arm in arms}
    train_means = np.zeros(len(selected))
    score_model = np.zeros(len(selected))
    score_linear = np.zeros(len(selected))
    target_t = np.zeros(len(selected))
    gene_rows: list[dict[str, object]] = []
    alphas = [float(a) for a in config["ridge_alphas"]]
    kernel_bp = float(config["gene_pool_kernel_bp"])
    n_dev_train = len(train_idx)

    for g, gene in enumerate(selected):
        item = eligible[gene]
        snp_idx = np.asarray(item["variant_indices"], dtype=np.int64)
        gene_positions = positions[snp_idx]
        panel_genotype = genotype[:, snp_idx]
        panel_af = af_all[snp_idx]
        panel_dosage = dosage_all[:, snp_idx]

        hidden = encode_window(encoder, panel_genotype, panel_af, gene_positions, block_size, False)
        hidden_masked = encode_window(
            encoder, panel_genotype, panel_af, gene_positions, block_size, MASK_EVERYTHING
        )
        pooled = gene_pool_features(hidden, gene_positions, int(item["tss_1based"]), kernel_bp)
        pooled_masked = gene_pool_features(
            hidden_masked, gene_positions, int(item["tss_1based"]), kernel_bp
        )
        flat = hidden.reshape(len(hidden), -1)
        flat_masked = hidden_masked.reshape(len(hidden_masked), -1)

        dosage_mean = panel_dosage[train_idx].mean(axis=0, dtype=np.float64)
        dosage_perturbed = np.broadcast_to(dosage_mean, panel_dosage.shape).astype(np.float32)

        features = {
            "A": covariates,
            "B": np.column_stack([covariates, panel_dosage]),
            "C_gene": np.column_stack([covariates, pooled]),
            "C_full": np.column_stack([covariates, flat]),
        }
        perturbed = {
            "B": np.column_stack([covariates, dosage_perturbed]),
            "C_gene": np.column_stack([covariates, pooled_masked]),
            "C_full": np.column_stack([covariates, flat_masked]),
        }

        y_all = np.asarray(item["values"])
        y_dev = y_all[dev_slice]
        y_train = y_dev[:n_dev_train]
        y_val = y_dev[n_dev_train:]
        y_eval = y_all[eval_slice]
        truth[g] = y_eval
        train_means[g] = float(y_train.mean())

        info_by_arm: dict[str, dict[str, object]] = {}
        for arm in arms:
            predict, info = fit_ridge(features[arm], y_train, train_idx, val_idx, alphas, y_val)
            predictions[arm][g] = predict(features[arm][eval_idx])
            info_by_arm[arm] = info
            if arm == "B":
                target_t[g] = r2_against_train_mean(
                    y_val, predict(features[arm][val_idx]), train_means[g]
                )
            if arm in perturbed:
                delta = predictions[arm][g] - predict(perturbed[arm][eval_idx])
                if arm == "B":
                    score_linear[g] = float(np.std(delta, ddof=1))
                elif arm == "C_gene":
                    score_model[g] = float(np.std(delta, ddof=1))

        gene_rows.append(
            {
                "gene_id": gene,
                "tss_1based": int(item["tss_1based"]),
                "cis_qc_variant_count": int(item["cis_qc_variant_count"]),
                "mean_maf": float(np.mean(np.minimum(panel_af, 1.0 - panel_af))),
                "dev_expression_variance": float(np.var(y_dev)),
                "train_mean": train_means[g],
                "target_T_g_dev_validation_r2_B": float(target_t[g]),
                "score_model_C_gene": float(score_model[g]),
                "score_linear_B": float(score_linear[g]),
                "alpha": {arm: info_by_arm[arm]["alpha"] for arm in arms},
                "input_features": {arm: info_by_arm[arm]["input_features"] for arm in arms},
            }
        )
        if (g + 1) % 10 == 0 or g + 1 == len(selected):
            print(f"  gene {g + 1}/{len(selected)} done", flush=True)

    # ---- estimands --------------------------------------------------------
    primary = two_way_clustered_bootstrap(
        truth,
        predictions,
        train_means,
        ("C_gene", "B"),
        int(config["bootstrap_seed"]),
        int(config["bootstrap_replicates"]),
    )
    secondary = {}
    for name in config["secondary_contrasts"]:
        left, right = name.split("-")
        secondary[name] = two_way_clustered_bootstrap(
            truth,
            predictions,
            train_means,
            (left, right),
            int(config["bootstrap_seed"]) + 1,
            int(config["bootstrap_replicates"]),
        )
        secondary[name].pop("per_gene_r2", None)

    strata_inputs = np.column_stack(
        [
            decile_strata(np.asarray([r["cis_qc_variant_count"] for r in gene_rows], dtype=float)),
            decile_strata(np.asarray([r["mean_maf"] for r in gene_rows], dtype=float)),
            decile_strata(np.asarray([r["dev_expression_variance"] for r in gene_rows], dtype=float)),
        ]
    )
    _, strata = np.unique(strata_inputs, axis=0, return_inverse=True)
    strata = strata.astype(np.int64)

    rehearsal = {
        "target_definition": config["gene_effect_rehearsal"]["target"],
        "model": matched_null_spearman(
            score_model,
            target_t,
            strata,
            int(config["permutation_seed"]),
            int(config["permutation_replicates"]),
        ),
        "linear": matched_null_spearman(
            score_linear,
            target_t,
            strata,
            int(config["permutation_seed"]) + 1,
            int(config["permutation_replicates"]),
        ),
        "model_minus_linear": paired_gene_bootstrap_difference(
            score_model,
            score_linear,
            target_t,
            int(config["bootstrap_seed"]) + 2,
            int(config["bootstrap_replicates"]),
        ),
        "spearman_model_vs_linear_score": spearman(score_model, score_linear),
        "is_burden_test": False,
    }

    # ---- artifacts --------------------------------------------------------
    with (result_dir / "GENE_TABLE.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "gene_id",
                "tss_1based",
                "cis_qc_variant_count",
                "mean_maf",
                "dev_expression_variance",
                "target_T_g",
                "score_model_C_gene",
                "score_linear_B",
                "r2_A",
                "r2_B",
                "r2_C_gene",
                "r2_C_full",
            ]
        )
        for g, row in enumerate(gene_rows):
            writer.writerow(
                [
                    row["gene_id"],
                    row["tss_1based"],
                    row["cis_qc_variant_count"],
                    f"{row['mean_maf']:.6f}",
                    f"{row['dev_expression_variance']:.6f}",
                    f"{row['target_T_g_dev_validation_r2_B']:.6f}",
                    f"{row['score_model_C_gene']:.6f}",
                    f"{row['score_linear_B']:.6f}",
                    *[f"{primary['per_gene_r2'][arm][g]:.6f}" for arm in arms],
                ]
            )

    evaluation_role = str(config["evaluation_role"])
    sealed = set(config["sealed_roles"])
    results = {
        "analysis_id": config["analysis_id"],
        "pilot_id": config["pilot_id"],
        "claim_class": config["claim_class"],
        "status": "VALID_COMPLETE",
        "evaluation_role": evaluation_role,
        "sealed_roles_opened": [evaluation_role] if evaluation_role in sealed else [],
        "split": {
            "dev_train": len(train_idx),
            "dev_validation": len(val_idx),
            "evaluation": len(eval_idx),
            "family_disjoint": True,
        },
        "qc": {
            "qc_variants": int(len(qc)),
            "maf_min": float(config["maf_min"]),
            "missingness_max": float(config["missingness_max"]),
            "scope": "dev_train_only",
        },
        "covariates": {
            "genotype_pc_count": int(pcs.shape[1]),
            "pc_panel_snps": int(len(pc_panel)),
            "pc_explained_variance_ratio": pc_ratio,
            "sex_coding": "female=1_male=0",
            "population_columns": POPULATIONS[1:],
        },
        "encoder": {
            "shared_across_genes": True,
            "config": model_cfg.__dict__,
            "pretrain_blocks": int(len(pretrain_idx) // block_size),
            "pretrain_snps": int(len(pretrain_idx)),
            "epochs": int(config["pretrain_epochs"]),
            "validation_metrics": pretrain_metrics,
            "checkpoint_sha256": sha256_file(checkpoint_path),
            "phenotype_seen_during_pretraining": False,
        },
        "genes": {
            "eligible": int(len(eligible)),
            "analysed": int(len(selected)),
            "selection_rule": config["gene_selection_rule"],
            "selection_uses_phenotype_signal": False,
            "cis_snp_selection_rule": config["cis_snp_selection_rule"],
            "snps_per_gene": snps_per_gene,
        },
        "primary_contrast": config["primary_contrast"],
        "primary": primary,
        "primary_support": bool(primary["pass"]),
        "secondary": secondary,
        "gene_effect_rehearsal": rehearsal,
        "pretrain_history": history,
        "interpretation_limits": config["interpretation_limits"],
    }
    results["primary"].pop("per_gene_r2", None)
    dump_json(result_dir / "RESULTS.json", results)
    return results


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--open-sealed-test",
        action="store_true",
        help="Open a sealed role. Requires OPEN_AUTHORIZATION.json and burns the asset.",
    )
    args = parser.parse_args()

    if (RESULT_DIR / "FINAL_STATUS.json").exists():
        raise SystemExit("this Pilot already has a final status and cannot be rerun")
    RESULT_DIR.mkdir(parents=True, exist_ok=True)

    config = json.loads(CONFIG_PATH.read_text(encoding="utf-8-sig"))
    binding = json.loads(BINDING_PATH.read_text(encoding="utf-8-sig"))
    set_reproducibility(int(config["seed"]))

    evaluation_role = str(config["evaluation_role"])
    sealed = set(config["sealed_roles"])
    if evaluation_role in sealed:
        if not args.open_sealed_test:
            raise SystemExit(
                f"evaluation_role '{evaluation_role}' is sealed. Rerun with --open-sealed-test "
                "and a written OPEN_AUTHORIZATION.json if you really mean to burn it."
            )
        if not AUTHORIZATION_PATH.exists():
            raise SystemExit("--open-sealed-test requires OPEN_AUTHORIZATION.json next to this script")

    paths = resolve_assets(binding)
    asset_audit = audit_assets(paths, binding)
    print(f"asset audit passed for {len(asset_audit)} files", flush=True)

    manifest_path = (PILOT_DIR / str(binding["sample_manifest_path"])).resolve()
    manifest_hash = sha256_file(manifest_path)
    want_manifest = binding.get("sample_manifest_expected_sha256")
    if want_manifest and want_manifest != manifest_hash:
        raise SystemExit("SAMPLE_MANIFEST.tsv hash mismatch: the frozen split changed")
    rows = load_tsv(manifest_path)

    psam_ids = [row["#IID"] for row in load_tsv(paths["psam"])]
    if [row["sample_id"] for row in rows] != psam_ids:
        raise SystemExit("manifest/PSAM sample order mismatch")
    if expression_header(paths["expression"]) != psam_ids:
        raise SystemExit("expression/PSAM sample order mismatch")

    train_idx, val_idx, eval_idx = resolve_roles(rows, evaluation_role)

    positions, _variant_ids = read_pvar(
        paths["pvar"], str(config["chromosome"]), int(binding["expected_variant_count"])
    )
    genotype = read_pgen(paths["pgen"], len(rows), len(positions))
    print(f"genotype loaded: {genotype.shape}", flush=True)

    def load_records(outcome_idx: np.ndarray) -> list[dict[str, object]]:
        return load_expression(
            paths["expression"], str(config["chromosome"]), [psam_ids[i] for i in outcome_idx]
        )

    results = run_analysis(
        config, rows, positions, genotype, load_records, train_idx, val_idx, eval_idx, RESULT_DIR
    )
    primary = results["primary"]
    rehearsal = results["gene_effect_rehearsal"]

    dump_json(
        RESULT_DIR / "RUN_BINDING.json",
        {
            "assets": asset_audit,
            "asset_dir": str(paths["pgen"].parent),
            "sample_manifest_sha256": manifest_hash,
            "config_sha256": sha256_file(CONFIG_PATH),
            "protocol_sha256": sha256_file(PROTOCOL_PATH) if PROTOCOL_PATH.exists() else None,
            "core_sha256": sha256_file(PILOT_DIR / "tqg1_core.py"),
            "runner_sha256": sha256_file(Path(__file__).resolve()),
            "runtime": {
                "python": platform.python_version(),
                "platform": platform.platform(),
                "numpy": np.__version__,
                "sklearn": sklearn.__version__,
                "torch": torch.__version__,
                "cuda": bool(torch.cuda.is_available()),
                "torch_threads": torch.get_num_threads(),
            },
        },
    )

    dump_json(
        RESULT_DIR / "FINAL_STATUS.json",
        {
            "analysis_id": config["analysis_id"],
            "status": "SUPPORTED" if primary["pass"] else "NOT-SUPPORTED",
            "primary_contrast": config["primary_contrast"],
            "macro_delta_r2": primary["macro_delta_r2"],
            "paired_bootstrap_ci95": primary["paired_bootstrap_ci95"],
            "n_genes": primary["n_genes"],
            "n_individuals": primary["n_individuals"],
            "evaluation_role": evaluation_role,
            "sealed_roles_still_sealed": sorted(sealed - {evaluation_role}),
            "rehearsal_model_spearman": rehearsal["model"]["observed_spearman"],
            "rehearsal_linear_spearman": rehearsal["linear"]["observed_spearman"],
            "rehearsal_model_beats_linear": rehearsal["model_minus_linear"]["pass"],
            "is_burden_test": False,
        },
    )
    print("\nTQ-G1 complete. Primary:", config["primary_contrast"],
          f"delta={primary['macro_delta_r2']:.5f}",
          f"CI95={primary['paired_bootstrap_ci95']}",
          "PASS" if primary["pass"] else "FAIL", flush=True)


if __name__ == "__main__":
    main()
