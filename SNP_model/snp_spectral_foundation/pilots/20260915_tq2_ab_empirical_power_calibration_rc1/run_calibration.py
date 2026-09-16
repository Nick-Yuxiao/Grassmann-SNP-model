from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
PACKAGE = HERE.parents[1]
SOURCE = PACKAGE / "pilots" / "20260915_geuvadis_chr18_stage2b_rc1"
RESULTS = HERE / "results"
CONFIG = json.loads((HERE / "CONFIG_FROZEN.json").read_text(encoding="utf-8"))


def load_stage2b():
    spec = importlib.util.spec_from_file_location("stage2b_calibration_source", SOURCE / "run_stage2b.py")
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot import Stage 2B source")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


S2 = load_stage2b()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


class FrozenRidgePath:
    def __init__(self, x: np.ndarray, train_idx: np.ndarray, val_idx: np.ndarray, test_idx: np.ndarray, alphas: list[float]):
        x_train = x[train_idx].astype(np.float64)
        mean = x_train.mean(axis=0)
        sd = x_train.std(axis=0)
        keep = sd >= 1e-8
        z_train = (x_train[:, keep] - mean[keep]) / sd[keep]
        z_val = (x[val_idx][:, keep].astype(np.float64) - mean[keep]) / sd[keep]
        z_test = (x[test_idx][:, keep].astype(np.float64) - mean[keep]) / sd[keep]
        values, vectors = np.linalg.eigh(z_train @ z_train.T)
        self.values = np.maximum(values, 0.0)
        self.vectors = vectors
        self.val_q = (z_val @ z_train.T) @ vectors
        self.test_q = (z_test @ z_train.T) @ vectors
        self.alphas = np.asarray(alphas, dtype=np.float64)

    def predict(self, y_train: np.ndarray, y_val: np.ndarray) -> tuple[np.ndarray, float]:
        y_mean = float(y_train.mean())
        y_sd = float(y_train.std(ddof=0))
        if y_sd <= 0:
            raise ValueError("constant simulated training phenotype")
        target = (y_train - y_mean) / y_sd
        qt = self.vectors.T @ target
        best = None
        for alpha in self.alphas:
            weights = qt / (self.values + alpha)
            val_prediction = y_mean + y_sd * (self.val_q @ weights)
            mse = float(np.mean((y_val - val_prediction) ** 2))
            if best is None or mse < best[0] - 1e-12 or (abs(mse - best[0]) <= 1e-12 and alpha > best[1]):
                best = (mse, float(alpha), weights.copy())
        assert best is not None
        return y_mean + y_sd * (self.test_q @ best[2]), best[1]


def standardize_train(values: np.ndarray, train_idx: np.ndarray) -> np.ndarray:
    mean = float(values[train_idx].mean())
    sd = float(values[train_idx].std(ddof=0))
    if sd <= 1e-12:
        raise ValueError("zero variance score")
    return (values - mean) / sd


def residualize_train(score: np.ndarray, covariates: np.ndarray, train_idx: np.ndarray) -> np.ndarray:
    design_train = np.column_stack([np.ones(len(train_idx)), covariates[train_idx].astype(np.float64)])
    beta = np.linalg.lstsq(design_train, score[train_idx], rcond=None)[0]
    design_all = np.column_stack([np.ones(len(score)), covariates.astype(np.float64)])
    return score - design_all @ beta


def gate(truths, pred_a, pred_b, train_means, bootstrap_indices):
    trait_a, trait_b = [], []
    for y, pa, pb, train_mean in zip(truths, pred_a, pred_b, train_means):
        denom = float(np.sum((y - train_mean) ** 2))
        trait_a.append(1.0 - float(np.sum((y - pa) ** 2)) / denom)
        trait_b.append(1.0 - float(np.sum((y - pb) ** 2)) / denom)
    point = float(np.mean(trait_b) - np.mean(trait_a))
    boot = []
    for sample in bootstrap_indices:
        values = []
        for y, pa, pb, train_mean in zip(truths, pred_a, pred_b, train_means):
            denom = float(np.sum((y[sample] - train_mean) ** 2))
            r2a = 1.0 - float(np.sum((y[sample] - pa[sample]) ** 2)) / denom
            r2b = 1.0 - float(np.sum((y[sample] - pb[sample]) ** 2)) / denom
            values.append(r2b - r2a)
        boot.append(float(np.mean(values)))
    ci = [float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975))]
    return point, ci, bool(point > 0 and ci[0] > 0)


def main() -> None:
    if RESULTS.exists():
        raise RuntimeError("single-use frozen calibration already has results")
    RESULTS.mkdir(parents=True)
    stage_config = json.loads((SOURCE / "CONFIG_FROZEN.json").read_text(encoding="utf-8"))
    trait_panel = json.loads((SOURCE / "results" / "TRAIT_PANEL.json").read_text(encoding="utf-8"))
    rows = S2.load_tsv(SOURCE / "SAMPLE_MANIFEST.tsv")
    train_idx = np.asarray([i for i, row in enumerate(rows) if row["development_role"] == "dev_train"], dtype=np.int64)
    val_idx = np.asarray([i for i, row in enumerate(rows) if row["development_role"] == "dev_validation"], dtype=np.int64)
    test_idx = np.asarray([i for i, row in enumerate(rows) if row["primary_role"] == "bridge_test"], dtype=np.int64)
    if (len(train_idx), len(val_idx), len(test_idx)) != (120, 30, 230):
        raise ValueError("frozen split mismatch")
    positions, _, _, _ = S2.read_pvar()
    genotype = S2.read_pgen(len(rows), len(positions))

    panels = []
    for item in trait_panel:
        raw = genotype[:, np.asarray(item["variant_indices"], dtype=np.int64)]
        _, _, dosage = S2.af_and_dosage(raw, train_idx)
        covariates, _ = S2.fixed_covariates(rows, dosage, train_idx, stage_config)
        features = {"A": covariates, "B": np.concatenate([covariates, dosage], axis=1)}
        panels.append({
            "gene_id": item["gene_id"],
            "dosage": dosage.astype(np.float64),
            "covariates": covariates.astype(np.float64),
            "ridge_A": FrozenRidgePath(features["A"], train_idx, val_idx, test_idx, list(stage_config["ridge_alphas"])),
            "ridge_B": FrozenRidgePath(features["B"], train_idx, val_idx, test_idx, list(stage_config["ridge_alphas"])),
        })

    n_seeds = int(CONFIG["simulation_seeds"])
    n_boot = int(CONFIG["bootstrap_replicates"])
    grid_results = []
    for grid_number, delta in enumerate(CONFIG["genetic_increment_grid"]):
        seed_results = []
        for simulation in range(n_seeds):
            rng = np.random.default_rng(int(CONFIG["seed"]) + 100000 * grid_number + simulation)
            truths, pred_a, pred_b, means, oracle_deltas = [], [], [], [], []
            alpha_a, alpha_b = [], []
            for trait_number, panel in enumerate(panels):
                causal = rng.choice(panel["dosage"].shape[1], size=int(CONFIG["causal_snps_per_trait"]), replace=False)
                weights = rng.normal(size=len(causal))
                genetic = panel["dosage"][:, causal] @ weights
                genetic = residualize_train(genetic, panel["covariates"], train_idx)
                genetic = standardize_train(genetic, train_idx)
                cov_beta = rng.normal(size=panel["covariates"].shape[1])
                cov_score = standardize_train(panel["covariates"] @ cov_beta, train_idx)
                noise = standardize_train(rng.normal(size=len(rows)), train_idx)
                cov_fraction = float(CONFIG["covariate_variance_fraction"])
                y = np.sqrt(cov_fraction) * cov_score + np.sqrt(float(delta)) * genetic + np.sqrt(1.0 - cov_fraction - float(delta)) * noise
                pa, aa = panel["ridge_A"].predict(y[train_idx], y[val_idx])
                pb, ab = panel["ridge_B"].predict(y[train_idx], y[val_idx])
                y_test = y[test_idx]
                truths.append(y_test)
                pred_a.append(pa)
                pred_b.append(pb)
                means.append(float(y[train_idx].mean()))
                alpha_a.append(aa)
                alpha_b.append(ab)
                oracle_a = np.sqrt(cov_fraction) * cov_score[test_idx]
                oracle_b = oracle_a + np.sqrt(float(delta)) * genetic[test_idx]
                denominator = float(np.sum((y_test - y[train_idx].mean()) ** 2))
                oracle_deltas.append(
                    (1.0 - float(np.sum((y_test - oracle_b) ** 2)) / denominator)
                    - (1.0 - float(np.sum((y_test - oracle_a) ** 2)) / denominator)
                )
            boot_rng = np.random.default_rng(int(CONFIG["bootstrap_seed"]) + 100000 * grid_number + simulation)
            boot_indices = boot_rng.integers(0, len(test_idx), size=(n_boot, len(test_idx)))
            point, ci, passed = gate(truths, pred_a, pred_b, means, boot_indices)
            seed_results.append({
                "simulation": simulation,
                "macro_delta_B_minus_A": point,
                "ci95": ci,
                "pass": passed,
                "mean_oracle_test_increment": float(np.mean(oracle_deltas)),
                "median_alpha_A": float(np.median(alpha_a)),
                "median_alpha_B": float(np.median(alpha_b)),
            })
        pass_rate = float(np.mean([x["pass"] for x in seed_results]))
        grid_results.append({
            "injected_increment": float(delta),
            "empirical_pass_rate": pass_rate,
            "mean_observed_B_minus_A": float(np.mean([x["macro_delta_B_minus_A"] for x in seed_results])),
            "median_observed_B_minus_A": float(np.median([x["macro_delta_B_minus_A"] for x in seed_results])),
            "mean_oracle_test_increment": float(np.mean([x["mean_oracle_test_increment"] for x in seed_results])),
            "seed_results": seed_results,
        })
        print(f"delta={delta:.2f} pass_rate={pass_rate:.3f}", flush=True)

    by_delta = {x["injected_increment"]: x for x in grid_results}
    null_fpr = by_delta[0.0]["empirical_pass_rate"]
    power_at_minimum = by_delta[float(CONFIG["minimum_scientific_increment"])]["empirical_pass_rate"]
    ready = null_fpr <= float(CONFIG["maximum_null_false_positive_rate"]) and power_at_minimum >= float(CONFIG["minimum_power"])
    final = {
        "analysis_id": CONFIG["analysis_id"],
        "status": "AB_DESIGN_POWER_READY" if ready else "AB_DESIGN_UNDERPOWERED",
        "null_false_positive_rate": null_fpr,
        "power_at_minimum_increment": power_at_minimum,
        "minimum_increment": float(CONFIG["minimum_scientific_increment"]),
        "required_power": float(CONFIG["minimum_power"]),
        "maximum_null_false_positive_rate": float(CONFIG["maximum_null_false_positive_rate"]),
        "real_phenotype_accessed": False,
        "encoder_or_grassmann_used": False,
        "grid": grid_results,
        "binding": {
            "config_sha256": sha256(HERE / "CONFIG_FROZEN.json"),
            "protocol_sha256": sha256(HERE / "FROZEN_PROTOCOL.zh-CN.md"),
            "source_trait_panel_sha256": sha256(SOURCE / "results" / "TRAIT_PANEL.json"),
            "source_manifest_sha256": sha256(SOURCE / "SAMPLE_MANIFEST.tsv"),
            "source_stage_config_sha256": sha256(SOURCE / "CONFIG_FROZEN.json"),
        },
    }
    dump(RESULTS / "CALIBRATION_RESULTS.json", final)
    dump(RESULTS / "FINAL_STATUS.json", {key: value for key, value in final.items() if key not in {"grid"}})
    print(json.dumps({key: value for key, value in final.items() if key not in {"grid", "binding"}}, indent=2))


if __name__ == "__main__":
    main()

