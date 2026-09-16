from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from sklearn.model_selection import StratifiedGroupKFold


HERE = Path(__file__).resolve().parent
PACKAGE = HERE.parents[1]
SOURCE = PACKAGE / "pilots" / "20260915_geuvadis_chr18_stage2b_rc1"
RC1 = PACKAGE / "pilots" / "20260915_tq2_ab_empirical_power_calibration_rc1"
RESULTS = HERE / "results"
CONFIG = json.loads((HERE / "CONFIG_FROZEN.json").read_text(encoding="utf-8"))


def load_stage2b():
    spec = importlib.util.spec_from_file_location("stage2b_rc2_source", SOURCE / "run_stage2b.py")
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


def standardize(values: np.ndarray, fit_idx: np.ndarray) -> np.ndarray:
    mean = float(values[fit_idx].mean())
    sd = float(values[fit_idx].std(ddof=0))
    if sd <= 1e-12:
        raise ValueError("zero variance score")
    return (values - mean) / sd


class ResidualFold:
    def __init__(self, raw: np.ndarray, rows, train_idx: np.ndarray, test_idx: np.ndarray, stage_config):
        _, _, dosage = S2.af_and_dosage(raw, train_idx)
        covariates, _ = S2.fixed_covariates(rows, dosage, train_idx, stage_config)
        self.a_train = np.column_stack([np.ones(len(train_idx)), covariates[train_idx].astype(np.float64)])
        self.a_test = np.column_stack([np.ones(len(test_idx)), covariates[test_idx].astype(np.float64)])
        beta_x = np.linalg.lstsq(self.a_train, dosage[train_idx].astype(np.float64), rcond=None)[0]
        x_train = dosage[train_idx].astype(np.float64) - self.a_train @ beta_x
        x_test = dosage[test_idx].astype(np.float64) - self.a_test @ beta_x
        mean = x_train.mean(axis=0)
        sd = x_train.std(axis=0)
        keep = sd >= 1e-8
        z_train = (x_train[:, keep] - mean[keep]) / sd[keep]
        z_test = (x_test[:, keep] - mean[keep]) / sd[keep]
        kernel = z_train @ z_train.T
        values, vectors = np.linalg.eigh(kernel)
        self.values = np.maximum(values, 0.0)
        self.vectors = vectors
        self.test_q = (z_test @ z_train.T) @ vectors
        self.leverage_basis = vectors * vectors
        self.alphas = np.asarray(stage_config["ridge_alphas"], dtype=np.float64)
        self.train_idx = train_idx
        self.test_idx = test_idx

    def predict(self, y: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
        beta_a = np.linalg.lstsq(self.a_train, y[self.train_idx], rcond=None)[0]
        pred_a_train = self.a_train @ beta_a
        pred_a_test = self.a_test @ beta_a
        residual = y[self.train_idx] - pred_a_train
        qt = self.vectors.T @ residual
        best = None
        for alpha in self.alphas:
            shrink = self.values / (self.values + alpha)
            fitted = self.vectors @ (shrink * qt)
            leverage = self.leverage_basis @ shrink
            loo = (residual - fitted) / np.maximum(1.0 - leverage, 1e-8)
            mse = float(np.mean(loo * loo))
            if best is None or mse < best[0] - 1e-12 or (abs(mse - best[0]) <= 1e-12 and alpha > best[1]):
                best = (mse, float(alpha), qt / (self.values + alpha))
        assert best is not None
        residual_test = self.test_q @ best[2]
        return pred_a_test, pred_a_test + residual_test, best[1]


def build_folds(rows):
    populations = np.asarray([row["pop"] for row in rows])
    families = np.asarray([row["family"] for row in rows])
    splitter = StratifiedGroupKFold(n_splits=int(CONFIG["outer_folds"]), shuffle=True, random_state=int(CONFIG["seed"]))
    folds = []
    seen = np.zeros(len(rows), dtype=np.int64)
    for train_idx, test_idx in splitter.split(np.zeros(len(rows)), populations, families):
        if set(families[train_idx]) & set(families[test_idx]):
            raise ValueError("family leakage in outer folds")
        seen[test_idx] += 1
        folds.append((train_idx.astype(np.int64), test_idx.astype(np.int64)))
    if not np.all(seen == 1):
        raise ValueError("outer folds do not partition individuals")
    return folds


def family_bootstrap_gate(truths, pred_a, pred_b, families, rng, replicates):
    unique = np.unique(families)
    members = [np.flatnonzero(families == family) for family in unique]

    def delta(index):
        values = []
        for y, pa, pb in zip(truths, pred_a, pred_b):
            center = float(y.mean())
            denominator = float(np.sum((y[index] - center) ** 2))
            r2a = 1.0 - float(np.sum((y[index] - pa[index]) ** 2)) / denominator
            r2b = 1.0 - float(np.sum((y[index] - pb[index]) ** 2)) / denominator
            values.append(r2b - r2a)
        return float(np.mean(values))

    all_idx = np.arange(len(families))
    point = delta(all_idx)
    boot = np.empty(replicates, dtype=np.float64)
    for b in range(replicates):
        chosen = rng.integers(0, len(unique), size=len(unique))
        boot[b] = delta(np.concatenate([members[i] for i in chosen]))
    ci = [float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975))]
    return point, ci, bool(point > 0 and ci[0] > 0)


def main() -> None:
    if RESULTS.exists():
        raise RuntimeError("single-use rc2 already has results")
    RESULTS.mkdir(parents=True)
    rc1 = json.loads((RC1 / "results" / "FINAL_STATUS.json").read_text(encoding="utf-8"))
    if rc1["status"] != "AB_DESIGN_UNDERPOWERED":
        raise RuntimeError("rc2 remediation requires the frozen rc1 failure")
    stage_config = json.loads((SOURCE / "CONFIG_FROZEN.json").read_text(encoding="utf-8"))
    trait_panel = json.loads((SOURCE / "results" / "TRAIT_PANEL.json").read_text(encoding="utf-8"))
    rows = S2.load_tsv(SOURCE / "SAMPLE_MANIFEST.tsv")
    positions, _, _, _ = S2.read_pvar()
    genotype = S2.read_pgen(len(rows), len(positions))
    folds = build_folds(rows)
    demographics = S2.base_demographics(rows, np.arange(len(rows)), np.arange(len(rows))).astype(np.float64)
    families = np.asarray([row["family"] for row in rows])

    panels = []
    for item in trait_panel:
        raw = genotype[:, np.asarray(item["variant_indices"], dtype=np.int64)]
        # Generator dosage uses a fixed phenotype-free all-sample imputation; every fitted fold recomputes train-only AF.
        _, _, generator_dosage = S2.af_and_dosage(raw, np.arange(len(rows)))
        panels.append({
            "gene_id": item["gene_id"],
            "generator_dosage": generator_dosage.astype(np.float64),
            "folds": [ResidualFold(raw, rows, train_idx, test_idx, stage_config) for train_idx, test_idx in folds],
        })

    grid_results = []
    for grid_number, delta in enumerate(CONFIG["genetic_increment_grid"]):
        seed_results = []
        for simulation in range(int(CONFIG["simulation_seeds"])):
            rng = np.random.default_rng(int(CONFIG["seed"]) + 100000 * grid_number + simulation)
            truths, predictions_a, predictions_b, oracle_deltas, selected_alphas = [], [], [], [], []
            for panel in panels:
                causal = rng.choice(panel["generator_dosage"].shape[1], size=int(CONFIG["causal_snps_per_trait"]), replace=False)
                genetic = panel["generator_dosage"][:, causal] @ rng.normal(size=len(causal))
                design = np.column_stack([np.ones(len(rows)), demographics])
                genetic = genetic - design @ np.linalg.lstsq(design, genetic, rcond=None)[0]
                genetic = standardize(genetic, np.arange(len(rows)))
                cov_score = standardize(demographics @ rng.normal(size=demographics.shape[1]), np.arange(len(rows)))
                noise = standardize(rng.normal(size=len(rows)), np.arange(len(rows)))
                cov_fraction = float(CONFIG["covariate_variance_fraction"])
                y = np.sqrt(cov_fraction) * cov_score + np.sqrt(float(delta)) * genetic + np.sqrt(1.0 - cov_fraction - float(delta)) * noise
                pa = np.empty(len(rows), dtype=np.float64)
                pb = np.empty(len(rows), dtype=np.float64)
                alphas = []
                for fold in panel["folds"]:
                    fa, fb, alpha = fold.predict(y)
                    pa[fold.test_idx] = fa
                    pb[fold.test_idx] = fb
                    alphas.append(alpha)
                truths.append(y)
                predictions_a.append(pa)
                predictions_b.append(pb)
                selected_alphas.extend(alphas)
                oracle_a = np.sqrt(cov_fraction) * cov_score
                oracle_b = oracle_a + np.sqrt(float(delta)) * genetic
                denom = float(np.sum((y - y.mean()) ** 2))
                oracle_deltas.append(
                    (1.0 - float(np.sum((y - oracle_b) ** 2)) / denom)
                    - (1.0 - float(np.sum((y - oracle_a) ** 2)) / denom)
                )
            boot_rng = np.random.default_rng(int(CONFIG["bootstrap_seed"]) + 100000 * grid_number + simulation)
            point, ci, passed = family_bootstrap_gate(
                truths, predictions_a, predictions_b, families, boot_rng, int(CONFIG["bootstrap_replicates"])
            )
            seed_results.append({
                "simulation": simulation,
                "macro_delta_B_minus_A": point,
                "ci95": ci,
                "pass": passed,
                "mean_oracle_increment": float(np.mean(oracle_deltas)),
                "median_residual_ridge_alpha": float(np.median(selected_alphas)),
            })
        pass_rate = float(np.mean([x["pass"] for x in seed_results]))
        grid_results.append({
            "injected_increment": float(delta),
            "empirical_pass_rate": pass_rate,
            "mean_observed_B_minus_A": float(np.mean([x["macro_delta_B_minus_A"] for x in seed_results])),
            "median_observed_B_minus_A": float(np.median([x["macro_delta_B_minus_A"] for x in seed_results])),
            "mean_oracle_increment": float(np.mean([x["mean_oracle_increment"] for x in seed_results])),
            "seed_results": seed_results,
        })
        print(f"delta={delta:.2f} pass_rate={pass_rate:.3f}", flush=True)

    by_delta = {x["injected_increment"]: x for x in grid_results}
    null_fpr = by_delta[0.0]["empirical_pass_rate"]
    power = by_delta[float(CONFIG["minimum_scientific_increment"])]["empirical_pass_rate"]
    ready = null_fpr <= float(CONFIG["maximum_null_false_positive_rate"]) and power >= float(CONFIG["minimum_power"])
    final = {
        "analysis_id": CONFIG["analysis_id"],
        "status": "AB_DESIGN_POWER_READY" if ready else "AB_DESIGN_UNDERPOWERED",
        "null_false_positive_rate": null_fpr,
        "power_at_minimum_increment": power,
        "minimum_increment": float(CONFIG["minimum_scientific_increment"]),
        "required_power": float(CONFIG["minimum_power"]),
        "maximum_null_false_positive_rate": float(CONFIG["maximum_null_false_positive_rate"]),
        "real_phenotype_accessed": False,
        "encoder_or_grassmann_used": False,
        "grid": grid_results,
        "folds": [
            {
                "fold": i,
                "train_n": int(len(train_idx)),
                "test_n": int(len(test_idx)),
                "train_families": int(len(set(families[train_idx]))),
                "test_families": int(len(set(families[test_idx]))),
            }
            for i, (train_idx, test_idx) in enumerate(folds)
        ],
        "binding": {
            "config_sha256": sha256(HERE / "CONFIG_FROZEN.json"),
            "protocol_sha256": sha256(HERE / "FROZEN_PROTOCOL.zh-CN.md"),
            "rc1_final_sha256": sha256(RC1 / "results" / "FINAL_STATUS.json"),
            "source_trait_panel_sha256": sha256(SOURCE / "results" / "TRAIT_PANEL.json"),
            "source_manifest_sha256": sha256(SOURCE / "SAMPLE_MANIFEST.tsv"),
        },
    }
    dump(RESULTS / "CALIBRATION_RESULTS.json", final)
    dump(RESULTS / "FINAL_STATUS.json", {key: value for key, value in final.items() if key != "grid"})
    print(json.dumps({key: value for key, value in final.items() if key not in {"grid", "binding", "folds"}}, indent=2))


if __name__ == "__main__":
    main()

