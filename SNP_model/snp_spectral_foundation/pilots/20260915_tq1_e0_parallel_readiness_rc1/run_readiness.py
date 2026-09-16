from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import math
import sys
from dataclasses import fields
from pathlib import Path

import numpy as np
import torch


HERE = Path(__file__).resolve().parent
PACKAGE = HERE.parents[1]
SOURCE = PACKAGE / "pilots" / "20260915_geuvadis_chr18_stage2b_rc1"
RESULTS = HERE / "results"
CONFIG = json.loads((HERE / "CONFIG_FROZEN.json").read_text(encoding="utf-8"))


def load_stage2b_module():
    spec = importlib.util.spec_from_file_location("stage2b_frozen_source", SOURCE / "run_stage2b.py")
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load frozen Stage 2B implementation")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


S2 = load_stage2b_module()


def dump(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def cluster_bootstrap_delta(
    per_person_left: np.ndarray,
    per_person_right: np.ndarray,
    families: np.ndarray,
    seed: int,
    replicates: int,
) -> tuple[float, list[float]]:
    """Inputs have shape [traits, people]; resampling unit is family."""
    if per_person_left.shape != per_person_right.shape or per_person_left.shape[1] != len(families):
        raise ValueError("cluster bootstrap shape mismatch")
    unique = np.unique(families)
    members = [np.flatnonzero(families == family) for family in unique]
    point = float(np.mean(per_person_left - per_person_right))
    rng = np.random.default_rng(seed)
    values = np.empty(replicates, dtype=np.float64)
    for b in range(replicates):
        chosen = rng.integers(0, len(unique), size=len(unique))
        index = np.concatenate([members[i] for i in chosen])
        values[b] = float(np.mean(per_person_left[:, index] - per_person_right[:, index]))
    return point, [float(np.quantile(values, 0.025)), float(np.quantile(values, 0.975))]


def tq1(rows, genotype, trait_panel, stage_config, psam_ids, train_idx, val_idx):
    bridge_idx = np.asarray([i for i, row in enumerate(rows) if row["primary_role"] == "bridge_test"], dtype=np.int64)
    if len(bridge_idx) != 230:
        raise ValueError("TQ1 requires the frozen 230-person pool")
    dev_order = np.concatenate([train_idx, val_idx])
    genes = {item["gene_id"] for item in trait_panel}
    dev_records = {x["gene_id"]: x for x in S2.load_expression_subset([psam_ids[i] for i in dev_order], genes)}
    predictions = {"A": [], "B": []}
    train_means = []
    decoder_audit = []
    for item in trait_panel:
        panel_idx = np.asarray(item["variant_indices"], dtype=np.int64)
        raw_panel = genotype[:, panel_idx]
        _, _, dosage = S2.af_and_dosage(raw_panel, train_idx)
        covariates, cov_info = S2.fixed_covariates(rows, dosage, train_idx, stage_config)
        features = {"A": covariates, "B": np.concatenate([covariates, dosage], axis=1)}
        y_dev = np.asarray(dev_records[item["gene_id"]]["values"], dtype=np.float64)
        arm_record = {"gene_id": item["gene_id"], "covariates": cov_info, "decoders": {}}
        for arm in ("A", "B"):
            pred, fitted = S2.fit_ridge_predictions(
                features[arm], y_dev, train_idx, val_idx, bridge_idx, list(stage_config["ridge_alphas"])
            )
            predictions[arm].append(pred)
            arm_record["decoders"][arm] = fitted
        train_means.append(float(y_dev[: len(train_idx)].mean()))
        decoder_audit.append(arm_record)

    # This is the only phenotype opening for the old 230-person pool.
    bridge_records = {
        x["gene_id"]: x for x in S2.load_expression_subset([psam_ids[i] for i in bridge_idx], genes)
    }
    truth = [np.asarray(bridge_records[item["gene_id"]]["values"], dtype=np.float64) for item in trait_panel]
    primary = S2.summarize_gate(
        truth,
        predictions,
        train_means,
        ("B", "A"),
        int(CONFIG["bootstrap_seed"]),
        int(CONFIG["bootstrap_replicates"]),
    )
    primary["status"] = "PASS" if primary["pass"] else "TASK-INELIGIBLE"

    # Diagnostic family-cluster CI on individual squared-error improvements, normalized per trait.
    a_gain, b_gain = [], []
    for k, y in enumerate(truth):
        denom = float(np.sum((y - train_means[k]) ** 2))
        a_gain.append(-((y - predictions["A"][k]) ** 2) * len(y) / denom)
        b_gain.append(-((y - predictions["B"][k]) ** 2) * len(y) / denom)
    families = np.asarray([rows[i]["family"] for i in bridge_idx])
    family_point, family_ci = cluster_bootstrap_delta(
        np.stack(b_gain), np.stack(a_gain), families, int(CONFIG["bootstrap_seed"]) + 1, int(CONFIG["bootstrap_replicates"])
    )
    primary["family_cluster_sensitivity"] = {
        "estimand": "mean normalized individual squared-error improvement, equal trait weight",
        "point": family_point,
        "ci95": family_ci,
        "changes_primary_decision": False,
        "families": int(len(np.unique(families))),
    }
    primary["pool_status_after_run"] = "CONSUMED_BY_TQ1_PROHIBITED_FOR_CDE_BRIDGE"
    return primary, decoder_audit


def population_matched_donor_indices(rows, holdout_idx: np.ndarray) -> np.ndarray:
    donors = np.empty(len(holdout_idx), dtype=np.int64)
    for pop in sorted({rows[i]["pop"] for i in holdout_idx}):
        local = np.flatnonzero(np.asarray([rows[i]["pop"] == pop for i in holdout_idx]))
        if len(local) < 2:
            raise ValueError(f"population {pop} has fewer than two E0 subjects")
        donors[local] = np.roll(local, 1)
    if np.any(donors == np.arange(len(holdout_idx))):
        raise ValueError("donor derangement failed")
    return donors


def e0(rows, genotype, trait_panel, stage_config, train_idx):
    from snp_spectral_foundation.config import ModelConfig
    from snp_spectral_foundation.model import LocalGenotypeEncoder, encode_genotypes

    holdout_idx = np.asarray(
        [i for i, row in enumerate(rows) if row["primary_role"] in {"task_gate", "bridge_test"}], dtype=np.int64
    )
    if len(holdout_idx) != 280:
        raise ValueError("E0 requires 280 non-development, non-RC2 individuals")
    donor_local = population_matched_donor_indices(rows, holdout_idx)
    families = np.asarray([rows[i]["family"] for i in holdout_idx])
    mask_reps = int(CONFIG["e0_mask_replicates"])
    mask_probability = float(CONFIG["e0_mask_probability"])
    smoothing = float(CONFIG["e0_empirical_smoothing"])
    panel_results = []
    per_person = {"empirical_ce": [], "real_ce": [], "donor_ce": [], "empirical_acc": [], "real_acc": [], "donor_acc": []}

    for panel_number, item in enumerate(trait_panel):
        variant_idx = np.asarray(item["variant_indices"], dtype=np.int64)
        raw = genotype[:, variant_idx]
        af, _, _ = S2.af_and_dosage(raw, train_idx)
        train = raw[train_idx]
        counts = np.stack([(train == g).sum(axis=0) for g in range(3)], axis=-1).astype(np.float64)
        probs = (counts + smoothing) / (counts.sum(axis=-1, keepdims=True) + 3.0 * smoothing)

        checkpoint_path = SOURCE / "results" / "checkpoints" / f"{item['gene_id']}.pt"
        checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
        expected_variants = list(item["variant_ids"])
        if checkpoint["gene_id"] != item["gene_id"] or checkpoint["variant_ids"] != expected_variants:
            raise ValueError("checkpoint/panel identity mismatch")
        cfg_values = checkpoint["config"]
        allowed = {field.name for field in fields(ModelConfig)}
        model = LocalGenotypeEncoder(ModelConfig(**{k: v for k, v in cfg_values.items() if k in allowed}))
        model.load_state_dict(checkpoint["state_dict"])
        model.eval()

        n = len(holdout_idx)
        sums = {name: np.zeros(n, dtype=np.float64) for name in ("empirical_ce", "real_ce", "donor_ce", "empirical_correct", "real_correct", "donor_correct", "count")}
        local_raw = raw[holdout_idx]
        donor_raw = local_raw[donor_local]
        blocks = local_raw.shape[1] // int(stage_config["block_size"])
        target_tensor = torch.from_numpy(local_raw.reshape(n, blocks, int(stage_config["block_size"])).astype(np.int64))
        donor_tensor = torch.from_numpy(donor_raw.reshape(n, blocks, int(stage_config["block_size"])).astype(np.int64))
        af_tensor = torch.from_numpy(af.reshape(blocks, int(stage_config["block_size"])).astype(np.float32))
        prob_tensor = torch.from_numpy(probs.reshape(blocks, int(stage_config["block_size"]), 3).astype(np.float32))

        for rep in range(mask_reps):
            generator = torch.Generator(device="cpu").manual_seed(int(CONFIG["seed"]) + 1000 * panel_number + rep)
            mask = (torch.rand(target_tensor.shape, generator=generator) < mask_probability) & (target_tensor >= 0)
            with torch.no_grad():
                _, real_logits = model(encode_genotypes(target_tensor, mask), af_tensor)
                _, donor_logits = model(encode_genotypes(donor_tensor, mask), af_tensor)
                real_logp = torch.log_softmax(real_logits, dim=-1)
                donor_logp = torch.log_softmax(donor_logits, dim=-1)
            target = target_tensor.clamp_min(0).long()
            empirical_logp = torch.log(prob_tensor.clamp_min(1e-12)).unsqueeze(0).expand(n, -1, -1, -1)
            for i in range(n):
                valid = mask[i]
                count = int(valid.sum())
                if count == 0:
                    continue
                idx = target[i][valid]
                sums["real_ce"][i] += float(-real_logp[i][valid].gather(1, idx[:, None]).sum())
                sums["donor_ce"][i] += float(-donor_logp[i][valid].gather(1, idx[:, None]).sum())
                sums["empirical_ce"][i] += float(-empirical_logp[i][valid].gather(1, idx[:, None]).sum())
                sums["real_correct"][i] += float((real_logits[i][valid].argmax(-1) == idx).sum())
                sums["donor_correct"][i] += float((donor_logits[i][valid].argmax(-1) == idx).sum())
                sums["empirical_correct"][i] += float((prob_tensor[valid].argmax(-1) == idx).sum())
                sums["count"][i] += count
        if np.any(sums["count"] == 0):
            raise ValueError("E0 produced an individual with no masked targets")
        metrics = {
            "empirical_ce": sums["empirical_ce"] / sums["count"],
            "real_ce": sums["real_ce"] / sums["count"],
            "donor_ce": sums["donor_ce"] / sums["count"],
            "empirical_acc": sums["empirical_correct"] / sums["count"],
            "real_acc": sums["real_correct"] / sums["count"],
            "donor_acc": sums["donor_correct"] / sums["count"],
        }
        for key, value in metrics.items():
            per_person[key].append(value)
        panel_results.append({
            "gene_id": item["gene_id"],
            "checkpoint_sha256": sha256(checkpoint_path),
            "n_individuals": n,
            "mean_masked_targets_per_individual": float(np.mean(sums["count"])),
            **{key: float(np.mean(value)) for key, value in metrics.items()},
            "ce_empirical_minus_real": float(np.mean(metrics["empirical_ce"] - metrics["real_ce"])),
            "ce_donor_minus_real": float(np.mean(metrics["donor_ce"] - metrics["real_ce"])),
        })

    stacked = {key: np.stack(value) for key, value in per_person.items()}
    point_marginal, ci_marginal = cluster_bootstrap_delta(
        stacked["empirical_ce"], stacked["real_ce"], families, int(CONFIG["bootstrap_seed"]) + 10, int(CONFIG["bootstrap_replicates"])
    )
    point_context, ci_context = cluster_bootstrap_delta(
        stacked["donor_ce"], stacked["real_ce"], families, int(CONFIG["bootstrap_seed"]) + 11, int(CONFIG["bootstrap_replicates"])
    )
    passed = point_marginal > 0 and ci_marginal[0] > 0 and point_context > 0 and ci_context[0] > 0
    return {
        "status": "PASS" if passed else "FAIL",
        "pass": passed,
        "n_individuals": int(len(holdout_idx)),
        "n_families": int(len(np.unique(families))),
        "panels": panel_results,
        "macro": {
            "empirical_ce": float(np.mean(stacked["empirical_ce"])),
            "real_context_ce": float(np.mean(stacked["real_ce"])),
            "donor_context_ce": float(np.mean(stacked["donor_ce"])),
            "empirical_accuracy": float(np.mean(stacked["empirical_acc"])),
            "real_context_accuracy": float(np.mean(stacked["real_acc"])),
            "donor_context_accuracy": float(np.mean(stacked["donor_acc"])),
        },
        "co_primary": {
            "ce_empirical_minus_real": {"point": point_marginal, "family_cluster_ci95": ci_marginal},
            "ce_donor_minus_real": {"point": point_context, "family_cluster_ci95": ci_context},
        },
        "claim_boundary": "current eight local encoders and frozen locus panels only; not genome-wide or phenotype evidence",
    }


def main() -> None:
    if RESULTS.exists():
        raise RuntimeError("results directory already exists; frozen run is single-use")
    RESULTS.mkdir(parents=True)
    stage_final = json.loads((SOURCE / "results" / "FINAL_STATUS.json").read_text(encoding="utf-8"))
    if stage_final["status"] != "TASK-INELIGIBLE" or stage_final["bridge_gate_opened"]:
        raise RuntimeError("source Stage 2B boundary is not the expected unopened bridge state")
    stage_config = json.loads((SOURCE / "CONFIG_FROZEN.json").read_text(encoding="utf-8"))
    trait_panel = json.loads((SOURCE / "results" / "TRAIT_PANEL.json").read_text(encoding="utf-8"))
    rows = S2.load_tsv(SOURCE / "SAMPLE_MANIFEST.tsv")
    psam_ids = [row["#IID"] for row in S2.load_tsv(S2.PSAM)]
    if [row["sample_id"] for row in rows] != psam_ids or S2.expression_header() != psam_ids:
        raise ValueError("sample order mismatch")
    train_idx = np.asarray([i for i, row in enumerate(rows) if row["development_role"] == "dev_train"], dtype=np.int64)
    val_idx = np.asarray([i for i, row in enumerate(rows) if row["development_role"] == "dev_validation"], dtype=np.int64)
    positions, _, _, _ = S2.read_pvar()
    genotype = S2.read_pgen(len(rows), len(positions))

    binding = {
        "config_sha256": sha256(HERE / "CONFIG_FROZEN.json"),
        "protocol_sha256": sha256(HERE / "FROZEN_PROTOCOL.zh-CN.md"),
        "source_final_sha256": sha256(SOURCE / "results" / "FINAL_STATUS.json"),
        "source_trait_panel_sha256": sha256(SOURCE / "results" / "TRAIT_PANEL.json"),
        "source_sample_manifest_sha256": sha256(SOURCE / "SAMPLE_MANIFEST.tsv"),
        "source_stage_config_sha256": sha256(SOURCE / "CONFIG_FROZEN.json"),
    }
    dump(RESULTS / "RUN_BINDING.json", binding)

    e0_result = e0(rows, genotype, trait_panel, stage_config, train_idx)
    dump(RESULTS / "E0_RESULTS.json", e0_result)
    tq1_result, decoder_audit = tq1(rows, genotype, trait_panel, stage_config, psam_ids, train_idx, val_idx)
    dump(RESULTS / "TQ1_RESULTS.json", tq1_result)
    dump(RESULTS / "TQ1_DECODER_AUDIT.json", decoder_audit)

    if tq1_result["pass"] and e0_result["pass"]:
        status = "READY_FOR_NEW_BRIDGE"
    elif not tq1_result["pass"] and not e0_result["pass"]:
        status = "NOT_READY_BOTH"
    elif not tq1_result["pass"]:
        status = "NOT_READY_TASK"
    else:
        status = "NOT_READY_ENCODER"
    final = {
        "analysis_id": CONFIG["analysis_id"],
        "status": status,
        "tq1_status": tq1_result["status"],
        "e0_status": e0_result["status"],
        "old_geuvadis_bridge_pool": "CONSUMED_BY_TQ1_PROHIBITED_FOR_CDE_BRIDGE",
        "new_bridge_authorized": status == "READY_FOR_NEW_BRIDGE",
        "new_bridge_requirement": "fresh phenotype individuals/cohort and a separately frozen contract",
        "grassmann_evaluated": False,
    }
    dump(RESULTS / "FINAL_STATUS.json", final)
    print(json.dumps(final, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()

