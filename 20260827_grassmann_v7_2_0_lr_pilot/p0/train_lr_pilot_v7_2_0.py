from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import torch
from torch.nn import functional as F

from masking_v7_1_11 import build_mask_bank, ld_macroblocks, stable_int
from profile_models_v7_1 import Architecture, MaskedGenotypeModel, parameter_count
from train_c0_real_v7_1_11 import evaluate


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def gradient_norm(model: torch.nn.Module) -> float:
    total = torch.zeros((), device=next(model.parameters()).device, dtype=torch.float64)
    for parameter in model.parameters():
        if parameter.grad is not None:
            total += parameter.grad.detach().double().square().sum()
    return float(total.sqrt().cpu())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", required=True)
    parser.add_argument("--mask", required=True)
    parser.add_argument("--mask-seed", type=int, required=True)
    parser.add_argument("--init-seed", type=int, required=True)
    parser.add_argument("--peak-lr", type=float, required=True)
    parser.add_argument("--steps", type=int, default=4000)
    parser.add_argument("--warmup-steps", type=int, default=500)
    parser.add_argument("--eval-interval", type=int, default=250)
    parser.add_argument("--mask-bank-size", type=int, default=256)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--firewall-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise SystemExit(f"refusing to overwrite: {args.output_dir}")
    args.output_dir.mkdir(parents=True, exist_ok=False)

    try:
        if args.peak_lr not in {0.0001, 0.0002, 0.0004, 0.0008}:
            raise ValueError("peak LR outside frozen grid")
        if args.steps != 4000 or args.warmup_steps != 500 or args.eval_interval != 250:
            raise ValueError("training horizon/schedule differs from frozen pilot")
        train = np.load(args.data_dir / "TRAIN_GT_SAMPLE_MAJOR.uint8.npy", mmap_mode="r")
        validation = np.load(args.data_dir / "VALIDATION_GT_SAMPLE_MAJOR.uint8.npy", mmap_mode="r")
        train_pc = np.load(args.data_dir / "TRAIN_PC16.float32.npy", mmap_mode="r")
        validation_pc = np.load(args.data_dir / "VALIDATION_PC16.float32.npy", mmap_mode="r")
        subset = np.load(args.data_dir / "SUBSET_100P_2247.indices.int64.npy", mmap_mode="r")
        if train.shape != (2247, 154850) or validation.shape != (249, 154850):
            raise ValueError("genotype shape mismatch")
        if train_pc.shape != (2247, 16) or validation_pc.shape != (249, 16):
            raise ValueError("PC shape mismatch")

        tuning_indices = [int(value) for value in (args.firewall_dir / "LR_TUNING.indices.txt").read_text().splitlines()]
        decision_indices = [int(value) for value in (args.firewall_dir / "DECISION_HOLDOUT.indices.txt").read_text().splitlines()]
        if not tuning_indices or set(tuning_indices) & set(decision_indices):
            raise ValueError("validation firewall invalid")
        if sorted(tuning_indices + decision_indices) != list(range(249)):
            raise ValueError("validation firewall does not partition all validation rows")

        blocks = ld_macroblocks(args.data_dir / "TRAIN_VARIANT_KEYS.txt", args.data_dir / "PCA_LD.prune.in")
        bank = build_mask_bank(args.mask, 154850, args.mask_seed, args.mask_bank_size, blocks)
        rates = bank.mean(axis=1)
        tolerance = 0.01 if args.mask == "ld_block_0p90" else 0.002
        if float(np.max(np.abs(rates - 0.90))) > tolerance:
            raise ValueError("mask-rate tolerance failed")

        torch.manual_seed(args.init_seed)
        torch.cuda.manual_seed_all(args.init_seed)
        if not torch.cuda.is_available():
            raise ValueError("CUDA unavailable")
        device = torch.device("cuda:0")
        model = MaskedGenotypeModel(args.model, Architecture()).to(device)
        optimizer = torch.optim.AdamW(model.parameters(), lr=args.peak_lr)
        config = {
            "schema_version": "1.0", "protocol_version": "v7.2.0",
            "role": "SHARED_LR_OPTIMIZATION_PILOT_ONLY",
            "model": args.model, "mask": args.mask, "mask_seed": args.mask_seed,
            "init_seed": args.init_seed, "peak_lr": args.peak_lr,
            "steps": args.steps, "warmup_steps": args.warmup_steps,
            "warmup_initial_fraction": 0.10, "post_warmup_schedule": "CONSTANT_PEAK",
            "eval_interval": args.eval_interval, "mask_bank_size": args.mask_bank_size,
            "tuning_validation_count": len(tuning_indices),
            "decision_holdout_read": False, "parameter_count": parameter_count(model),
            "cuda_visible_devices": os.environ.get("CUDA_VISIBLE_DEVICES"),
            "data_manifest_sha256": sha256(args.data_dir / "PREPROCESS_MANIFEST.v7.1.10.sha256"),
            "firewall_manifest_sha256": sha256(args.firewall_dir / "VALIDATION_FIREWALL.v7.2.0.sha256"),
            "hgdp_used": False,
        }
        (args.output_dir / "RUN_CONFIG.json").write_text(json.dumps(config, indent=2, sort_keys=True) + "\n")

        curve = args.output_dir / "PILOT_CURVE.jsonl"
        interval_loss_sum = 0.0
        interval_tokens = 0
        interval_started = time.perf_counter()
        model.train()
        torch.cuda.reset_peak_memory_stats(device)
        for step in range(1, args.steps + 1):
            fraction = 0.10 + 0.90 * min(step, args.warmup_steps) / args.warmup_steps
            lr = args.peak_lr * fraction if step <= args.warmup_steps else args.peak_lr
            for group in optimizer.param_groups:
                group["lr"] = lr
            row = int(subset[stable_int("train_sample", args.mask_seed, step) % len(subset)])
            target = torch.from_numpy(np.array(train[row], copy=True, dtype=np.int64)).unsqueeze(0).to(device)
            pc = torch.from_numpy(np.array(train_pc[row], copy=True, dtype=np.float32)).unsqueeze(0).to(device)
            mask = torch.from_numpy(np.array(bank[(step - 1) % len(bank)], copy=True)).unsqueeze(0).to(device)
            tokens = target.clone(); tokens[mask] = 3
            optimizer.zero_grad(set_to_none=True)
            with torch.autocast(device_type="cuda", dtype=torch.bfloat16):
                logits = model(tokens, pc)
                loss = F.cross_entropy(logits[mask], target[mask])
            loss_value = float(loss.detach().cpu())
            if not math.isfinite(loss_value):
                raise RuntimeError(f"nonfinite training loss at step {step}")
            if step > args.warmup_steps and loss_value > 5.0:
                raise RuntimeError(f"training loss exceeded hard ceiling at step {step}: {loss_value}")
            masked = int(mask.sum().item())
            interval_loss_sum += loss_value * masked
            interval_tokens += masked
            loss.backward()
            grad = gradient_norm(model) if step % args.eval_interval == 0 else None
            if grad is not None and not math.isfinite(grad):
                raise RuntimeError(f"nonfinite gradient norm at step {step}")
            optimizer.step()

            if step % args.eval_interval == 0:
                elapsed = time.perf_counter() - interval_started
                validation_nll, validation_tokens = evaluate(
                    model, validation, validation_pc, tuning_indices, bank, device
                )
                record = {
                    "step": step,
                    "train_masked_nll_interval": interval_loss_sum / interval_tokens,
                    "train_masked_tokens_interval": interval_tokens,
                    "validation_masked_nll": validation_nll,
                    "validation_masked_tokens": validation_tokens,
                    "validation_samples": len(tuning_indices),
                    "learning_rate": lr,
                    "gradient_norm": grad,
                    "train_masked_tokens_per_second_interval": interval_tokens / elapsed,
                    "peak_allocated_bytes": torch.cuda.max_memory_allocated(device),
                    "peak_reserved_bytes": torch.cuda.max_memory_reserved(device),
                }
                with curve.open("a", encoding="utf-8") as handle:
                    handle.write(json.dumps(record, sort_keys=True) + "\n")
                print(json.dumps(record, sort_keys=True), flush=True)
                interval_loss_sum = 0.0
                interval_tokens = 0
                interval_started = time.perf_counter()
                model.train()

        result = config | {
            "status": "PASS", "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "curve_sha256": sha256(curve), "final_step": args.steps,
            "final_peak_allocated_bytes": torch.cuda.max_memory_allocated(device),
            "final_peak_reserved_bytes": torch.cuda.max_memory_reserved(device),
        }
        (args.output_dir / "RESULT.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    except Exception as exc:
        failure = {
            "status": "FAIL", "error_type": type(exc).__name__, "error": str(exc),
            "traceback": traceback.format_exc(), "created_at_utc": datetime.now(timezone.utc).isoformat(),
        }
        (args.output_dir / "FAILURE.json").write_text(json.dumps(failure, indent=2, sort_keys=True) + "\n")
        raise


if __name__ == "__main__":
    main()
