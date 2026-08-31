from __future__ import annotations

import argparse
import gc
import json
import math
import os
import time
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import torch
from torch import nn
from torch.nn import functional as F


@dataclass(frozen=True)
class Architecture:
    model_dim: int = 448
    ff_dim: int = 1280
    layers: int = 4
    heads: int = 8
    chunk_size: int = 512
    pc_dim: int = 16
    grassmann_rank: int = 32
    genotype_states: int = 3


class ChunkedAttentionBlock(nn.Module):
    def __init__(self, cfg: Architecture, ff_dim: int | None = None) -> None:
        super().__init__()
        self.chunk_size = cfg.chunk_size
        effective_ff_dim = ff_dim or cfg.ff_dim
        self.norm1 = nn.LayerNorm(cfg.model_dim)
        self.attn = nn.MultiheadAttention(
            cfg.model_dim, cfg.heads, dropout=0.0, batch_first=True
        )
        self.norm2 = nn.LayerNorm(cfg.model_dim)
        self.ff = nn.Sequential(
            nn.Linear(cfg.model_dim, effective_ff_dim),
            nn.GELU(),
            nn.Linear(effective_ff_dim, cfg.model_dim),
        )

    def forward(self, hidden: torch.Tensor) -> torch.Tensor:
        batch, length, width = hidden.shape
        padding = (-length) % self.chunk_size
        if padding:
            hidden = F.pad(hidden, (0, 0, 0, padding))
        chunks = hidden.reshape(batch * (hidden.shape[1] // self.chunk_size), self.chunk_size, width)
        normalized = self.norm1(chunks)
        mixed, _ = self.attn(normalized, normalized, normalized, need_weights=False)
        chunks = chunks + mixed
        chunks = chunks + self.ff(self.norm2(chunks))
        return chunks.reshape(batch, -1, width)[:, :length]


class GrassmannBlockMixer(nn.Module):
    """A block-summary exterior-product channel; token attention remains local."""

    def __init__(self, cfg: Architecture) -> None:
        super().__init__()
        self.chunk_size = cfg.chunk_size
        self.reduce = nn.Linear(cfg.model_dim, cfg.grassmann_rank)
        i, j = torch.triu_indices(cfg.grassmann_rank, cfg.grassmann_rank, offset=1)
        self.register_buffer("wedge_i", i)
        self.register_buffer("wedge_j", j)
        self.project = nn.Linear(len(i), cfg.model_dim)
        self.gate = nn.Linear(2 * cfg.model_dim, cfg.model_dim)

    def forward(self, hidden: torch.Tensor) -> torch.Tensor:
        batch, length, width = hidden.shape
        padding = (-length) % self.chunk_size
        padded = F.pad(hidden, (0, 0, 0, padding)) if padding else hidden
        blocks = padded.reshape(batch, -1, self.chunk_size, width).mean(dim=2)
        reduced = self.reduce(blocks)
        context = reduced.mean(dim=1, keepdim=True)
        wedge = (
            reduced[..., self.wedge_i] * context[..., self.wedge_j]
            - reduced[..., self.wedge_j] * context[..., self.wedge_i]
        )
        wedge = F.normalize(wedge, p=2, dim=-1, eps=1e-8)
        geometry = self.project(wedge).repeat_interleave(self.chunk_size, dim=1)[:, :length]
        alpha = torch.sigmoid(self.gate(torch.cat([hidden, geometry], dim=-1)))
        return alpha * hidden + (1.0 - alpha) * geometry


class MaskedGenotypeModel(nn.Module):
    def __init__(self, kind: str, cfg: Architecture) -> None:
        super().__init__()
        if kind not in {"local_attn", "local_attn_gpc", "grassmann_full"}:
            raise ValueError(kind)
        self.kind = kind
        self.cfg = cfg
        self.token = nn.Embedding(cfg.genotype_states + 1, cfg.model_dim)
        self.local_position = nn.Embedding(cfg.chunk_size, cfg.model_dim)
        self.pc_projection = nn.Linear(cfg.pc_dim, cfg.model_dim) if kind != "local_attn" else None
        self.grassmann = GrassmannBlockMixer(cfg) if kind == "grassmann_full" else None
        # The Grassmann channel has its own parameters, so its FF width is reduced
        # to keep all three candidates near the frozen 8M matched-parameter budget.
        self.effective_ff_dim = 1104 if kind == "grassmann_full" else cfg.ff_dim
        self.blocks = nn.ModuleList([
            ChunkedAttentionBlock(cfg, ff_dim=self.effective_ff_dim) for _ in range(cfg.layers)
        ])
        self.final_norm = nn.LayerNorm(cfg.model_dim)
        self.head = nn.Linear(cfg.model_dim, cfg.genotype_states)

    def forward(self, tokens: torch.Tensor, pcs: torch.Tensor) -> torch.Tensor:
        positions = torch.arange(tokens.shape[1], device=tokens.device) % self.cfg.chunk_size
        hidden = self.token(tokens) + self.local_position(positions)
        if self.pc_projection is not None:
            hidden = hidden + self.pc_projection(pcs).unsqueeze(1)
        if self.grassmann is not None:
            hidden = self.grassmann(hidden)
        for block in self.blocks:
            hidden = block(hidden)
        return self.head(self.final_norm(hidden))


def parameter_count(model: nn.Module) -> int:
    return sum(parameter.numel() for parameter in model.parameters() if parameter.requires_grad)


def make_batch(length: int, batch_size: int, pc_dim: int, seed: int) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
    generator = torch.Generator(device="cpu").manual_seed(seed)
    target = torch.randint(0, 3, (batch_size, length), generator=generator, dtype=torch.long)
    mask = torch.rand((batch_size, length), generator=generator) < 0.15
    if not bool(mask.any()):
        mask[0, 0] = True
    tokens = target.clone()
    tokens[mask] = 3
    pcs = torch.randn((batch_size, pc_dim), generator=generator)
    return tokens, pcs, target, mask


def clear_cuda() -> None:
    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
        torch.cuda.reset_peak_memory_stats()


def one_step(model: nn.Module, optimizer: torch.optim.Optimizer, length: int, batch_size: int, device: torch.device, seed: int) -> tuple[float, int]:
    tokens, pcs, target, mask = make_batch(length, batch_size, model.cfg.pc_dim, seed)
    tokens = tokens.to(device, non_blocking=True)
    pcs = pcs.to(device, non_blocking=True)
    target = target.to(device, non_blocking=True)
    mask = mask.to(device, non_blocking=True)
    optimizer.zero_grad(set_to_none=True)
    autocast_enabled = device.type == "cuda"
    with torch.autocast(device_type=device.type, dtype=torch.bfloat16, enabled=autocast_enabled):
        logits = model(tokens, pcs)
        loss = F.cross_entropy(logits[mask], target[mask])
    if not bool(torch.isfinite(loss).item()):
        raise RuntimeError("non-finite loss")
    loss.backward()
    finite_gradients = all(
        parameter.grad is None or bool(torch.isfinite(parameter.grad).all().item())
        for parameter in model.parameters()
    )
    if not finite_gradients:
        raise RuntimeError("non-finite gradients")
    optimizer.step()
    return float(loss.detach().cpu()), int(mask.sum().item())


def fits_one_step(kind: str, cfg: Architecture, length: int, batch_size: int, device: torch.device) -> tuple[bool, str | None]:
    clear_cuda()
    try:
        model = MaskedGenotypeModel(kind, cfg).to(device)
        optimizer = torch.optim.AdamW(model.parameters(), lr=1e-4)
        one_step(model, optimizer, length, batch_size, device, seed=9917)
        if device.type == "cuda":
            torch.cuda.synchronize()
        del optimizer, model
        clear_cuda()
        return True, None
    except torch.OutOfMemoryError:
        clear_cuda()
        return False, "cuda_oom"
    except RuntimeError as exc:
        if "out of memory" in str(exc).lower():
            clear_cuda()
            return False, "runtime_oom"
        raise


def find_oom_boundary(kind: str, cfg: Architecture, base_length: int, limit: int, batch_size: int, device: torch.device) -> dict[str, int | None]:
    if device.type != "cuda":
        return {"largest_tested_fit": base_length, "first_tested_oom": None, "search_limit": limit}
    largest_fit: int | None = None
    first_oom: int | None = None
    length = base_length
    while length <= limit:
        fits, _ = fits_one_step(kind, cfg, length, batch_size, device)
        if not fits:
            first_oom = length
            break
        largest_fit = length
        length *= 2
    return {"largest_tested_fit": largest_fit, "first_tested_oom": first_oom, "search_limit": limit}


def profile_cell(kind: str, cfg: Architecture, length: int, steps: int, batch_size: int, device: torch.device, output_dir: Path) -> dict[str, object]:
    clear_cuda()
    model = MaskedGenotypeModel(kind, cfg).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-4)
    losses: list[float] = []
    masked_tokens = 0

    # Warm-up is excluded from throughput but recorded by the synchronized peak allocator.
    one_step(model, optimizer, length, batch_size, device, seed=1000)
    if device.type == "cuda":
        torch.cuda.synchronize()
        torch.cuda.reset_peak_memory_stats()
    started = time.perf_counter()
    for step in range(steps):
        loss, count = one_step(model, optimizer, length, batch_size, device, seed=2000 + step)
        losses.append(loss)
        masked_tokens += count
    if device.type == "cuda":
        torch.cuda.synchronize()
    training_seconds = time.perf_counter() - started

    checkpoint = output_dir / f"checkpoint_{kind}_L{length}.pt"
    checkpoint_started = time.perf_counter()
    torch.save({"model": model.state_dict(), "optimizer": optimizer.state_dict(), "steps": steps}, checkpoint)
    checkpoint_seconds = time.perf_counter() - checkpoint_started
    allocated = torch.cuda.max_memory_allocated() if device.type == "cuda" else None
    reserved = torch.cuda.max_memory_reserved() if device.type == "cuda" else None
    result: dict[str, object] = {
        "model": kind,
        "sequence_length": length,
        "batch_size": batch_size,
        "steps": steps,
        "parameter_count": parameter_count(model),
        "effective_ff_dim": model.effective_ff_dim,
        "mean_loss": sum(losses) / len(losses),
        "training_seconds": training_seconds,
        "seconds_per_step": training_seconds / steps,
        "masked_tokens": masked_tokens,
        "masked_tokens_per_second": masked_tokens / training_seconds,
        "checkpoint_seconds": checkpoint_seconds,
        "checkpoint_bytes": checkpoint.stat().st_size,
        "peak_allocated_bytes": allocated,
        "peak_reserved_bytes": reserved,
        "status": "PASS",
    }
    del optimizer, model
    clear_cuda()
    return result


def parse_lengths(value: str) -> list[int]:
    lengths = sorted({int(item) for item in value.split(",") if item.strip()})
    if not lengths or any(length <= 0 for length in lengths):
        raise argparse.ArgumentTypeError("lengths must be positive comma-separated integers")
    return lengths


def main() -> None:
    parser = argparse.ArgumentParser(description="T03 synchronized end-to-end profiler for the three frozen V7 model families.")
    parser.add_argument("--models", default="local_attn,local_attn_gpc,grassmann_full")
    parser.add_argument("--lengths", type=parse_lengths, default=parse_lengths("8192,131072,262144"))
    parser.add_argument("--steps", type=int, default=100)
    parser.add_argument("--batch-size", type=int, default=1)
    parser.add_argument("--oom-search-limit", type=int, default=524288)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--allow-cpu-dry-run", action="store_true")
    args = parser.parse_args()
    if args.steps <= 0 or args.batch_size <= 0:
        raise SystemExit("steps and batch size must be positive")

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if device.type != "cuda" and not args.allow_cpu_dry_run:
        raise SystemExit("CUDA is required for a valid T03 profile; use --allow-cpu-dry-run only for code checks")
    if device.type == "cuda":
        props = torch.cuda.get_device_properties(0)
        if (props.major, props.minor) < (12, 0):
            raise SystemExit(f"T03 requires sm_120 or newer, found sm_{props.major}{props.minor}")

    models = [item.strip() for item in args.models.split(",") if item.strip()]
    cfg = Architecture()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    boundaries: dict[str, dict[str, int | None]] = {}
    effective_lengths = args.lengths
    if device.type == "cpu":
        effective_lengths = [min(args.lengths[0], 128)]

    for kind in models:
        boundaries[kind] = find_oom_boundary(
            kind, cfg, min(effective_lengths), args.oom_search_limit,
            args.batch_size, device,
        )
        largest_fit = boundaries[kind]["largest_tested_fit"]
        selected = list(effective_lengths)
        if isinstance(largest_fit, int) and largest_fit not in selected:
            selected.append(largest_fit)
        for length in sorted(set(selected)):
            if isinstance(largest_fit, int) and length > largest_fit:
                rows.append({"model": kind, "sequence_length": length, "status": "SKIPPED_ABOVE_OOM_BOUNDARY"})
                continue
            rows.append(profile_cell(kind, cfg, length, args.steps, args.batch_size, device, args.output_dir))

    payload = {
        "schema_version": "1.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "valid_t03_measurement": device.type == "cuda" and args.steps == 100,
        "device_type": device.type,
        "cuda_visible_devices_set": "CUDA_VISIBLE_DEVICES" in os.environ,
        "torch_version": torch.__version__,
        "torch_cuda_runtime": torch.version.cuda,
        "architecture": asdict(cfg),
        "oom_boundaries": boundaries,
        "profiles": rows,
    }
    output = args.output_dir / "PROFILE_REPORT.json"
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(output), "valid_t03_measurement": payload["valid_t03_measurement"], "profiles": len(rows)}, indent=2))


if __name__ == "__main__":
    main()
