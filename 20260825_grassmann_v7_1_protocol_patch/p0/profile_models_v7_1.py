from __future__ import annotations

import argparse
import gc
import json
import os
import time
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path

import torch
from torch import nn
from torch.nn import functional as F


MODEL_KINDS = (
    "local_attn_8m_w256",
    "local_attn_gpc_8m_w256",
    "grassmann_full_8m_w256",
)


@dataclass(frozen=True)
class Architecture:
    model_dim: int = 448
    ff_dim: int = 1280
    layers: int = 4
    heads: int = 8
    attention_window: int = 256
    query_block: int = 128
    pc_dim: int = 16
    grassmann_rank: int = 32
    genotype_states: int = 3

    def validate(self) -> None:
        if self.attention_window != 256:
            raise ValueError("v7.1.0 freezes attention_window=256")
        if self.query_block <= 0 or self.attention_window < self.query_block:
            raise ValueError("query_block must be positive and <= attention_window")
        if (self.attention_window - self.query_block) % 2:
            raise ValueError("attention context must split evenly around query blocks")
        if self.model_dim % self.heads:
            raise ValueError("model_dim must be divisible by heads")


class BlockCenteredLocalAttention(nn.Module):
    """Each query attends to one block-centered context of exactly <=256 keys.

    This is an auditable bounded-local P0 implementation. It is not represented as
    token-by-token sliding-window equivalence to SNPBag.
    """

    def __init__(self, cfg: Architecture, ff_dim: int | None = None) -> None:
        super().__init__()
        self.window = cfg.attention_window
        self.query_block = cfg.query_block
        self.left_context = (self.window - self.query_block) // 2
        self.right_context = self.window - self.query_block - self.left_context
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
        tail = (-length) % self.query_block
        normalized = self.norm1(hidden)

        queries = F.pad(normalized, (0, 0, 0, tail))
        n_blocks = queries.shape[1] // self.query_block
        queries = queries.reshape(batch * n_blocks, self.query_block, width)

        context = F.pad(
            normalized,
            (0, 0, self.left_context, self.right_context + tail),
        )
        context = context.unfold(1, self.window, self.query_block)
        context = context.permute(0, 1, 3, 2).reshape(
            batch * n_blocks, self.window, width
        )

        valid = torch.ones(length, dtype=torch.bool, device=hidden.device)
        valid = F.pad(
            valid,
            (self.left_context, self.right_context + tail),
            value=False,
        )
        valid = valid.unfold(0, self.window, self.query_block)
        key_padding_mask = (~valid).unsqueeze(0).expand(batch, -1, -1)
        key_padding_mask = key_padding_mask.reshape(batch * n_blocks, self.window)

        mixed, _ = self.attn(
            queries,
            context,
            context,
            key_padding_mask=key_padding_mask,
            need_weights=False,
        )
        mixed = mixed.reshape(batch, n_blocks * self.query_block, width)[:, :length]
        hidden = hidden + mixed
        return hidden + self.ff(self.norm2(hidden))


class GrassmannBlockMixer(nn.Module):
    """Global block-summary exterior-product channel; token attention stays local."""

    def __init__(self, cfg: Architecture) -> None:
        super().__init__()
        self.block_size = cfg.attention_window
        self.reduce = nn.Linear(cfg.model_dim, cfg.grassmann_rank)
        i, j = torch.triu_indices(cfg.grassmann_rank, cfg.grassmann_rank, offset=1)
        self.register_buffer("wedge_i", i)
        self.register_buffer("wedge_j", j)
        self.project = nn.Linear(len(i), cfg.model_dim)
        self.gate = nn.Linear(2 * cfg.model_dim, cfg.model_dim)

    def forward(self, hidden: torch.Tensor) -> torch.Tensor:
        batch, length, width = hidden.shape
        tail = (-length) % self.block_size
        padded = F.pad(hidden, (0, 0, 0, tail)) if tail else hidden
        blocks = padded.reshape(batch, -1, self.block_size, width).mean(dim=2)
        reduced = self.reduce(blocks)
        context = reduced.mean(dim=1, keepdim=True)
        wedge = (
            reduced[..., self.wedge_i] * context[..., self.wedge_j]
            - reduced[..., self.wedge_j] * context[..., self.wedge_i]
        )
        wedge = F.normalize(wedge, p=2, dim=-1, eps=1e-8)
        geometry = self.project(wedge).repeat_interleave(self.block_size, dim=1)
        geometry = geometry[:, :length]
        alpha = torch.sigmoid(self.gate(torch.cat([hidden, geometry], dim=-1)))
        return alpha * hidden + (1.0 - alpha) * geometry


class MaskedGenotypeModel(nn.Module):
    def __init__(self, kind: str, cfg: Architecture) -> None:
        super().__init__()
        cfg.validate()
        if kind not in MODEL_KINDS:
            raise ValueError(kind)
        self.kind = kind
        self.cfg = cfg
        self.token = nn.Embedding(cfg.genotype_states + 1, cfg.model_dim)
        self.local_position = nn.Embedding(cfg.attention_window, cfg.model_dim)
        uses_pc = kind != "local_attn_8m_w256"
        self.pc_projection = nn.Linear(cfg.pc_dim, cfg.model_dim) if uses_pc else None
        self.grassmann = (
            GrassmannBlockMixer(cfg) if kind == "grassmann_full_8m_w256" else None
        )
        self.effective_ff_dim = (
            1104 if kind == "grassmann_full_8m_w256" else cfg.ff_dim
        )
        self.blocks = nn.ModuleList(
            BlockCenteredLocalAttention(cfg, ff_dim=self.effective_ff_dim)
            for _ in range(cfg.layers)
        )
        self.final_norm = nn.LayerNorm(cfg.model_dim)
        self.head = nn.Linear(cfg.model_dim, cfg.genotype_states)

    def forward(self, tokens: torch.Tensor, pcs: torch.Tensor) -> torch.Tensor:
        positions = torch.arange(tokens.shape[1], device=tokens.device)
        positions = positions % self.cfg.attention_window
        hidden = self.token(tokens) + self.local_position(positions)
        if self.pc_projection is not None:
            hidden = hidden + self.pc_projection(pcs).unsqueeze(1)
        if self.grassmann is not None:
            hidden = self.grassmann(hidden)
        for block in self.blocks:
            hidden = block(hidden)
        return self.head(self.final_norm(hidden))


def parameter_count(model: nn.Module) -> int:
    return sum(p.numel() for p in model.parameters() if p.requires_grad)


def make_batch(
    length: int,
    batch_size: int,
    pc_dim: int,
    seed: int,
    mask_rate: float,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
    generator = torch.Generator(device="cpu").manual_seed(seed)
    target = torch.randint(
        0, 3, (batch_size, length), generator=generator, dtype=torch.long
    )
    mask = torch.rand((batch_size, length), generator=generator) < mask_rate
    if not bool(mask.any()):
        mask[0, 0] = True
    tokens = target.clone()
    tokens[mask] = 3
    pcs = torch.randn((batch_size, pc_dim), generator=generator)
    return tokens, pcs, target, mask


def clear_device() -> None:
    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
        torch.cuda.reset_peak_memory_stats()


def one_step(
    model: MaskedGenotypeModel,
    optimizer: torch.optim.Optimizer,
    length: int,
    batch_size: int,
    device: torch.device,
    seed: int,
    mask_rate: float,
) -> tuple[float, int]:
    tokens, pcs, target, mask = make_batch(
        length, batch_size, model.cfg.pc_dim, seed, mask_rate
    )
    tokens = tokens.to(device, non_blocking=True)
    pcs = pcs.to(device, non_blocking=True)
    target = target.to(device, non_blocking=True)
    mask = mask.to(device, non_blocking=True)
    optimizer.zero_grad(set_to_none=True)
    with torch.autocast(
        device_type=device.type,
        dtype=torch.bfloat16,
        enabled=device.type == "cuda",
    ):
        logits = model(tokens, pcs)
        loss = F.cross_entropy(logits[mask], target[mask])
    if not bool(torch.isfinite(loss).item()):
        raise RuntimeError("non-finite loss")
    loss.backward()
    if not all(
        p.grad is None or bool(torch.isfinite(p.grad).all().item())
        for p in model.parameters()
    ):
        raise RuntimeError("non-finite gradients")
    optimizer.step()
    return float(loss.detach().cpu()), int(mask.sum().item())


def profile_cell(
    kind: str,
    cfg: Architecture,
    length: int,
    steps: int,
    batch_size: int,
    device: torch.device,
    output_dir: Path,
    mask_rate: float,
) -> dict[str, object]:
    clear_device()
    model = MaskedGenotypeModel(kind, cfg).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-4)
    losses: list[float] = []
    masked_tokens = 0
    one_step(model, optimizer, length, batch_size, device, 1000, mask_rate)
    if device.type == "cuda":
        torch.cuda.synchronize()
        torch.cuda.reset_peak_memory_stats()
    started = time.perf_counter()
    for step in range(steps):
        loss, count = one_step(
            model, optimizer, length, batch_size, device, 2000 + step, mask_rate
        )
        losses.append(loss)
        masked_tokens += count
    if device.type == "cuda":
        torch.cuda.synchronize()
    elapsed = time.perf_counter() - started
    checkpoint = output_dir / f"checkpoint_{kind}_L{length}.pt"
    checkpoint_started = time.perf_counter()
    torch.save(
        {"model": model.state_dict(), "optimizer": optimizer.state_dict(), "steps": steps},
        checkpoint,
    )
    checkpoint_seconds = time.perf_counter() - checkpoint_started
    row: dict[str, object] = {
        "model": kind,
        "sequence_length": length,
        "batch_size": batch_size,
        "steps": steps,
        "profile_mask_rate": mask_rate,
        "parameter_count": parameter_count(model),
        "effective_ff_dim": model.effective_ff_dim,
        "mean_loss": sum(losses) / len(losses),
        "training_seconds": elapsed,
        "seconds_per_step": elapsed / steps,
        "masked_tokens": masked_tokens,
        "masked_tokens_per_second": masked_tokens / elapsed,
        "checkpoint_seconds": checkpoint_seconds,
        "checkpoint_bytes": checkpoint.stat().st_size,
        "peak_allocated_bytes": (
            torch.cuda.max_memory_allocated() if device.type == "cuda" else None
        ),
        "peak_reserved_bytes": (
            torch.cuda.max_memory_reserved() if device.type == "cuda" else None
        ),
        "status": "PASS",
    }
    del optimizer, model
    clear_device()
    return row


def parse_lengths(value: str) -> list[int]:
    values = sorted({int(x.strip()) for x in value.split(",") if x.strip()})
    if not values or any(x <= 0 for x in values):
        raise argparse.ArgumentTypeError("lengths must be positive comma-separated integers")
    return values


def main() -> None:
    parser = argparse.ArgumentParser(description="V7.1.0 fixed-window T03 profiler")
    parser.add_argument("--models", default=",".join(MODEL_KINDS))
    parser.add_argument("--lengths", type=parse_lengths, required=True)
    parser.add_argument("--steps", type=int, default=100)
    parser.add_argument("--batch-size", type=int, default=1)
    parser.add_argument("--profile-mask-rate", type=float, default=0.90)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--allow-cpu-dry-run", action="store_true")
    args = parser.parse_args()
    if args.steps <= 0 or args.batch_size <= 0:
        raise SystemExit("steps and batch size must be positive")
    if not 0.0 < args.profile_mask_rate < 1.0:
        raise SystemExit("profile mask rate must be in (0,1)")

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if device.type != "cuda" and not args.allow_cpu_dry_run:
        raise SystemExit("CUDA is required for valid T03; CPU is a code-check only")
    if device.type == "cuda":
        props = torch.cuda.get_device_properties(0)
        if (props.major, props.minor) < (12, 0):
            raise SystemExit(f"T03 requires sm_120+, found sm_{props.major}{props.minor}")

    models = tuple(x.strip() for x in args.models.split(",") if x.strip())
    unknown = sorted(set(models) - set(MODEL_KINDS))
    if unknown:
        raise SystemExit(f"unknown models: {unknown}")
    cfg = Architecture()
    cfg.validate()
    args.output_dir.mkdir(parents=True, exist_ok=False)
    effective_lengths = args.lengths if device.type == "cuda" else [min(args.lengths[0], 256)]
    rows: list[dict[str, object]] = []
    for kind in models:
        for length in effective_lengths:
            try:
                rows.append(
                    profile_cell(
                        kind,
                        cfg,
                        length,
                        args.steps,
                        args.batch_size,
                        device,
                        args.output_dir,
                        args.profile_mask_rate,
                    )
                )
            except (torch.OutOfMemoryError, RuntimeError) as exc:
                if isinstance(exc, torch.OutOfMemoryError) or "out of memory" in str(exc).lower():
                    clear_device()
                    rows.append(
                        {
                            "model": kind,
                            "sequence_length": length,
                            "status": "OOM",
                            "error_type": type(exc).__name__,
                        }
                    )
                    continue
                raise

    visible = os.environ.get("CUDA_VISIBLE_DEVICES")
    payload = {
        "schema_version": "1.1",
        "protocol_version": "v7.1.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "valid_t03_measurement": (
            device.type == "cuda"
            and args.steps == 100
            and cfg.attention_window == 256
            and visible is not None
            and visible.strip() not in {"", "0"}
        ),
        "device_type": device.type,
        "cuda_visible_devices": visible,
        "logical_cuda_device": 0 if device.type == "cuda" else None,
        "torch_version": torch.__version__,
        "torch_cuda_runtime": torch.version.cuda,
        "attention_semantics": "block_centered_local_attention_w256_q128_not_exact_snpbag_sliding",
        "architecture": asdict(cfg),
        "profiles": rows,
    }
    output = args.output_dir / "PROFILE_REPORT.v7.1.0.json"
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(output), "valid_t03_measurement": payload["valid_t03_measurement"], "profiles": len(rows)}, indent=2))
    if any(row.get("status") != "PASS" for row in rows):
        raise SystemExit(4)


if __name__ == "__main__":
    main()
