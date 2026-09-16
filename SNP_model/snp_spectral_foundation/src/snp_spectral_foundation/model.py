from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import torch
from torch import Tensor, nn
import torch.nn.functional as F

from .config import ModelConfig


MISSING_TOKEN = 3
MASK_TOKEN = 4
MemoryMode = Literal["hybrid", "spectral_only", "residual_only", "aligned_control"]


def encode_genotypes(genotype: Tensor, masked_positions: Tensor | None = None) -> Tensor:
    """Map dosage {-1,0,1,2} to model tokens and optionally apply SSL masks."""
    if genotype.ndim != 3:
        raise ValueError("genotype must have shape [batch, blocks, snps]")
    invalid = (genotype < -1) | (genotype > 2)
    if bool(invalid.any()):
        raise ValueError("genotype entries must be -1 (missing), 0, 1, or 2")
    tokens = genotype.long().clone()
    tokens[tokens < 0] = MISSING_TOKEN
    if masked_positions is not None:
        if masked_positions.shape != genotype.shape:
            raise ValueError("masked_positions must match genotype shape")
        tokens[masked_positions.bool() & (genotype >= 0)] = MASK_TOKEN
    return tokens


def grassmann_probe_features(u: Tensor, probes: Tensor) -> Tensor:
    """Compact basis-invariant measurements of the projector ``U U^T``.

    Each output is ``a_i^T U U^T a_i``. Consequently ``U`` and ``UQ`` give
    the same result for every orthogonal rank-by-rank matrix ``Q``.
    """
    if u.ndim != 4 or probes.ndim != 2 or u.shape[-2] != probes.shape[0]:
        raise ValueError("expected U [..., d_model, rank] and probes [d_model, n]")
    coordinates = torch.einsum("bndr,dp->bnrp", u, probes)
    return coordinates.square().sum(dim=-2)


@dataclass
class SpectralState:
    basis: Tensor
    eigenvalues: Tensor
    hidden_residual: Tensor


class LocalGenotypeEncoder(nn.Module):
    """Shared short-range encoder used for masked-genotype pretraining."""

    def __init__(self, cfg: ModelConfig) -> None:
        super().__init__()
        self.cfg = cfg
        self.genotype_embedding = nn.Embedding(5, cfg.d_model)
        if cfg.local_position_mode == "legacy_lookup":
            self.position_embedding: nn.Module | None = nn.Embedding(cfg.snps_per_block, cfg.d_model)
            self.continuous_position_projection: nn.Module | None = None
        else:
            self.position_embedding = None
            self.continuous_position_projection = nn.Sequential(
                nn.Linear(3, cfg.d_model), nn.GELU(), nn.Linear(cfg.d_model, cfg.d_model)
            )
        self.block_embedding: nn.Module | None = (
            nn.Embedding(cfg.max_blocks, cfg.d_model) if cfg.use_local_block_embedding else None
        )
        self.frequency_projection = nn.Sequential(
            nn.Linear(2, cfg.d_model), nn.GELU(), nn.Linear(cfg.d_model, cfg.d_model)
        )
        layer = nn.TransformerEncoderLayer(
            d_model=cfg.d_model,
            nhead=cfg.n_heads,
            dim_feedforward=cfg.local_ff_dim,
            dropout=cfg.dropout,
            activation="gelu",
            batch_first=True,
            norm_first=True,
        )
        self.encoder = nn.TransformerEncoder(layer, cfg.local_layers)
        self.norm = nn.LayerNorm(cfg.d_model)
        self.genotype_head = nn.Linear(cfg.d_model, 3)

    def forward(
        self,
        tokens: Tensor,
        allele_frequency: Tensor,
        block_ids: Tensor | None = None,
        *,
        genomic_positions: Tensor | None = None,
    ) -> tuple[Tensor, Tensor]:
        batch, blocks, snps = tokens.shape
        if snps != self.cfg.snps_per_block:
            raise ValueError(f"expected {self.cfg.snps_per_block} SNPs per block, got {snps}")
        if allele_frequency.shape not in ((blocks, snps), (batch, blocks, snps)):
            raise ValueError("allele_frequency must have shape [blocks, snps] or [batch, blocks, snps]")
        if bool(((allele_frequency < 0) | (allele_frequency > 1)).any()):
            raise ValueError("allele_frequency entries must be in [0, 1]")
        if allele_frequency.ndim == 2:
            allele_frequency = allele_frequency.unsqueeze(0).expand(batch, -1, -1)
        minor_frequency = torch.minimum(allele_frequency, 1.0 - allele_frequency)
        frequency_features = torch.stack((allele_frequency, minor_frequency), dim=-1)
        x = self.genotype_embedding(tokens)
        if self.cfg.local_position_mode == "legacy_lookup":
            if genomic_positions is not None:
                raise ValueError("legacy_lookup does not accept genomic_positions")
            positions = torch.arange(snps, device=tokens.device)
            if self.position_embedding is None:
                raise RuntimeError("legacy position embedding is missing")
            x = x + self.position_embedding(positions)[None, None, :, :]
        else:
            if genomic_positions is None:
                raise ValueError("relative_continuous requires genomic_positions")
            if genomic_positions.shape == (blocks, snps):
                genomic_positions = genomic_positions.unsqueeze(0).expand(batch, -1, -1)
            if genomic_positions.shape != (batch, blocks, snps):
                raise ValueError(
                    "genomic_positions must have shape [blocks, snps] or [batch, blocks, snps]"
                )
            if torch.is_floating_point(genomic_positions):
                raise ValueError("genomic_positions must use an integer dtype")
            # Genomic bp coordinates can be hundreds of millions. Subtract in
            # integer space before converting differences to the model dtype;
            # absolute float32 coordinates would quantize short marker gaps.
            integer_positions = genomic_positions.to(device=tokens.device, dtype=torch.int64)
            gaps_integer = integer_positions[..., 1:] - integer_positions[..., :-1]
            if bool((gaps_integer < 0).any()):
                raise ValueError("genomic_positions must be nondecreasing within each block")
            span_integer = (
                integer_positions[..., -1:] - integer_positions[..., :1]
            ).clamp_min(1)
            offsets = (integer_positions - integer_positions[..., :1]).to(x.dtype)
            span = span_integer.to(x.dtype)
            relative = offsets / span
            leading_zero = torch.zeros_like(integer_positions[..., :1])
            gaps = torch.cat((leading_zero, gaps_integer), dim=-1).to(x.dtype)
            scale = torch.as_tensor(
                self.cfg.continuous_position_scale_bp, device=tokens.device, dtype=x.dtype
            )
            log_denominator = torch.log1p(scale)
            log_gap = torch.log1p(gaps) / log_denominator
            log_span = (torch.log1p(span) / log_denominator).expand_as(relative)
            position_features = torch.stack((relative, log_gap, log_span), dim=-1)
            if self.continuous_position_projection is None:
                raise RuntimeError("continuous position projection is missing")
            x = x + self.continuous_position_projection(position_features)

        if self.block_embedding is not None:
            if block_ids is None:
                block_ids = torch.arange(blocks, device=tokens.device)[None, :].expand(batch, -1)
            elif block_ids.ndim == 1:
                block_ids = block_ids[None, :].expand(batch, -1)
            if block_ids.shape != (batch, blocks):
                raise ValueError("block_ids must have shape [blocks] or [batch, blocks]")
            if bool(((block_ids < 0) | (block_ids >= self.cfg.max_blocks)).any()):
                raise ValueError("block_ids are outside the configured vocabulary")
            x = x + self.block_embedding(block_ids)[:, :, None, :]
        elif block_ids is not None:
            raise ValueError("block_ids are prohibited when use_local_block_embedding is false")
        x = x + self.frequency_projection(frequency_features.to(x.dtype))
        x = x.reshape(batch * blocks, snps, self.cfg.d_model)
        hidden = self.norm(self.encoder(x)).reshape(batch, blocks, snps, self.cfg.d_model)
        return hidden, self.genotype_head(hidden)


class BlockSpectralMemory(nn.Module):
    """Convert local SNP states into spectral, residual, and control memories."""

    def __init__(self, cfg: ModelConfig) -> None:
        super().__init__()
        self.cfg = cfg
        probes = torch.randn(cfg.d_model, cfg.grassmann_probes)
        probes = F.normalize(probes, dim=0)
        self.register_buffer("grassmann_probes", probes)
        sketch_width = cfg.spectral_rank + cfg.spectral_oversample
        sketch = torch.randn(cfg.d_model, sketch_width)
        sketch, _ = torch.linalg.qr(sketch, mode="reduced")
        self.register_buffer("spectral_sketch", sketch)

        spectral_input = cfg.grassmann_probes + cfg.spectral_rank + 4
        self.spectral_projector = nn.Sequential(
            nn.Linear(spectral_input, cfg.memory_dim), nn.GELU(), nn.LayerNorm(cfg.memory_dim)
        )
        # This independent sidecar is position-aware and does not quotient away
        # signed individual deviations from the train-only allele frequency.
        self.residual_projector = nn.Sequential(
            nn.Linear(2 * cfg.snps_per_block, cfg.memory_dim),
            nn.GELU(),
            nn.LayerNorm(cfg.memory_dim),
        )
        # A strong orientation-preserving control demanded by the previous gates.
        self.aligned_projector = nn.Sequential(
            nn.Linear(cfg.snps_per_block, cfg.memory_dim),
            nn.GELU(),
            nn.LayerNorm(cfg.memory_dim),
        )
        self.hybrid_fusion = nn.Sequential(
            nn.Linear(2 * cfg.memory_dim, cfg.memory_dim),
            nn.GELU(),
            nn.LayerNorm(cfg.memory_dim),
        )

    def decompose(self, hidden: Tensor, observed: Tensor) -> SpectralState:
        weights = observed.to(hidden.dtype).unsqueeze(-1)
        count = weights.sum(dim=-2, keepdim=True).clamp_min(1.0)
        mean = (hidden * weights).sum(dim=-2, keepdim=True) / count
        centered = (hidden - mean) * weights
        denominator = (count.squeeze(-2).unsqueeze(-1) - 1.0).clamp_min(1.0)

        def covariance_multiply(vectors: Tensor) -> Tensor:
            token_coordinates = torch.einsum("bnsd,bndk->bnsk", centered, vectors)
            return torch.einsum("bnsd,bnsk->bndk", centered, token_coordinates) / denominator

        # Randomized subspace iteration avoids materializing/eigendecomposing a
        # d_model x d_model covariance for every person-block. The expensive
        # operations are linear in d_model, SNP count, and requested rank.
        sketch = self.spectral_sketch.to(dtype=hidden.dtype)
        sketch = sketch[None, None, :, :].expand(hidden.shape[0], hidden.shape[1], -1, -1)
        candidate, _ = torch.linalg.qr(covariance_multiply(sketch), mode="reduced")
        for _ in range(self.cfg.spectral_power_iterations):
            candidate, _ = torch.linalg.qr(covariance_multiply(candidate), mode="reduced")
        covariance_candidate = covariance_multiply(candidate)
        small = torch.einsum("bndk,bndl->bnkl", candidate, covariance_candidate)
        eye = torch.eye(small.shape[-1], device=hidden.device, dtype=hidden.dtype)
        values, rotation = torch.linalg.eigh(small + self.cfg.spectral_jitter * eye)
        values = values[..., -self.cfg.spectral_rank :].flip(-1).clamp_min(0.0)
        rotation = rotation[..., -self.cfg.spectral_rank :].flip(-1)
        basis = torch.einsum("bndk,bnkr->bndr", candidate, rotation)
        projected = torch.einsum("bnsd,bndr->bnsr", centered, basis)
        reconstruction = torch.einsum("bnsr,bndr->bnsd", projected, basis)
        return SpectralState(basis=basis, eigenvalues=values, hidden_residual=centered - reconstruction)

    @staticmethod
    def _maf_summary(allele_frequency: Tensor) -> Tensor:
        maf = torch.minimum(allele_frequency, 1.0 - allele_frequency)
        return torch.stack(
            (maf.mean(-1), maf.std(-1, unbiased=False), maf.amin(-1), maf.amax(-1)), dim=-1
        )

    def forward(
        self,
        hidden: Tensor,
        genotype: Tensor,
        allele_frequency: Tensor,
        mode: MemoryMode = "hybrid",
    ) -> dict[str, Tensor]:
        batch, blocks, snps, _ = hidden.shape
        if allele_frequency.ndim == 2:
            allele_frequency = allele_frequency.unsqueeze(0).expand(batch, -1, -1)
        observed = genotype >= 0
        state = self.decompose(hidden, observed)
        grass = grassmann_probe_features(state.basis, self.grassmann_probes)
        spectrum = torch.log1p(state.eigenvalues)
        spectrum = spectrum / spectrum.sum(dim=-1, keepdim=True).clamp_min(1e-6)
        spectral = self.spectral_projector(
            torch.cat((grass, spectrum, self._maf_summary(allele_frequency)), -1)
        )

        dosage = genotype.clamp_min(0).to(hidden.dtype)
        # Residualization needs oriented ALT-allele frequency, not MAF. MAF alone
        # cannot distinguish an ALT-major site from a REF-major site.
        dosage_residual = (dosage - 2.0 * allele_frequency) * observed.to(hidden.dtype)
        sidecar_input = torch.cat((dosage_residual, (~observed).to(hidden.dtype)), dim=-1)
        residual = self.residual_projector(sidecar_input)
        aligned = self.aligned_projector(dosage_residual)

        if mode == "hybrid":
            memory = self.hybrid_fusion(torch.cat((spectral, residual), -1))
        elif mode == "spectral_only":
            memory = spectral
        elif mode == "residual_only":
            memory = residual
        elif mode == "aligned_control":
            memory = aligned
        else:
            raise ValueError(f"unknown memory mode: {mode}")
        return {
            "memory": memory,
            "spectral": spectral,
            "residual": residual,
            "aligned": aligned,
            "basis": state.basis,
            "eigenvalues": state.eigenvalues,
            "hidden_residual": state.hidden_residual,
        }


class ScalableGlobalMixer(nn.Module):
    """Linear-in-block-count latent-slot mixer instead of dense block attention."""

    def __init__(self, cfg: ModelConfig) -> None:
        super().__init__()
        self.cfg = cfg
        self.block_positions = nn.Embedding(cfg.max_blocks, cfg.memory_dim)
        self.slots = nn.Parameter(torch.randn(cfg.global_slots, cfg.memory_dim) * 0.02)
        self.slots_read_blocks = nn.MultiheadAttention(
            cfg.memory_dim, cfg.n_heads, dropout=cfg.dropout, batch_first=True
        )
        layer = nn.TransformerEncoderLayer(
            cfg.memory_dim,
            cfg.n_heads,
            dim_feedforward=2 * cfg.memory_dim,
            dropout=cfg.dropout,
            activation="gelu",
            batch_first=True,
            norm_first=True,
        )
        self.slot_mixer = nn.TransformerEncoder(layer, cfg.global_layers)
        self.blocks_read_slots = nn.MultiheadAttention(
            cfg.memory_dim, cfg.n_heads, dropout=cfg.dropout, batch_first=True
        )
        self.norm = nn.LayerNorm(cfg.memory_dim)

    def forward(
        self,
        block_memory: Tensor,
        block_mask: Tensor | None = None,
        block_ids: Tensor | None = None,
    ) -> tuple[Tensor, Tensor]:
        batch, blocks, _ = block_memory.shape
        if blocks > self.cfg.max_blocks:
            raise ValueError(f"received {blocks} blocks; max_blocks={self.cfg.max_blocks}")
        if block_ids is None:
            block_ids = torch.arange(blocks, device=block_memory.device)[None, :].expand(batch, -1)
        elif block_ids.ndim == 1:
            block_ids = block_ids[None, :].expand(batch, -1)
        if block_ids.shape != (batch, blocks):
            raise ValueError("block_ids must have shape [blocks] or [batch, blocks]")
        x = block_memory + self.block_positions(block_ids)
        slots = self.slots.unsqueeze(0).expand(batch, -1, -1)
        slots, _ = self.slots_read_blocks(slots, x, x, key_padding_mask=block_mask)
        slots = self.slot_mixer(slots)
        global_context, _ = self.blocks_read_slots(x, slots, slots)
        return self.norm(x + global_context), slots


class SNPFoundationModel(nn.Module):
    """End-to-end implementation of the proposed SNP architecture.

    The phenotype condition is a *trait identity*, never an observed outcome.
    This lets one encoder serve multiple endpoints without label leakage.
    """

    def __init__(self, cfg: ModelConfig) -> None:
        super().__init__()
        self.cfg = cfg
        self.local_encoder = LocalGenotypeEncoder(cfg)
        self.block_memory = BlockSpectralMemory(cfg)
        self.global_mixer = ScalableGlobalMixer(cfg)
        self.trait_embedding = nn.Embedding(cfg.num_traits, cfg.trait_dim)
        self.trait_query = nn.Linear(cfg.trait_dim, cfg.memory_dim)
        self.effect_projector = nn.Sequential(
            nn.Linear(cfg.memory_dim + cfg.trait_dim, cfg.effect_dim),
            nn.GELU(),
            nn.LayerNorm(cfg.effect_dim),
        )
        self.residual_to_effect = nn.Linear(cfg.memory_dim, cfg.effect_dim)
        self.fusion_gate = nn.Linear(2 * cfg.effect_dim + cfg.trait_dim, cfg.effect_dim)
        self.prediction_head = nn.Sequential(
            nn.Linear(cfg.effect_dim, cfg.effect_dim), nn.GELU(), nn.Linear(cfg.effect_dim, cfg.output_dim)
        )

    def forward_pretrain(
        self,
        genotype: Tensor,
        allele_frequency: Tensor,
        masked_positions: Tensor,
        block_ids: Tensor | None = None,
        *,
        genomic_positions: Tensor | None = None,
    ) -> dict[str, Tensor]:
        tokens = encode_genotypes(genotype, masked_positions)
        hidden, logits = self.local_encoder(
            tokens, allele_frequency, block_ids, genomic_positions=genomic_positions
        )
        return {"hidden": hidden, "logits": logits, "mask": masked_positions.bool() & (genotype >= 0)}

    def forward(
        self,
        genotype: Tensor,
        allele_frequency: Tensor,
        trait_id: Tensor,
        *,
        memory_mode: MemoryMode = "hybrid",
        block_mask: Tensor | None = None,
        use_global_mixer: bool = True,
        block_ids: Tensor | None = None,
        genomic_positions: Tensor | None = None,
    ) -> dict[str, Tensor]:
        batch = genotype.shape[0]
        if trait_id.shape != (batch,):
            raise ValueError("trait_id must have shape [batch]")
        hidden, pretrain_logits = self.local_encoder(
            encode_genotypes(genotype),
            allele_frequency,
            block_ids,
            genomic_positions=genomic_positions,
        )
        memory = self.block_memory(hidden, genotype, allele_frequency, memory_mode)
        if block_mask is not None and bool(block_mask.all(dim=-1).any()):
            raise ValueError("each person must retain at least one block")
        if use_global_mixer:
            mixed_blocks, global_slots = self.global_mixer(memory["memory"], block_mask, block_ids)
        else:
            mixed_blocks = memory["memory"]
            global_slots = memory["memory"].new_empty(batch, 0, self.cfg.memory_dim)

        trait = self.trait_embedding(trait_id)
        query = self.trait_query(trait)
        scores = torch.einsum("bnd,bd->bn", mixed_blocks, query) / self.cfg.memory_dim**0.5
        if block_mask is not None:
            scores = scores.masked_fill(block_mask, torch.finfo(scores.dtype).min)
        weights = scores.softmax(dim=-1)
        spectral_summary = torch.einsum("bn,bnd->bd", weights, mixed_blocks)
        residual_summary = torch.einsum("bn,bnd->bd", weights, memory["residual"])
        effect = self.effect_projector(torch.cat((spectral_summary, trait), -1))
        residual_effect = self.residual_to_effect(residual_summary)
        gate = torch.sigmoid(self.fusion_gate(torch.cat((effect, residual_effect, trait), -1)))
        fused = gate * effect + (1.0 - gate) * residual_effect
        prediction = self.prediction_head(fused)
        return {
            "prediction": prediction,
            "effect": effect,
            "fused_effect": fused,
            "block_weights": weights,
            "global_slots": global_slots,
            "pretrain_logits": pretrain_logits,
            **memory,
        }
