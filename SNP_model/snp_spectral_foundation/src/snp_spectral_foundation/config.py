from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ModelConfig:
    """Configuration for the executable research prototype.

    ``max_blocks`` and ``snps_per_block`` describe the padded tensor layout. A
    production reader may stream any number of chromosomes through this same
    block interface.
    """

    # Base-scale defaults: 16,384 x 64 covers just over one million variants.
    max_blocks: int = 16_384
    snps_per_block: int = 64
    d_model: int = 128
    n_heads: int = 8
    local_layers: int = 4
    local_ff_dim: int = 512
    spectral_rank: int = 8
    spectral_oversample: int = 4
    spectral_power_iterations: int = 1
    grassmann_probes: int = 24
    memory_dim: int = 128
    global_slots: int = 64
    global_layers: int = 4
    num_traits: int = 256
    trait_dim: int = 64
    effect_dim: int = 64
    output_dim: int = 1
    dropout: float = 0.1
    spectral_jitter: float = 1e-5
    # ``legacy_lookup`` preserves compatibility with frozen historical
    # checkpoints. E0-X must use ``relative_continuous`` together with
    # ``use_local_block_embedding=False`` so unseen loci have no lookup-only
    # parameters.
    local_position_mode: str = "legacy_lookup"
    use_local_block_embedding: bool = True
    continuous_position_scale_bp: float = 1_000_000.0

    def __post_init__(self) -> None:
        if self.d_model % self.n_heads:
            raise ValueError("d_model must be divisible by n_heads")
        if not 1 <= self.spectral_rank <= min(self.d_model, self.snps_per_block - 1):
            raise ValueError("spectral_rank must not exceed d_model or snps_per_block - 1")
        if self.spectral_oversample < 0:
            raise ValueError("spectral_oversample must be non-negative")
        if self.spectral_rank + self.spectral_oversample > self.d_model:
            raise ValueError("spectral rank plus oversampling must not exceed d_model")
        if self.spectral_power_iterations < 0:
            raise ValueError("spectral_power_iterations must be non-negative")
        if self.max_blocks < 1 or self.snps_per_block < 2:
            raise ValueError("at least one block and two SNPs per block are required")
        if self.global_slots < 1 or self.num_traits < 1:
            raise ValueError("global_slots and num_traits must be positive")
        if self.local_position_mode not in {"legacy_lookup", "relative_continuous"}:
            raise ValueError("local_position_mode must be legacy_lookup or relative_continuous")
        if self.continuous_position_scale_bp <= 0:
            raise ValueError("continuous_position_scale_bp must be positive")
