from __future__ import annotations

import json
from pathlib import Path
import sys

import torch
import torch.nn.functional as F

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from snp_spectral_foundation import ModelConfig, SNPFoundationModel  # noqa: E402
from snp_spectral_foundation.training import (  # noqa: E402
    contextual_gate_metrics,
    masked_genotype_loss,
    sample_ssl_mask,
)


def synthetic_genotypes(n: int, blocks: int, snps: int, seed: int) -> torch.Tensor:
    """Small LD-correlated smoke fixture; not scientific evidence."""
    generator = torch.Generator().manual_seed(seed)
    latent = torch.randint(0, 3, (n, blocks, 1), generator=generator)
    noise_mask = torch.rand((n, blocks, snps), generator=generator) < 0.18
    noise = torch.randint(0, 3, (n, blocks, snps), generator=generator)
    return torch.where(noise_mask, noise, latent.expand(-1, -1, snps))


def main() -> None:
    torch.manual_seed(17)
    cfg = ModelConfig(
        max_blocks=4,
        snps_per_block=8,
        d_model=16,
        n_heads=4,
        local_layers=1,
        local_ff_dim=32,
        spectral_rank=3,
        grassmann_probes=6,
        memory_dim=16,
        global_slots=2,
        num_traits=2,
        trait_dim=8,
        effect_dim=8,
        dropout=0.0,
    )
    genotype = synthetic_genotypes(24, 4, 8, 11)
    train = genotype[:18]
    validation = genotype[18:]
    # Train-only MAF. Validation information is never used to construct it.
    allele_frequency = train.float().mean(0) / 2.0
    model = SNPFoundationModel(cfg)

    optimizer = torch.optim.AdamW(model.local_encoder.parameters(), lr=2e-3)
    model.train()
    for step in range(40):
        mask = sample_ssl_mask(train, 0.25, torch.Generator().manual_seed(100 + step))
        result = model.forward_pretrain(train, allele_frequency, mask)
        loss = masked_genotype_loss(result["logits"], train, result["mask"])
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

    model.eval()
    val_mask = sample_ssl_mask(validation, 0.25, torch.Generator().manual_seed(999))
    with torch.no_grad():
        pretrain = model.forward_pretrain(validation, allele_frequency, val_mask)
        gate = contextual_gate_metrics(
            pretrain["logits"], validation, val_mask, allele_frequency
        )

    # Demonstrate the supervised path without claiming a valid phenotype result.
    trait_id = torch.zeros(len(train), dtype=torch.long)
    phenotype = (train[:, 0, :].float().mean(-1) + 0.2 * train[:, 3, 0].float()).unsqueeze(-1)
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3)
    model.train()
    supervised_loss = torch.tensor(float("nan"))
    for _ in range(3):
        prediction = model(train, allele_frequency, trait_id)["prediction"]
        supervised_loss = F.mse_loss(prediction, phenotype)
        optimizer.zero_grad()
        supervised_loss.backward()
        optimizer.step()

    output = {
        "status": "SMOKE_ONLY_NOT_SCIENTIFIC_EVIDENCE",
        "contextual_gate": gate,
        "supervised_loss_after_3_steps": float(supervised_loss.detach()),
        "parameters": sum(parameter.numel() for parameter in model.parameters()),
    }
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
