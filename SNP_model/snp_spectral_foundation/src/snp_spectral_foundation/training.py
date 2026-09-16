from __future__ import annotations

import torch
from torch import Tensor
import torch.nn.functional as F


def sample_ssl_mask(genotype: Tensor, probability: float, generator: torch.Generator | None = None) -> Tensor:
    """Sample masks only over observed genotypes."""
    if not 0.0 < probability < 1.0:
        raise ValueError("probability must be between zero and one")
    random = torch.rand(genotype.shape, device=genotype.device, generator=generator)
    return (random < probability) & (genotype >= 0)


def masked_genotype_loss(logits: Tensor, genotype: Tensor, mask: Tensor) -> Tensor:
    valid = mask.bool() & (genotype >= 0)
    if not bool(valid.any()):
        raise ValueError("mask contains no observed target")
    return F.cross_entropy(logits[valid], genotype.long()[valid])


@torch.no_grad()
def contextual_gate_metrics(
    logits: Tensor, genotype: Tensor, mask: Tensor, allele_frequency: Tensor
) -> dict[str, float]:
    """Compare SSL accuracy with the train-only Hardy-Weinberg modal baseline."""
    valid = mask.bool() & (genotype >= 0)
    if allele_frequency.ndim == 2:
        allele_frequency = allele_frequency.unsqueeze(0).expand_as(genotype)
    model_call = logits.argmax(-1)
    p = allele_frequency.clamp(0.0, 1.0)
    hwe = torch.stack(((1 - p).square(), 2 * p * (1 - p), p.square()), dim=-1)
    marginal_call = hwe.argmax(-1)
    target = genotype.long()
    model_accuracy = (model_call[valid] == target[valid]).float().mean().item()
    marginal_accuracy = (marginal_call[valid] == target[valid]).float().mean().item()
    return {
        "model_accuracy": model_accuracy,
        "marginal_accuracy": marginal_accuracy,
        "contextual_lift": model_accuracy - marginal_accuracy,
        "n_masked": int(valid.sum().item()),
    }
