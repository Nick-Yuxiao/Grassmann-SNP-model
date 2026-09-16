from __future__ import annotations

import unittest

import torch

from snp_spectral_foundation import ModelConfig, SNPFoundationModel
from snp_spectral_foundation.model import LocalGenotypeEncoder, encode_genotypes, grassmann_probe_features
from snp_spectral_foundation.training import masked_genotype_loss, sample_ssl_mask


class ModelTest(unittest.TestCase):
    def setUp(self) -> None:
        torch.manual_seed(7)
        self.cfg = ModelConfig(
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
            num_traits=3,
            trait_dim=8,
            effect_dim=8,
            dropout=0.0,
        )
        self.model = SNPFoundationModel(self.cfg)
        self.genotype = torch.randint(0, 3, (3, 4, 8))
        self.genotype[0, 0, 0] = -1
        self.allele_frequency = torch.rand(4, 8) * 0.9 + 0.05

    def test_end_to_end_modes_are_finite(self) -> None:
        for mode in ("hybrid", "spectral_only", "residual_only", "aligned_control"):
            out = self.model(
                self.genotype, self.allele_frequency, torch.tensor([0, 1, 2]), memory_mode=mode
            )
            self.assertEqual(out["prediction"].shape, (3, 1))
            self.assertEqual(out["block_weights"].shape, (3, 4))
            self.assertTrue(torch.isfinite(out["prediction"]).all())
            torch.testing.assert_close(out["block_weights"].sum(-1), torch.ones(3))

    def test_pretraining_loss_backpropagates(self) -> None:
        mask = sample_ssl_mask(self.genotype, 0.4, torch.Generator().manual_seed(9))
        out = self.model.forward_pretrain(self.genotype, self.allele_frequency, mask)
        loss = masked_genotype_loss(out["logits"], self.genotype, out["mask"])
        loss.backward()
        self.assertTrue(torch.isfinite(loss))
        self.assertIsNotNone(self.model.local_encoder.genotype_embedding.weight.grad)

    def test_global_mixer_can_be_ablated(self) -> None:
        out = self.model(
            self.genotype,
            self.allele_frequency,
            torch.tensor([0, 1, 2]),
            use_global_mixer=False,
        )
        self.assertEqual(out["global_slots"].shape, (3, 0, self.cfg.memory_dim))

    def test_grassmann_features_are_basis_invariant(self) -> None:
        u, _ = torch.linalg.qr(torch.randn(2, 3, 10, 4))
        q, _ = torch.linalg.qr(torch.randn(2, 3, 4, 4))
        probes = torch.randn(10, 7)
        original = grassmann_probe_features(u, probes)
        rotated = grassmann_probe_features(torch.einsum("bndr,bnrs->bnds", u, q), probes)
        torch.testing.assert_close(original, rotated, atol=2e-5, rtol=2e-5)

    def test_transferable_local_encoder_has_no_block_lookup(self) -> None:
        cfg = ModelConfig(
            max_blocks=4,
            snps_per_block=8,
            d_model=16,
            n_heads=4,
            local_layers=1,
            local_ff_dim=32,
            spectral_rank=3,
            dropout=0.0,
            local_position_mode="relative_continuous",
            use_local_block_embedding=False,
        )
        encoder = LocalGenotypeEncoder(cfg).eval()
        tokens = encode_genotypes(self.genotype)
        base = torch.tensor(
            [
                [100, 120, 155, 210, 290, 410, 570, 800],
                [1000, 1020, 1055, 1110, 1190, 1310, 1470, 1700],
                [2000, 2020, 2055, 2110, 2190, 2310, 2470, 2700],
                [3000, 3020, 3055, 3110, 3190, 3310, 3470, 3700],
            ],
            dtype=torch.int64,
        )
        hidden, logits = encoder(tokens, self.allele_frequency, genomic_positions=base)
        shifted_hidden, shifted_logits = encoder(
            tokens, self.allele_frequency, genomic_positions=base + 1_000_000_000
        )
        self.assertEqual(hidden.shape, (3, 4, 8, 16))
        self.assertEqual(logits.shape, (3, 4, 8, 3))
        torch.testing.assert_close(hidden, shifted_hidden, atol=1e-6, rtol=1e-6)
        torch.testing.assert_close(logits, shifted_logits, atol=1e-6, rtol=1e-6)
        with self.assertRaisesRegex(ValueError, "block_ids are prohibited"):
            encoder(tokens, self.allele_frequency, torch.arange(4), genomic_positions=base)

    def test_transferable_local_encoder_requires_coordinates(self) -> None:
        cfg = ModelConfig(
            max_blocks=1,
            snps_per_block=8,
            d_model=16,
            n_heads=4,
            local_layers=1,
            local_ff_dim=32,
            spectral_rank=3,
            local_position_mode="relative_continuous",
            use_local_block_embedding=False,
        )
        encoder = LocalGenotypeEncoder(cfg)
        with self.assertRaisesRegex(ValueError, "requires genomic_positions"):
            encoder(encode_genotypes(self.genotype[:, :1]), self.allele_frequency[:1])
        with self.assertRaisesRegex(ValueError, "integer dtype"):
            encoder(
                encode_genotypes(self.genotype[:, :1]),
                self.allele_frequency[:1],
                genomic_positions=torch.arange(8, dtype=torch.float32).reshape(1, 8),
            )


if __name__ == "__main__":
    unittest.main()
