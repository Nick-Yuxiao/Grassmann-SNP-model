from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import torch


HERE = Path(__file__).resolve().parent
PACKAGE_ROOT = HERE.parents[2]
sys.path.insert(0, str(PACKAGE_ROOT / "src"))

from snp_spectral_foundation.config import ModelConfig  # noqa: E402
from snp_spectral_foundation.model import LocalGenotypeEncoder, encode_genotypes  # noqa: E402
from snp_spectral_foundation.training import masked_genotype_loss, sample_ssl_mask  # noqa: E402


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def dosage(gt: tuple[int | None, ...] | None) -> int | None:
    if gt is None or len(gt) != 2 or any(value not in (0, 1) for value in gt):
        return None
    return int(gt[0]) + int(gt[1])  # type: ignore[arg-type]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Development-only real-VCF E1 gradient smoke test.")
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--family-manifest", type=Path, required=True)
    parser.add_argument("--block-manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--samples", type=int, default=4)
    parser.add_argument("--updates", type=int, default=2)
    parser.add_argument("--seed", type=int, default=20260916)
    return parser.parse_args()


def run(args: argparse.Namespace) -> dict[str, object]:
    if args.output.exists():
        raise FileExistsError(f"refusing existing output: {args.output}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    torch.set_num_threads(4)
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    family_rows = read_tsv(args.family_manifest)
    train_ids = [row["sample_id"] for row in family_rows if row["split"] == "train"]
    if len(train_ids) < args.samples:
        raise ValueError("not enough train samples")
    smoke_ids = train_ids[: args.samples]
    train_blocks = [row for row in read_tsv(args.block_manifest) if row["split"] == "train"]
    if not train_blocks:
        raise ValueError("no train block")

    try:
        import pysam
    except ImportError as error:
        raise RuntimeError("pysam==0.24.0 is required") from error

    selected_block: dict[str, str] | None = None
    selected_genotypes: list[list[int]] = []
    selected_positions: list[int] = []
    selected_af: list[float] = []
    with pysam.VariantFile(str(args.vcf)) as variants:
        available = set(variants.header.samples)
        if any(sample_id not in available for sample_id in train_ids):
            raise ValueError("train manifest contains samples absent from VCF")
        for block in train_blocks:
            genotypes: list[list[int]] = []
            positions: list[int] = []
            allele_frequencies: list[float] = []
            for record in variants.fetch(
                block["chrom"], int(block["core_start0"]), int(block["core_end0"])
            ):
                if len(record.alleles or ()) != 2 or len(record.ref) != 1 or len(record.alts[0]) != 1:
                    continue
                smoke_dosages = [dosage(record.samples[sample_id].get("GT")) for sample_id in smoke_ids]
                if any(value is None for value in smoke_dosages):
                    continue
                train_dosages = [dosage(record.samples[sample_id].get("GT")) for sample_id in train_ids]
                called = [value for value in train_dosages if value is not None]
                if not called:
                    continue
                genotypes.append([int(value) for value in smoke_dosages if value is not None])
                positions.append(int(record.pos))
                allele_frequencies.append(float(sum(called) / (2.0 * len(called))))
                if len(positions) == 256:
                    break
            if len(positions) == 256:
                selected_block = block
                selected_genotypes = genotypes
                selected_positions = positions
                selected_af = allele_frequencies
                break
    if selected_block is None:
        raise ValueError("no train block supplied 256 eligible SNPs")

    genotype = torch.tensor(np.asarray(selected_genotypes, dtype=np.int64).T[:, None, :])
    positions = torch.tensor(np.asarray(selected_positions, dtype=np.int64)[None, :])
    allele_frequency = torch.tensor(np.asarray(selected_af, dtype=np.float32)[None, :])
    cfg = ModelConfig(
        max_blocks=1,
        snps_per_block=256,
        d_model=128,
        n_heads=8,
        local_layers=4,
        local_ff_dim=512,
        spectral_rank=8,
        dropout=0.1,
        local_position_mode="relative_continuous",
        use_local_block_embedding=False,
    )
    model = LocalGenotypeEncoder(cfg)
    if any("block_embedding" in name for name in model.state_dict()):
        raise AssertionError("E1 unexpectedly contains a block embedding")
    optimizer = torch.optim.AdamW(model.parameters(), lr=3e-4, weight_decay=0.01)
    history: list[dict[str, float | int]] = []
    started = time.perf_counter()
    model.train()
    for update in range(args.updates):
        generator = torch.Generator(device="cpu").manual_seed(args.seed + update)
        mask = sample_ssl_mask(genotype, 0.20, generator)
        optimizer.zero_grad(set_to_none=True)
        _, logits = model(
            encode_genotypes(genotype, mask),
            allele_frequency,
            genomic_positions=positions,
        )
        loss = masked_genotype_loss(logits, genotype, mask)
        loss.backward()
        grad_norm = torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()
        history.append(
            {
                "update": update + 1,
                "loss": float(loss.detach()),
                "gradient_norm_before_clip": float(grad_norm),
                "masked_targets": int(mask.sum()),
            }
        )
    elapsed = time.perf_counter() - started
    if any(not np.isfinite(item["loss"]) for item in history):
        raise ValueError("non-finite smoke loss")

    result: dict[str, object] = {
        "status": "PASS",
        "classification": "DEVELOPMENT_ONLY_NON_EVIDENCE",
        "formal_use_allowed": False,
        "checkpoint_written": False,
        "data_scope": "train_families_x_train_block_only",
        "train_sample_count_for_AF": len(train_ids),
        "smoke_sample_count": len(smoke_ids),
        "block_id": selected_block["block_id"],
        "variant_count": 256,
        "position1_min": min(selected_positions),
        "position1_max": max(selected_positions),
        "updates": args.updates,
        "history": history,
        "elapsed_seconds": elapsed,
        "model": {
            "recipe_id": "E1_PRIMARY_LOCAL_TRANSFORMER",
            "parameter_count": sum(parameter.numel() for parameter in model.parameters()),
            "local_position_mode": cfg.local_position_mode,
            "use_local_block_embedding": cfg.use_local_block_embedding,
        },
        "inputs": {
            "vcf_sha256": sha256(args.vcf),
            "family_manifest_sha256": sha256(args.family_manifest),
            "block_manifest_sha256": sha256(args.block_manifest),
        },
    }
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return result


def main() -> None:
    print(json.dumps(run(parse_args()), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
