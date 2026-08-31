from __future__ import annotations

"""V7 CPU learning-smoke: proves the frozen Grassmann architecture can train and
learn on a *structured, learnable* synthetic task.

Classification: ENGINEERING_NON_EVIDENCE.

This binds to the FROZEN model in ``profile_models_v7_1.py`` (it does not
re-implement any architecture). Unlike the T03 profiler -- whose ``make_batch``
draws i.i.d. uniform genotypes, so masked cross-entropy stays pinned at
``ln(3) ~= 1.0986`` and never falls -- this smoke feeds a task with a real,
recoverable signal so a working model must drive the loss below that baseline.

It answers exactly one question: "does the model run forward/backward and reduce
loss on learnable structure?" It does NOT compare architectures, touch real 1KGP
data, read any holdout, or authorize any GPU / A1 / evidence-chain work.

Learnable DGP (per sample b, site j):

    target[b, j] = ( pos_pattern[j mod period] + shift[b] ) mod 3

- ``pos_pattern`` is a fixed random vector over one attention window
  (period = attention_window = 256), so the model's periodic position embedding
  can represent it.
- ``shift[b]`` in {0,1,2} is a per-sample global offset. It is exposed cleanly to
  the PC arm through ``pcs[:, 0]`` and is also inferable from any unmasked token
  by the local and grassmann arms. So all three arms have a learnable path.

Two deliberate choices, both documented in SMOKE_PROTOCOL.md:
- length defaults to 512 (>= 2 blocks of 256) so the GrassmannBlockMixer's global
  wedge channel is non-degenerate; at L <= 256 it is identically zero.
- mask_rate defaults to 0.15 (not the protocol's 0.90) to leave ample context for
  a clean, fast "it learns" signal. This is a smoke, not the protocol horizon.
"""

import argparse
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path


def _import_frozen_model(model_dir: str | None):
    """Import the frozen architecture module without copying it.

    Resolution order: explicit --model-dir / $V7_MODEL_DIR, then this script's
    own directory (works when the smoke is dropped next to profile_models_v7_1.py
    inside a release p0/).
    """
    candidates: list[Path] = []
    if model_dir:
        candidates.append(Path(model_dir))
    env_dir = os.environ.get("V7_MODEL_DIR")
    if env_dir:
        candidates.append(Path(env_dir))
    candidates.append(Path(__file__).resolve().parent)
    for cand in candidates:
        if (cand / "profile_models_v7_1.py").is_file():
            sys.path.insert(0, str(cand))
            import profile_models_v7_1 as frozen  # type: ignore

            return frozen
    searched = ", ".join(str(c) for c in candidates)
    raise SystemExit(
        "cannot locate profile_models_v7_1.py; pass --model-dir or set V7_MODEL_DIR "
        f"(searched: {searched})"
    )


def build_pos_pattern(period: int, seed: int, torch):
    generator = torch.Generator(device="cpu").manual_seed(seed)
    return torch.randint(0, 3, (period,), generator=generator, dtype=torch.long)


def structured_batch(
    *,
    torch,
    length: int,
    batch_size: int,
    pc_dim: int,
    seed: int,
    mask_rate: float,
    pos_pattern,
    period: int,
    device,
):
    """Return (tokens, pcs, target, mask) for the learnable DGP above.

    All tensors are built on CPU (pos_pattern is a CPU buffer) and moved to
    ``device`` only at return, so the CUDA path never mixes devices.
    """
    generator = torch.Generator(device="cpu").manual_seed(seed)
    shift = torch.randint(0, 3, (batch_size,), generator=generator, dtype=torch.long)
    positions = torch.arange(length) % period
    base = pos_pattern[positions]  # [length]
    target = (base.unsqueeze(0) + shift.unsqueeze(1)) % 3  # [batch, length]

    mask = torch.rand((batch_size, length), generator=generator) < mask_rate
    if not bool(mask.any()):
        mask[0, 0] = True

    tokens = target.clone()
    tokens[mask] = 3  # mask token == genotype_states

    # Expose the global shift to the PC arm on channel 0 (centered), plus mild noise
    # on all channels so the projection is not a trivial constant.
    pcs = 0.10 * torch.randn((batch_size, pc_dim), generator=generator)
    pcs[:, 0] = shift.float() - 1.0

    return (
        tokens.to(device),
        pcs.to(device),
        target.to(device),
        mask.to(device),
    )


def train_one_arm(
    *,
    frozen,
    torch,
    kind: str,
    cfg,
    length: int,
    steps: int,
    batch_size: int,
    lr: float,
    mask_rate: float,
    seed: int,
    pos_pattern,
    period: int,
    device,
    eval_every: int,
):
    F = torch.nn.functional
    torch.manual_seed(seed)
    model = frozen.MaskedGenotypeModel(kind, cfg).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=lr)

    # Fixed held-out eval batch (never trained on) for a clean learning curve.
    eval_tokens, eval_pcs, eval_target, eval_mask = structured_batch(
        torch=torch,
        length=length,
        batch_size=max(batch_size, 8),
        pc_dim=cfg.pc_dim,
        seed=seed + 777_777,
        mask_rate=mask_rate,
        pos_pattern=pos_pattern,
        period=period,
        device=device,
    )

    def eval_loss() -> float:
        model.eval()
        with torch.no_grad():
            logits = model(eval_tokens, eval_pcs)
            loss = F.cross_entropy(logits[eval_mask], eval_target[eval_mask])
        model.train()
        return float(loss.detach().cpu())

    curve: list[dict[str, float]] = []
    model.train()
    initial_eval = eval_loss()
    curve.append({"step": 0, "eval_loss": initial_eval})

    train_losses: list[float] = []
    for step in range(1, steps + 1):
        tokens, pcs, target, mask = structured_batch(
            torch=torch,
            length=length,
            batch_size=batch_size,
            pc_dim=cfg.pc_dim,
            seed=seed + step,
            mask_rate=mask_rate,
            pos_pattern=pos_pattern,
            period=period,
            device=device,
        )
        optimizer.zero_grad(set_to_none=True)
        logits = model(tokens, pcs)
        loss = F.cross_entropy(logits[mask], target[mask])
        if not bool(torch.isfinite(loss).item()):
            raise RuntimeError(f"{kind}: non-finite training loss at step {step}")
        loss.backward()
        optimizer.step()
        train_losses.append(float(loss.detach().cpu()))
        if step % eval_every == 0 or step == steps:
            curve.append({"step": step, "eval_loss": eval_loss()})

    final_eval = curve[-1]["eval_loss"]
    min_eval = min(point["eval_loss"] for point in curve)
    return {
        "model": kind,
        "parameter_count": frozen.parameter_count(model),
        "initial_eval_loss": initial_eval,
        "final_eval_loss": final_eval,
        "min_eval_loss": min_eval,
        "final_train_loss": train_losses[-1] if train_losses else None,
        "eval_curve": curve,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="V7 CPU learning-smoke (non-evidence)")
    parser.add_argument("--model-dir", default=None, help="dir containing profile_models_v7_1.py")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--length", type=int, default=512)
    parser.add_argument("--steps", type=int, default=300)
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--mask-rate", type=float, default=0.15)
    parser.add_argument("--seed", type=int, default=70101)
    parser.add_argument("--eval-every", type=int, default=25)
    parser.add_argument("--device", default="cpu", choices=["cpu", "cuda"])
    parser.add_argument("--success-threshold", type=float, default=0.30)
    parser.add_argument("--models", default=None, help="comma list; default all three frozen kinds")
    args = parser.parse_args()

    frozen = _import_frozen_model(args.model_dir)
    import torch  # noqa: E402  (imported after path resolution / frozen import)

    if args.length < 512:
        print(
            "WARNING: length < 512 makes the GrassmannBlockMixer global channel "
            "degenerate (single 256-block => zero wedge).",
            file=sys.stderr,
        )
    if args.device == "cuda" and not torch.cuda.is_available():
        raise SystemExit("--device cuda requested but CUDA is not available")
    device = torch.device(args.device)

    cfg = frozen.Architecture()
    cfg.validate()
    period = cfg.attention_window
    pos_pattern = build_pos_pattern(period, args.seed, torch)  # CPU buffer; batches move to device

    kinds = (
        tuple(x.strip() for x in args.models.split(",") if x.strip())
        if args.models
        else frozen.MODEL_KINDS
    )
    unknown = sorted(set(kinds) - set(frozen.MODEL_KINDS))
    if unknown:
        raise SystemExit(f"unknown models: {unknown}")

    baseline = float(__import__("math").log(cfg.genotype_states))
    results = []
    for kind in kinds:
        row = train_one_arm(
            frozen=frozen,
            torch=torch,
            kind=kind,
            cfg=cfg,
            length=args.length,
            steps=args.steps,
            batch_size=args.batch_size,
            lr=args.lr,
            mask_rate=args.mask_rate,
            seed=args.seed,
            pos_pattern=pos_pattern,
            period=period,
            device=device,
            eval_every=args.eval_every,
        )
        row["learned"] = bool(row["final_eval_loss"] < args.success_threshold)
        row["below_random_baseline"] = bool(row["final_eval_loss"] < baseline - 0.05)
        results.append(row)
        print(
            f"{kind}: init={row['initial_eval_loss']:.4f} "
            f"final={row['final_eval_loss']:.4f} min={row['min_eval_loss']:.4f} "
            f"learned={row['learned']}"
        )

    all_learned = all(r["learned"] for r in results)
    status = "LEARNING_SMOKE_PASS" if all_learned else "LEARNING_SMOKE_INCOMPLETE"
    payload = {
        "schema_version": "1.0",
        "classification": "ENGINEERING_NON_EVIDENCE",
        "purpose": "prove frozen architecture trains and reduces loss on learnable synthetic structure",
        "does_not_authorize": [
            "architecture_comparison_or_ranking",
            "a1_efficiency_or_capacity_claim",
            "gpu_or_server_evidence_run",
            "real_data_or_holdout_read",
        ],
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "device_type": device.type,
        "torch_version": torch.__version__,
        "random_baseline_cross_entropy": baseline,
        "config": {
            "length": args.length,
            "steps": args.steps,
            "batch_size": args.batch_size,
            "lr": args.lr,
            "mask_rate": args.mask_rate,
            "seed": args.seed,
            "success_threshold": args.success_threshold,
        },
        "architecture": frozen.asdict(cfg) if hasattr(frozen, "asdict") else None,
        "results": results,
        "status": status,
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    out = args.output_dir / "LEARNING_SMOKE_REPORT.json"
    out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(out), "status": status}, indent=2))
    if not all_learned:
        raise SystemExit(4)


if __name__ == "__main__":
    main()
