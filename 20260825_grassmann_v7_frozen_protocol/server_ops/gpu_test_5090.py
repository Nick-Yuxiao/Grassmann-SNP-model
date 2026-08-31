from __future__ import annotations

import argparse
import json
import math
import os
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import torch


def pair(version: str | None) -> tuple[int, int] | None:
    if not version:
        return None
    fields = version.split("+")[0].split(".")
    try:
        return int(fields[0]), int(fields[1])
    except (IndexError, ValueError):
        return None


def timed_matmul(size: int, iterations: int, dtype: torch.dtype) -> dict[str, object]:
    left = torch.randn((size, size), device="cuda", dtype=dtype, requires_grad=True)
    right = torch.randn((size, size), device="cuda", dtype=dtype, requires_grad=True)
    torch.cuda.synchronize()
    started = time.perf_counter()
    last_loss = None
    for _ in range(iterations):
        product = left @ right
        loss = product.float().square().mean()
        loss.backward()
        if left.grad is None or right.grad is None:
            raise RuntimeError("missing_gradient")
        left.grad = None
        right.grad = None
        last_loss = loss.detach()
    torch.cuda.synchronize()
    seconds = time.perf_counter() - started
    finite = bool(last_loss is not None and torch.isfinite(last_loss).item())
    return {
        "dtype": str(dtype).replace("torch.", ""),
        "matrix_size": size,
        "iterations": iterations,
        "seconds": seconds,
        "iterations_per_second": iterations / seconds,
        "finite": finite,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Short RTX 5090 allocation/forward/backward test.")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--matrix-size", type=int, default=2048)
    parser.add_argument("--iterations", type=int, default=10)
    parser.add_argument("--physical-gpu-index", required=True)
    args = parser.parse_args()
    errors: list[str] = []
    tests: list[dict[str, object]] = []
    device: dict[str, object] = {}

    if sys.version_info < (3, 10):
        errors.append("python_below_3.10")
    torch_pair = pair(torch.__version__)
    if torch_pair is None or torch_pair < (2, 7):
        errors.append("torch_below_2.7")
    cuda_pair = pair(torch.version.cuda)
    if cuda_pair is None or cuda_pair < (12, 8):
        errors.append("torch_cuda_runtime_below_12.8")
    if not torch.cuda.is_available():
        errors.append("cuda_unavailable")

    if torch.cuda.is_available():
        props = torch.cuda.get_device_properties(0)
        device = {
            "logical_index": 0,
            "physical_index_selected_by_wrapper": args.physical_gpu_index,
            "name": props.name,
            "compute_capability": [props.major, props.minor],
            "total_memory_bytes": props.total_memory,
            "architectures_in_torch_build": torch.cuda.get_arch_list(),
        }
        if (props.major, props.minor) < (12, 0):
            errors.append("compute_capability_below_sm120")
        try:
            torch.cuda.reset_peak_memory_stats()
            tests.append(timed_matmul(args.matrix_size, args.iterations, torch.float32))
            if torch.cuda.is_bf16_supported():
                tests.append(timed_matmul(args.matrix_size, args.iterations, torch.bfloat16))
            else:
                errors.append("bf16_not_supported")
            if not all(bool(test["finite"]) for test in tests):
                errors.append("nonfinite_forward_or_gradient")
        except torch.OutOfMemoryError:
            errors.append("unexpected_cuda_oom")
        except RuntimeError as exc:
            errors.append(f"runtime_error:{type(exc).__name__}")

    payload = {
        "schema_version": "1.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PASS" if not errors else "FAIL",
        "python": sys.version,
        "platform": platform.platform(),
        "torch": torch.__version__,
        "torch_cuda_runtime": torch.version.cuda,
        "cuda_visible_devices_was_set": "CUDA_VISIBLE_DEVICES" in os.environ,
        "device": device,
        "tests": tests,
        "peak_allocated_bytes": torch.cuda.max_memory_allocated() if torch.cuda.is_available() else None,
        "peak_reserved_bytes": torch.cuda.max_memory_reserved() if torch.cuda.is_available() else None,
        "errors": errors,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
