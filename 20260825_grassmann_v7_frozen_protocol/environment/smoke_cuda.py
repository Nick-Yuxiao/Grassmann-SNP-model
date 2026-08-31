from __future__ import annotations

import argparse
import json
import math
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import torch


def version_pair(value: str) -> tuple[int, int]:
    major, minor, *_ = value.split(".")
    return int(major), int(minor)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    errors: list[str] = []
    if sys.version_info < (3, 10):
        errors.append("python_below_3.10")
    if version_pair(torch.__version__.split("+")[0]) < (2, 7):
        errors.append("torch_below_2.7")
    if not torch.cuda.is_available():
        errors.append("cuda_unavailable")

    device_payload: dict[str, object] = {}
    backward_ok = False
    finite = False
    if torch.cuda.is_available():
        props = torch.cuda.get_device_properties(0)
        device_payload = {
            "index": 0,
            "name": props.name,
            "compute_capability": [props.major, props.minor],
            "total_memory_bytes": props.total_memory,
        }
        if (props.major, props.minor) < (12, 0):
            errors.append("compute_capability_below_sm120")
        x = torch.randn(1024, 1024, device="cuda", requires_grad=True)
        loss = (x.square().mean() + x[:64, :64].matmul(x[:64, :64].T).mean())
        loss.backward()
        torch.cuda.synchronize()
        backward_ok = x.grad is not None
        finite = bool(math.isfinite(float(loss.detach().cpu()))) and bool(torch.isfinite(x.grad).all().item())
        if not backward_ok:
            errors.append("backward_missing_gradient")
        if not finite:
            errors.append("nonfinite_forward_or_gradient")

    payload = {
        "schema_version": "1.0",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "torch": torch.__version__,
        "torch_cuda_runtime": torch.version.cuda,
        "cuda_available": torch.cuda.is_available(),
        "device": device_payload,
        "backward_ok": backward_ok,
        "finite": finite,
        "status": "PASS" if not errors else "FAIL",
        "errors": errors,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
