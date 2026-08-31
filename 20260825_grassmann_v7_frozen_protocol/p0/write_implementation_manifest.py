from __future__ import annotations

import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
P0 = ROOT / "p0"
SERVER_OPS = ROOT / "server_ops"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    files = [ROOT / "METRIC_DEFINITIONS.md"]
    files.extend(
        path for path in P0.rglob("*")
        if path.is_file()
        and path.name != "P0_IMPLEMENTATION.sha256"
        and "__pycache__" not in path.parts
        and "runtime_verdict" not in path.parts
        and path.suffix != ".pyc"
    )
    files.extend(
        path for path in SERVER_OPS.rglob("*")
        if path.is_file()
        and "__pycache__" not in path.parts
        and path.suffix != ".pyc"
    )
    lines = [f"{sha256(path)}  {path.relative_to(ROOT).as_posix()}" for path in sorted(files)]
    output = P0 / "P0_IMPLEMENTATION.sha256"
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {len(lines)} entries to {output}")


if __name__ == "__main__":
    main()
