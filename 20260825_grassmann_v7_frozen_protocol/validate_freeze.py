from __future__ import annotations

import argparse
import csv
import hashlib
import json
import py_compile
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parent
TRACKER = ROOT / "V7_任务追踪表_v7.0.2.tsv"
GRID = ROOT / "GRID_SPEC.v7.0.1.yaml"
DECISIONS = ROOT / "DECISIONS.v7.0.1.tsv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-manifest", action="store_true")
    args = parser.parse_args()

    tracker = read_tsv(TRACKER)
    expected = [f"T{i:02d}" for i in range(16)]
    assert [row["task_id"] for row in tracker] == expected
    assert len({row["task_id"] for row in tracker}) == len(tracker)
    assert tracker[0]["status"] == "IN_PROGRESS_REMOTE_ACCESS_PENDING"
    for row in tracker:
        if row["start_date"] != "TBD" and row["end_date"] != "TBD":
            assert date.fromisoformat(row["start_date"]) <= date.fromisoformat(row["end_date"])

    grid = json.loads(GRID.read_text(encoding="utf-8"))
    decisions = {row["name"]: row["value"] for row in read_tsv(DECISIONS)}
    for name in ("delta_min", "delta_NI", "delta_LD", "overfit_thr", "pc_control_thr"):
        assert float(decisions[name]) == float(grid["thresholds"][name])
    assert grid["repeats"]["data_seeds"] == [1, 2, 3, 4, 5]
    assert grid["repeats"]["init_seeds"] == [1, 2]
    assert grid["stages"]["A3"]["held_out_length"] == max(grid["stages"]["A3"]["sequence_lengths"])
    assert grid["stages"]["A2"]["primary_comparison"] == ["op_wedge_norm", "op_bilinear_learned"]
    assert len(grid["stages"]["A2"]["operators"]) == 6

    for path in (ROOT / "build_tracker.py", ROOT / "validate_freeze.py", ROOT / "environment" / "smoke_cuda.py", ROOT / "environment" / "probe_server.py"):
        py_compile.compile(str(path), doraise=True)

    files = sorted(
        path for path in ROOT.rglob("*")
        if path.is_file()
        and path.name not in {"MANIFEST.sha256"}
        and "__pycache__" not in path.parts
        and path.suffix != ".pyc"
    )
    if args.write_manifest:
        lines = [f"{sha256(path)}  {path.relative_to(ROOT).as_posix()}" for path in files]
        (ROOT / "MANIFEST.sha256").write_text("\n".join(lines) + "\n", encoding="utf-8")

    payload = {
        "status": "PASS",
        "tracker_rows": len(tracker),
        "tracker_columns": len(tracker[0]),
        "grid_protocol_version": grid["protocol_version"],
        "a2_operator_count": len(grid["stages"]["A2"]["operators"]),
        "manifest_entries": len(files),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
