#!/usr/bin/env python3
"""B1: find which UKB fields live in which SAS file, without a scientific stack.

Two methods are tried in order. If pyreadstat is importable the file metadata is
read authoritatively. Otherwise column names are recovered by scanning the head
of the file for SAS column-name text, which is where sas7bdat keeps its column
metadata pages. The fallback is heuristic and is labelled as such in the output.

Only column NAMES are read. No participant row is opened.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_split import sha256_file  # noqa: E402

NAME_PATTERNS = [
    re.compile(rb"n_\d{1,6}_\d{1,3}_\d{1,3}"),
    re.compile(rb"f\.\d{1,6}\.\d{1,3}\.\d{1,3}"),
    re.compile(rb"x\d{1,6}_\d{1,3}_\d{1,3}"),
]

# Fields Task B cares about. Keys are UKB field ids.
FIELDS_OF_INTEREST = {
    "22021": "genetic kinship to other participants (relatedness axis)",
    "22027": "outlier for heterozygosity or missing rate",
    "22019": "sex chromosome aneuploidy",
    "22020": "used in genetic principal components",
    "22006": "genetic ethnic grouping",
    "22001": "genetic sex",
    "22000": "genotype measurement batch",
    "22009": "genetic principal components",
    "31": "self-reported sex",
    "21003": "age at assessment",
    "54": "assessment centre",
    "30010": "red blood cell (erythrocyte) count",
    "30020": "haemoglobin concentration",
    "30040": "mean corpuscular volume",
    "30050": "mean corpuscular haemoglobin",
    "30070": "red blood cell distribution width",
    "30080": "platelet count",
    "30000": "white blood cell (leukocyte) count",
    "30690": "cholesterol",
    "30760": "HDL cholesterol",
    "30780": "LDL direct",
    "30870": "triglycerides",
    "30750": "glycated haemoglobin HbA1c",
}

# Continuous, well-measured, high-heritability traits preferred for qualification.
PREFERRED_TRAITS = ["30020", "30010", "30040", "30780", "30760", "30870", "30690"]


def field_id(column: str) -> str | None:
    match = re.match(r"^(?:n_|x)(\d+)_\d+_\d+$", column)
    if match:
        return match.group(1)
    match = re.match(r"^f\.(\d+)\.\d+\.\d+$", column)
    return match.group(1) if match else None


def columns_via_pyreadstat(path: Path) -> tuple[list[str], int | None] | None:
    try:
        import pyreadstat  # type: ignore
    except Exception:  # noqa: BLE001 - optional dependency
        return None
    try:
        _, meta = pyreadstat.read_sas7bdat(str(path), metadataonly=True)
    except Exception as error:  # noqa: BLE001 - corrupt or unsupported file
        raise RuntimeError(f"pyreadstat failed on {path.name}: {error}") from error
    return list(meta.column_names), meta.number_rows


def columns_via_byte_scan(path: Path, scan_bytes: int) -> list[str]:
    found: set[str] = set()
    remaining = scan_bytes
    with path.open("rb") as handle:
        previous = b""
        while remaining > 0:
            chunk = handle.read(min(1 << 22, remaining))
            if not chunk:
                break
            remaining -= len(chunk)
            window = previous + chunk
            for pattern in NAME_PATTERNS:
                for match in pattern.finditer(window):
                    found.add(match.group().decode("ascii"))
            previous = window[-64:]
    return sorted(found)


def summarize(columns: list[str]) -> dict[str, list[str]]:
    by_field: dict[str, list[str]] = {}
    for column in columns:
        identifier = field_id(column)
        if identifier is None:
            continue
        by_field.setdefault(identifier, []).append(column)
    for value in by_field.values():
        value.sort()
    return by_field


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sas-dir", type=Path, help="Directory of .sas7bdat files to scan.")
    parser.add_argument("--sas", type=Path, action="append", default=[],
                        help="Individual .sas7bdat file; repeatable.")
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--scan-bytes", type=int, default=300 << 20,
                        help="Bytes of file head scanned when pyreadstat is unavailable.")
    parser.add_argument("--hash-files", action="store_true",
                        help="Also sha256 each scanned file; slow on very large inputs.")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    candidates = list(args.sas)
    if args.sas_dir is not None:
        candidates.extend(sorted(args.sas_dir.glob("*.sas7bdat")))
    if not candidates:
        raise SystemExit("No SAS files given; pass --sas-dir or --sas")

    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    report_path = out_dir / "FIELD_SCAN.json"
    if report_path.exists() and not args.overwrite:
        raise SystemExit(f"Refusing to overwrite: {report_path}")

    files: list[dict[str, object]] = []
    for path in candidates:
        path = path.resolve()
        entry: dict[str, object] = {
            "path": str(path),
            "name": path.name,
            "size_bytes": path.stat().st_size,
        }
        if args.hash_files:
            entry["sha256"] = sha256_file(path)
        try:
            authoritative = columns_via_pyreadstat(path)
        except RuntimeError as error:
            entry["error"] = str(error)
            authoritative = None
        if authoritative is not None:
            columns, rows = authoritative
            entry["method"] = "pyreadstat_metadata"
            entry["row_count"] = rows
        else:
            columns = columns_via_byte_scan(path, args.scan_bytes)
            entry["method"] = "byte_scan_heuristic"
            entry["row_count"] = None
            entry["scan_bytes"] = args.scan_bytes
        by_field = summarize(columns)
        entry["column_count_seen"] = len(columns)
        entry["distinct_field_ids"] = len(by_field)
        entry["fields_of_interest"] = {
            identifier: {
                "description": description,
                "present": identifier in by_field,
                "columns": by_field.get(identifier, [])[:6],
                "column_count": len(by_field.get(identifier, [])),
            }
            for identifier, description in FIELDS_OF_INTEREST.items()
        }
        files.append(entry)

    relatedness = [f for f in files if f["fields_of_interest"]["22021"]["present"]]  # type: ignore[index]
    trait_ready = []
    for identifier in PREFERRED_TRAITS:
        holders = [f["name"] for f in files if f["fields_of_interest"][identifier]["present"]]  # type: ignore[index]
        if holders:
            trait_ready.append({
                "field": identifier,
                "description": FIELDS_OF_INTEREST[identifier],
                "files": holders,
            })

    decision: dict[str, object] = {
        "relatedness_field_22021_found": bool(relatedness),
        "relatedness_source_files": [f["name"] for f in relatedness],
        "qualifying_trait_candidates": trait_ready,
        "next_step": (
            "freeze the unrelated subset with b3_build_unrelated_split.py"
            if relatedness
            else "22021 absent: obtain the official kinship file, or derive components with "
                 "LD-pruned KING (see KING_FALLBACK.zh-CN.md). Do not start the bridge on a "
                 "random individual split."
        ),
    }

    report = {
        "classification": "TASK_B_FIELD_SCAN",
        "identifiers_included": False,
        "phenotype_values_read": False,
        "files": files,
        "decision": decision,
    }
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    lines = ["# Task B 字段扫描", ""]
    for entry in files:
        marks = []
        for identifier in ["22021", "22006", "22027", *PREFERRED_TRAITS]:
            info = entry["fields_of_interest"][identifier]  # type: ignore[index]
            marks.append(f"{identifier}{'✅' if info['present'] else '❌'}")
        lines.append(
            f"- `{entry['name']}` · {int(entry['size_bytes'])/(1<<30):.1f} GiB · "
            f"{entry['method']} · 字段 {entry['distinct_field_ids']} 个 · " + " ".join(marks)
        )
    lines += ["", f"**22021 是否存在：{'是' if relatedness else '否'}**", "", decision["next_step"], ""]
    (out_dir / "FIELD_SCAN.md").write_text("\n".join(lines), encoding="utf-8")

    print(json.dumps({"status": "OK", "report": str(report_path), **decision},
                     ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
