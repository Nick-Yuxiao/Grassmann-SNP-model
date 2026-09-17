#!/usr/bin/env python3
"""B2b: audit the observed coding of a flag column before freezing any filter.

A field being present is not the same as a filter being correct. This prints the
full observed value histogram for UKB field 22021 and any QC flags, separates
"unrelated" from "relatedness unknown", and reports how far the flag file's
participants overlap the covariate/genotype file, which is the other way a
filter silently goes wrong.

Nothing is filtered here and no split is written. Standard library only.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_split import read_tsv, sha256_file  # noqa: E402

# What each documented value of field 22021 means for a relatedness-disjoint design.
FIELD_22021_SEMANTICS = {
    "0": {"meaning": "no kinship found", "disposition": "keep", "relatedness": "known_unrelated"},
    "1": {"meaning": "at least one relative identified", "disposition": "drop",
          "relatedness": "known_related"},
    "10": {"meaning": "ten or more third-degree relatives identified", "disposition": "drop",
           "relatedness": "known_related"},
    "-1": {"meaning": "excluded from the kinship inference process", "disposition": "drop",
           "relatedness": "UNKNOWN"},
    "<blank>": {"meaning": "no value recorded", "disposition": "drop", "relatedness": "UNKNOWN"},
}


def normalize(value: str | None) -> str:
    text = (value or "").strip()
    if not text:
        return "<blank>"
    try:
        return str(int(float(text)))
    except ValueError:
        return text


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", type=Path, required=True,
                        help="TSV produced by b2_extract_sas_columns.py.")
    parser.add_argument("--id-column", default="n_eid")
    parser.add_argument("--column", action="append", required=True,
                        help="Flag column to audit; repeatable. Put 22021 first.")
    parser.add_argument("--cross-reference", type=Path, default=None,
                        help="Covariates TSV, to measure participant overlap.")
    parser.add_argument("--cross-id-column", default="eid")
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    report_path = out_dir / "FLAG_AUDIT.json"
    if report_path.exists() and not args.overwrite:
        raise SystemExit(f"Refusing to overwrite: {report_path}")

    histograms: dict[str, Counter[str]] = {name: Counter() for name in args.column}
    identifiers: set[str] = set()
    duplicate_rows = 0
    rows = 0
    for row in read_tsv(args.table):
        rows += 1
        identifier = (row.get(args.id_column) or "").strip()
        if identifier:
            if identifier in identifiers:
                duplicate_rows += 1
            identifiers.add(identifier)
        for name in args.column:
            if name not in row:
                raise SystemExit(f"Column not present in {args.table.name}: {name}")
            histograms[name][normalize(row.get(name))] += 1

    columns_report: dict[str, object] = {}
    for name, histogram in histograms.items():
        total = sum(histogram.values())
        documented = FIELD_22021_SEMANTICS if "22021" in name else {}
        values = []
        for value, count in histogram.most_common():
            entry = {
                "value": value,
                "count": count,
                "share": round(count / total, 6) if total else 0.0,
            }
            if documented:
                known = documented.get(value)
                entry["documented"] = bool(known)
                entry["meaning"] = known["meaning"] if known else "NOT IN THE DOCUMENTED CODING"
                entry["disposition"] = known["disposition"] if known else "drop (undocumented)"
                entry["relatedness"] = known["relatedness"] if known else "UNKNOWN"
            values.append(entry)
        report_column: dict[str, object] = {"rows": total, "distinct_values": len(histogram),
                                            "values": values}
        if documented:
            keep = histogram.get("0", 0)
            known_related = histogram.get("1", 0) + histogram.get("10", 0)
            unknown = sum(
                count for value, count in histogram.items()
                if documented.get(value, {}).get("relatedness", "UNKNOWN") == "UNKNOWN"
            )
            undocumented = sorted(value for value in histogram if value not in documented)
            report_column["disposition_summary"] = {
                "keep_known_unrelated": keep,
                "drop_known_related": known_related,
                "drop_relatedness_unknown": unknown,
                "undocumented_values_present": undocumented,
                "keep_share": round(keep / total, 6) if total else 0.0,
            }
        columns_report[name] = report_column

    overlap: dict[str, object] = {}
    if args.cross_reference is not None:
        other: set[str] = set()
        for row in read_tsv(args.cross_reference):
            value = (row.get(args.cross_id_column) or "").strip()
            if value:
                other.add(value)
        overlap = {
            "cross_reference": str(args.cross_reference),
            "flag_file_participants": len(identifiers),
            "cross_reference_participants": len(other),
            "in_both": len(identifiers & other),
            "flag_only": len(identifiers - other),
            "cross_reference_only": len(other - identifiers),
            "coverage_of_cross_reference": (
                round(len(identifiers & other) / len(other), 6) if other else None
            ),
        }

    report = {
        "classification": "TASK_B_FLAG_CODING_AUDIT",
        "identifiers_included": False,
        "table": str(args.table),
        "table_sha256": sha256_file(args.table),
        "rows_read": rows,
        "duplicate_identifier_rows": duplicate_rows,
        "columns": columns_report,
        "participant_overlap": overlap,
        "filter_rule": (
            "keep only field 22021 == 0. Values 1 and 10 are known related; -1, blank, "
            "absent and any undocumented value mean relatedness is UNKNOWN and are dropped "
            "rather than assumed unrelated, because an unknown-relatedness participant is "
            "exactly the leakage the design is built to exclude."
        ),
        "verify_against": "UKB Showcase field 22021 and its data-coding page",
    }
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    lines = ["# Flag coding audit", "", f"- 行数 {rows}；重复标识行 {duplicate_rows}", ""]
    for name, column in columns_report.items():
        lines.append(f"## `{name}`")
        lines.append("")
        lines.append("| value | count | share | 含义 | 处置 |")
        lines.append("| --- | --- | --- | --- | --- |")
        for entry in column["values"]:  # type: ignore[index]
            lines.append(
                f"| `{entry['value']}` | {entry['count']} | {entry['share']:.4f} | "
                f"{entry.get('meaning', '-')} | {entry.get('disposition', '-')} |"
            )
        summary = column.get("disposition_summary")  # type: ignore[union-attr]
        if summary:
            lines += [
                "",
                f"- 保留（已知无亲缘）：**{summary['keep_known_unrelated']}**"
                f"（{summary['keep_share']:.1%}）",
                f"- 剔除（已知有亲缘）：{summary['drop_known_related']}",
                f"- 剔除（亲缘未知）：{summary['drop_relatedness_unknown']}",
                f"- 未记载的取值：{summary['undocumented_values_present'] or '无'}",
            ]
        lines.append("")
    if overlap:
        lines += [
            "## 参与者覆盖",
            "",
            f"- flag 文件 {overlap['flag_file_participants']}；对照文件 "
            f"{overlap['cross_reference_participants']}",
            f"- 交集 {overlap['in_both']}；仅 flag {overlap['flag_only']}；"
            f"仅对照 {overlap['cross_reference_only']}",
            f"- 对照文件被覆盖比例 {overlap['coverage_of_cross_reference']}",
            "",
        ]
    lines += ["> 取值含义须与 UKB Showcase field 22021 的 data-coding 页面核对后才算冻结。", ""]
    (out_dir / "FLAG_AUDIT.md").write_text("\n".join(lines), encoding="utf-8")

    printable = {
        "status": "OK",
        "report": str(report_path),
        "rows_read": rows,
    }
    for name, column in columns_report.items():
        if "disposition_summary" in column:  # type: ignore[operator]
            printable[name] = column["disposition_summary"]  # type: ignore[index]
    if overlap:
        printable["participant_overlap"] = overlap
    print(json.dumps(printable, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
