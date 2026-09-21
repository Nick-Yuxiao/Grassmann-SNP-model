#!/usr/bin/env python3
"""B2: extract only the named columns from a SAS file into a TSV.

Nothing but the requested columns leaves the SAS file. Trait columns extracted
here belong to Task B only; the E0 genotype protocol must not read them.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_split import sha256_file  # noqa: E402


def format_value(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        if math.isnan(value):
            return ""
        if value.is_integer() and abs(value) < 1e15:
            return str(int(value))
        return repr(value)
    text = str(value).strip()
    return "" if text.lower() in {"nan", "none"} else text


def iter_chunks(
    path: Path,
    columns: list[str],
    chunk_size: int,
    row_limit: int | None,
    read_mode: str = "single",
):
    """Yield (rows, column_names) batches from a SAS file.

    read_mode "single" issues one read_sas7bdat call. That is the default because
    pyreadstat's chunked reader re-parses the file from the start for every chunk,
    including a metadata pass and an out-of-range probe past the last chunk. On a
    52 GB UKB basket that turns a ~15 minute job into hours while producing no
    output at all until the end. With usecols the whole result is only a handful
    of columns, so a single read stays small in memory.
    """
    try:
        import pyreadstat  # type: ignore
    except ImportError:
        pyreadstat = None  # type: ignore

    if pyreadstat is not None:
        if read_mode == "single":
            kwargs: dict[str, object] = {"usecols": columns}
            if row_limit is not None:
                kwargs["row_limit"] = row_limit
            frame, _ = pyreadstat.read_sas7bdat(str(path), **kwargs)
            yield (
                [list(record) for record in frame.itertuples(index=False, name=None)],
                list(frame.columns),
            )
            return
        reader = pyreadstat.read_file_in_chunks(
            pyreadstat.read_sas7bdat, str(path), chunksize=chunk_size, usecols=columns
        )
        seen = 0
        for frame, _ in reader:
            yield [list(record) for record in frame.itertuples(index=False, name=None)], list(frame.columns)
            seen += len(frame)
            if row_limit is not None and seen >= row_limit:
                return
        return

    try:
        import pandas as pd  # type: ignore
    except ImportError as error:
        raise SystemExit(
            "Neither pyreadstat nor pandas is importable. Install one of them "
            "(conda install -c conda-forge pyreadstat) before running B2."
        ) from error

    seen = 0
    for frame in pd.read_sas(str(path), format="sas7bdat", chunksize=chunk_size, iterator=True):
        frame.columns = [
            name.decode() if isinstance(name, bytes) else str(name) for name in frame.columns
        ]
        missing = [name for name in columns if name not in frame.columns]
        if missing:
            raise SystemExit(f"Columns not present in {path.name}: {missing}")
        subset = frame[columns]
        yield [list(record) for record in subset.itertuples(index=False, name=None)], columns
        seen += len(subset)
        if row_limit is not None and seen >= row_limit:
            return


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sas", type=Path, required=True)
    parser.add_argument("--column", action="append", required=True,
                        help="Exact SAS column name; repeatable. Include the identifier column.")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--read-mode", default="single", choices=["single", "chunked"],
                        help="single issues one read call; chunked re-parses the file per "
                             "chunk and is only useful when memory forces it.")
    parser.add_argument("--chunk-size", type=int, default=100_000,
                        help="Rows per chunk in chunked mode; ignored in single mode.")
    parser.add_argument("--row-limit", type=int, default=None,
                        help="Stop after this many rows; use for a smoke run.")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.out.exists() and not args.overwrite:
        raise SystemExit(f"Refusing to overwrite: {args.out}")
    args.out.parent.mkdir(parents=True, exist_ok=True)

    columns = list(dict.fromkeys(args.column))
    written = 0
    header_written = False
    with args.out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for rows, frame_columns in iter_chunks(
            args.sas, columns, args.chunk_size, args.row_limit, args.read_mode
        ):
            if not header_written:
                writer.writerow(frame_columns)
                header_written = True
            for row in rows:
                if args.row_limit is not None and written >= args.row_limit:
                    break
                writer.writerow([format_value(value) for value in row])
                written += 1
            if args.row_limit is not None and written >= args.row_limit:
                break

    summary = {
        "classification": "TASK_B_COLUMN_EXTRACT",
        "identifiers_included": True,
        "source": str(args.sas),
        "columns": columns,
        "read_mode": args.read_mode,
        "rows_written": written,
        "output": str(args.out),
        "output_sha256": sha256_file(args.out),
        "note": "Trait columns are Task B only; E0 must not read this file.",
    }
    args.out.with_suffix(args.out.suffix + ".summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({k: v for k, v in summary.items() if k != "columns"},
                     ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
