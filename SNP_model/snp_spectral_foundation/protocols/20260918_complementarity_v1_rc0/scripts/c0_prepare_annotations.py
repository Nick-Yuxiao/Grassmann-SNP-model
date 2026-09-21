#!/usr/bin/env python3
"""C0: fetch erythroid-relevant, phenotype-free annotation and write c2's manifest.

The panel is GRCh37/hg19, so every track here is hg19 and no lift-over happens. A
GRCh38 track would join to nothing and c2 would stop; that is the intended failure,
not something to work around.

Phenotype-free is a hard rail, not a preference. A track built from trait
associations (GWAS catalogue, same-trait eQTL, an existing PRS weight set) would
put the outcome into a predictor and make the F arm meaningless. The catalogue
below is cell-biology only: chromatin state, accessibility, conservation, gene
structure. Any --local track whose name looks association-derived is refused
outright.

Offline machines: run with --offline after copying the files in by hand, and point
each one at --local NAME=PATH. Nothing here needs network at analysis time.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_arms import read_tsv, sha256_file  # noqa: E402

# Names that would smuggle the outcome into a predictor. Checked on --local tracks.
# Word boundaries are not enough here: "ukb_PRS_weights" has no \b before the P,
# because an underscore is a word character. Underscores are separators in practice.
ASSOCIATION_DERIVED = re.compile(
    r"gwas|catalog|catalogue|(?<![a-z])prs(?![a-z])|(?<![a-z])pgs(?![a-z])|"
    r"eqtl|sqtl|pqtl|twas|finemap|heritab|sumstat|ldsc|magma",
    re.IGNORECASE,
)

CATALOGUE: dict[str, dict[str, object]] = {
    "k562_chromhmm_broad": {
        "url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/"
               "wgEncodeBroadHmmK562HMM.txt.gz",
        "format": "ucsc_bed_bin",
        "split_by_name": True,
        "erythroid": True,
        "why_phenotype_free": (
            "Broad ChromHMM states called from histone ChIP-seq in K562, an erythroleukemia "
            "line. Built from chromatin marks alone; no trait or association enters it."
        ),
        "caveat": "K562 is a cancer line standing in for erythroid chromatin, not primary cells.",
    },
    "roadmap_e123_k562_chromhmm": {
        "url": "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/"
               "ChmmModels/coreMarks/jointModel/final/E123_15_coreMarks_dense.bed.gz",
        "format": "bed",
        "split_by_name": True,
        "erythroid": True,
        "why_phenotype_free": (
            "Roadmap 15-state core-marks segmentation for E123 (K562), from five histone "
            "marks. Chromatin only."
        ),
        "caveat": "Same cell-line caveat as the Broad states; the two overlap heavily.",
    },
    "dnase_clusters": {
        "url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/"
               "wgEncodeRegDnaseClusteredV3.txt.gz",
        "format": "ucsc_bed_bin",
        "split_by_name": False,
        "erythroid": False,
        "why_phenotype_free": "DNase accessibility clustered over ENCODE cell types.",
        "caveat": "Pan-tissue, not erythroid-specific. Keep it as the generic-regulatory "
                  "contrast against the K562 states.",
    },
    "phastcons100way": {
        "url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/"
               "phastConsElements100way.txt.gz",
        "format": "ucsc_bed_bin",
        "split_by_name": False,
        "erythroid": False,
        "why_phenotype_free": (
            "Cross-species conservation elements. Computed from genome alignments; it has "
            "never seen a human phenotype."
        ),
        "caveat": "Conservation is the null-ish functional axis: strong, generic, not erythroid.",
    },
    "tss_proximity": {
        "url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz",
        "format": "refgene_tss",
        "split_by_name": False,
        "erythroid": False,
        "why_phenotype_free": "Gene models only. Distance to the nearest transcription start.",
        "caveat": "Structural, not regulatory; included because it is nearly free and "
                  "anchors how much of F is just gene density.",
    },
}

DEFAULT_TRACKS = ["k562_chromhmm_broad", "phastcons100way", "tss_proximity"]
TSS_DECAY_BP = 10000


def download(url: str, destination: Path, timeout: int) -> None:
    if shutil.which("curl") is None:
        raise SystemExit("curl not found; fetch the files by hand and rerun with --offline")
    result = subprocess.run(
        ["curl", "-sS", "--fail", "--location", "--max-time", str(timeout), "-o",
         str(destination), url],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise SystemExit(
            f"Download failed for {url}\n{result.stderr.strip()}\n"
            "If this machine has no outbound access, fetch the file elsewhere and pass "
            "--offline --local <track>=<path>."
        )


def read_rows(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            yield line.rstrip("\n").split("\t")


def to_intervals(path: Path, layout: str) -> list[tuple[str, int, int, str, float]]:
    """Normalise a source file to (chrom, start0, end, name, score)."""
    intervals: list[tuple[str, int, int, str, float]] = []
    if layout == "refgene_tss":
        seen: set[tuple[str, int]] = set()
        for fields in read_rows(path):
            if len(fields) < 7:
                continue
            chrom, strand = fields[2].replace("chr", ""), fields[3]
            try:
                tss = int(fields[4]) if strand == "+" else int(fields[5])
            except ValueError:
                continue
            if (chrom, tss) in seen:
                continue
            seen.add((chrom, tss))
            # A window either side; the score decays with distance inside c0, not c2.
            intervals.append((chrom, max(tss - 5 * TSS_DECAY_BP, 0), tss + 5 * TSS_DECAY_BP,
                              "tss", float(tss)))
        return intervals

    offset = 1 if layout == "ucsc_bed_bin" else 0
    for fields in read_rows(path):
        if len(fields) < offset + 3:
            continue
        chrom = fields[offset].replace("chr", "")
        try:
            start, end = int(fields[offset + 1]), int(fields[offset + 2])
        except ValueError:
            continue
        name = fields[offset + 3] if len(fields) > offset + 3 else "interval"
        score = 1.0
        if len(fields) > offset + 4:
            try:
                score = float(fields[offset + 4])
            except ValueError:
                score = 1.0
        intervals.append((chrom, start, end, name, score))
    return intervals


def safe_name(text: str) -> str:
    return re.sub(r"[^0-9A-Za-z]+", "_", text).strip("_").lower() or "state"


def write_bed(path: Path, rows: list[tuple[str, int, int, str, float]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        for chrom, start, end, name, score in sorted(rows, key=lambda r: (r[0], r[1])):
            handle.write(f"{chrom}\t{start}\t{end}\t{name}\t{score:g}\n")


def write_tss_points(path: Path, intervals: list[tuple[str, int, int, str, float]],
                     variants: list[dict[str, str]]) -> int:
    """Per-variant proximity score, exp(-distance / decay), written as a points table."""
    by_chrom: dict[str, list[int]] = defaultdict(list)
    for chrom, _, _, _, tss in intervals:
        by_chrom[chrom].append(int(tss))
    for chrom in by_chrom:
        by_chrom[chrom].sort()
    from bisect import bisect_left
    import math

    annotated = 0
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("chrom\tpos\tscore\n")
        for row in variants:
            chrom, position = row["chrom"].replace("chr", ""), int(row["pos"])
            sites = by_chrom.get(chrom)
            if not sites:
                continue
            slot = bisect_left(sites, position)
            candidates = [sites[i] for i in (slot - 1, slot) if 0 <= i < len(sites)]
            if not candidates:
                continue
            distance = min(abs(position - site) for site in candidates)
            score = math.exp(-distance / TSS_DECAY_BP)
            if score > 1e-6:
                annotated += 1
            handle.write(f"{chrom}\t{position}\t{score:.6g}\n")
    return annotated


def coverage_of_bed(rows: list[tuple[str, int, int, str, float]],
                    variants: list[dict[str, str]]) -> int:
    from bisect import bisect_right

    by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for chrom, start, end, _, _ in rows:
        by_chrom[chrom].append((start, end))
    for chrom in by_chrom:
        by_chrom[chrom].sort()
    starts_by_chrom = {chrom: [span[0] for span in spans] for chrom, spans in by_chrom.items()}
    hits = 0
    for row in variants:
        chrom = row["chrom"].replace("chr", "")
        spans = by_chrom.get(chrom)
        if not spans:
            continue
        position0 = int(row["pos"]) - 1
        index = bisect_right(starts_by_chrom[chrom], position0) - 1
        # Intervals can overlap, so walk back a bounded number of them.
        for candidate in range(index, max(index - 32, -1), -1):
            if candidate < 0:
                break
            if spans[candidate][0] <= position0 < spans[candidate][1]:
                hits += 1
                break
    return hits


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--panel-dir", type=Path, required=True,
                        help="Panel directory; panel_variants.tsv drives the coverage report.")
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--track", action="append", default=[],
                        choices=sorted(CATALOGUE), help=f"Default: {' '.join(DEFAULT_TRACKS)}")
    parser.add_argument("--local", action="append", default=[],
                        metavar="NAME=PATH",
                        help="Use an already-downloaded file for a catalogue track, or add a "
                             "new one with --local-format.")
    parser.add_argument("--local-format", default="bed",
                        choices=["bed", "ucsc_bed_bin", "refgene_tss"])
    parser.add_argument("--local-split-by-name", action="store_true")
    parser.add_argument("--offline", action="store_true",
                        help="Never call the network; every track must come from --local.")
    parser.add_argument("--timeout", type=int, default=1800)
    parser.add_argument("--min-coverage", type=float, default=0.002,
                        help="Drop a track annotating a smaller share of panel variants.")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = out_dir / "annotation_manifest.tsv"
    if manifest_path.exists() and not args.overwrite:
        raise SystemExit(f"Refusing to overwrite: {manifest_path}")

    local: dict[str, Path] = {}
    for item in args.local:
        if "=" not in item:
            raise SystemExit(f"--local wants NAME=PATH, got {item}")
        name, path = item.split("=", 1)
        if ASSOCIATION_DERIVED.search(name) or ASSOCIATION_DERIVED.search(path):
            raise SystemExit(
                f"Refusing track {name!r}: the name or path looks association-derived. "
                "A GWAS, eQTL or PRS track would put the outcome into the F arm."
            )
        local[name] = Path(path)

    wanted = args.track or list(DEFAULT_TRACKS)
    for name in local:
        if name not in wanted:
            wanted.append(name)

    variants = list(read_tsv(args.panel_dir.resolve() / "panel_variants.tsv"))
    if not variants:
        raise SystemExit("panel_variants.tsv is empty")

    manifest_rows: list[dict[str, object]] = []
    track_reports: list[dict[str, object]] = []
    raw_dir = out_dir / "raw"
    raw_dir.mkdir(exist_ok=True)

    for name in wanted:
        entry = CATALOGUE.get(name)
        if entry is None:
            entry = {
                "url": None,
                "format": args.local_format,
                "split_by_name": args.local_split_by_name,
                "erythroid": False,
                "why_phenotype_free": "declared by the operator when passing --local",
                "caveat": "not from the reviewed catalogue",
            }
        if name in local:
            source = local[name]
            if not source.exists():
                raise SystemExit(f"--local file missing for {name}: {source}")
        else:
            if args.offline:
                raise SystemExit(f"--offline given but {name} has no --local file")
            suffix = ".gz" if str(entry["url"]).endswith(".gz") else ""
            source = raw_dir / f"{name}{suffix}"
            if not source.exists():
                download(str(entry["url"]), source, args.timeout)

        layout = str(entry["format"])
        intervals = to_intervals(source, layout)
        if not intervals:
            track_reports.append({"track": name, "status": "EMPTY_SOURCE",
                                  "source": str(source)})
            continue

        if layout == "refgene_tss":
            points_path = out_dir / f"{name}.points.tsv"
            annotated = write_tss_points(points_path, intervals, variants)
            share = annotated / len(variants)
            track_reports.append({
                "track": name, "status": "OK", "kind": "points",
                "transcription_starts": len(intervals),
                "variants_annotated": annotated, "coverage": round(share, 6),
                "source_sha256": sha256_file(source),
            })
            if share >= args.min_coverage:
                manifest_rows.append({
                    "track": name, "path": str(points_path), "kind": "points",
                    "value_column": "score", "phenotype_free": "true",
                    "source": str(entry["url"] or source),
                    "justification": str(entry["why_phenotype_free"]),
                })
            continue

        groups: dict[str, list] = defaultdict(list)
        if entry["split_by_name"]:
            for row in intervals:
                groups[safe_name(row[3])].append(row)
        else:
            groups[name] = intervals

        for suffix, rows in sorted(groups.items()):
            track_name = name if suffix == name else f"{name}__{suffix}"
            bed_path = out_dir / f"{track_name}.bed"
            write_bed(bed_path, rows)
            annotated = coverage_of_bed(rows, variants)
            share = annotated / len(variants)
            track_reports.append({
                "track": track_name, "status": "OK", "kind": "bed",
                "intervals": len(rows), "variants_annotated": annotated,
                "coverage": round(share, 6), "retained": share >= args.min_coverage,
                "source_sha256": sha256_file(source),
            })
            if share >= args.min_coverage:
                manifest_rows.append({
                    "track": track_name, "path": str(bed_path), "kind": "bed",
                    "value_column": "score", "phenotype_free": "true",
                    "source": str(entry["url"] or source),
                    "justification": str(entry["why_phenotype_free"]),
                })

    if not manifest_rows:
        raise SystemExit(
            "No track clears --min-coverage. The usual cause is a genome-build mismatch: "
            "the panel is GRCh37 and a GRCh38 track hits nothing."
        )

    fields = ["track", "path", "kind", "value_column", "phenotype_free", "source",
              "justification"]
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv_writer(handle, fields)
        writer.writeheader()
        writer.writerows(manifest_rows)

    summary = {
        "classification": "COMPLEMENTARITY_ANNOTATION_PREP",
        "genome_build": "GRCh37/hg19, no lift-over performed",
        "phenotype_free_rail": (
            "every track is chromatin, accessibility, conservation or gene structure; "
            "association-derived names are refused outright"
        ),
        "panel_variants": len(variants),
        "tracks_requested": wanted,
        "tracks_in_manifest": [row["track"] for row in manifest_rows],
        "min_coverage": args.min_coverage,
        "reports": track_reports,
        "caveats": {name: CATALOGUE[name]["caveat"] for name in wanted if name in CATALOGUE},
        "manifest": str(manifest_path),
        "manifest_sha256": sha256_file(manifest_path),
    }
    (out_dir / "ANNOTATION_SUMMARY.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "status": "OK",
        "tracks_in_manifest": len(manifest_rows),
        "coverage": {r["track"]: r.get("coverage") for r in track_reports},
        "manifest": str(manifest_path),
    }, ensure_ascii=False, indent=2))
    return 0


def csv_writer(handle, fields):
    import csv

    return csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")


if __name__ == "__main__":
    raise SystemExit(main())
