#!/usr/bin/env python3
"""R1 read-only UKB asset inventory for the E0 foundation-v0 protocol.

This script never writes inside a scanned root, never reads phenotype values, and
never emits participant identifiers. It answers exactly the R1 questions:
cohort/build/phasing candidates, sample counts, variant counts, kinship source and
connected-component structure, disk paths, sizes and a hash plan.

Formal binding still requires a human to write the confirmed fields into
CONTRACT.json; this script only produces evidence and a draft.
"""

from __future__ import annotations

import argparse
import getpass
import gzip
import hashlib
import json
import os
import platform
import shutil
import socket
import struct
import subprocess
import sys
import time
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Iterable, Iterator

SCHEMA_VERSION = "r1_inventory/1.0.0"
PROTOCOL_ID = "20260916_e0_foundation_v0_rc0"

HEAD_TAIL_BYTES = 64 << 20
READ_CHUNK = 1 << 20

# GRCh37 vs GRCh38 contig lengths, used only to emit a build *hint*.
CONTIG_LENGTHS = {
    "GRCh37": {
        "1": 249250621, "2": 243199373, "3": 198022430, "4": 191154276,
        "5": 180915260, "6": 171115067, "7": 159138663, "8": 146364022,
        "9": 141213431, "10": 135534747, "11": 135006516, "12": 133851895,
        "13": 115169878, "14": 107349540, "15": 102531392, "16": 90354753,
        "17": 81195210, "18": 78077248, "19": 59128983, "20": 63025520,
        "21": 48129895, "22": 51304566, "X": 155270560, "Y": 59373566,
    },
    "GRCh38": {
        "1": 248956422, "2": 242193529, "3": 198295559, "4": 190214555,
        "5": 181538259, "6": 170805979, "7": 159345973, "8": 145138636,
        "9": 138394717, "10": 133797422, "11": 135086622, "12": 133275309,
        "13": 114364328, "14": 107043718, "15": 101991189, "16": 90338345,
        "17": 83257441, "18": 80373285, "19": 58617616, "20": 64444167,
        "21": 46709983, "22": 50818468, "X": 156040895, "Y": 57227415,
    },
}

KINSHIP_DEGREE_THRESHOLDS = {
    "third_degree_or_closer_0.0442": 0.0442,
    "second_degree_or_closer_0.0884": 0.0884,
    "first_degree_or_closer_0.177": 0.177,
    "mz_or_duplicate_0.354": 0.354,
}


# --------------------------------------------------------------------------- io


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8-sig", errors="replace", newline="")


def sha256_full(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(READ_CHUNK), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_head_tail(path: Path, size: int, window: int = HEAD_TAIL_BYTES) -> str:
    digest = hashlib.sha256()
    digest.update(str(size).encode("ascii"))
    with path.open("rb") as handle:
        remaining = min(window, size)
        while remaining > 0:
            chunk = handle.read(min(READ_CHUNK, remaining))
            if not chunk:
                break
            digest.update(chunk)
            remaining -= len(chunk)
        if size > window:
            handle.seek(max(size - window, window))
            while True:
                chunk = handle.read(READ_CHUNK)
                if not chunk:
                    break
                digest.update(chunk)
    return digest.hexdigest()


def digest_file(path: Path, size: int, full_hash_max_bytes: int) -> dict[str, object]:
    started = time.monotonic()
    if size <= full_hash_max_bytes:
        value, algorithm, complete = sha256_full(path), "sha256_full", True
    else:
        value = sha256_head_tail(path, size)
        algorithm, complete = "sha256_head_tail_size", False
    return {
        "algorithm": algorithm,
        "value": value,
        "covers_whole_file": complete,
        "window_bytes": None if complete else HEAD_TAIL_BYTES,
        "seconds": round(time.monotonic() - started, 3),
    }


def count_lines(path: Path, max_lines: int) -> tuple[int, bool]:
    count = 0
    with open_text(path) as handle:
        for _ in handle:
            count += 1
            if count >= max_lines:
                return count, True
    return count, False


# ------------------------------------------------------------------ classifying


GENOTYPE_GROUPS = {
    ".bed": ("plink1_bed", ".bed"),
    ".bim": ("plink1_bed", ".bed"),
    ".fam": ("plink1_bed", ".bed"),
    ".pgen": ("plink2_pgen", ".pgen"),
    ".pvar": ("plink2_pgen", ".pgen"),
    ".psam": ("plink2_pgen", ".pgen"),
    ".bgen": ("bgen", ".bgen"),
    ".sample": ("bgen", ".bgen"),
}

VCF_SUFFIXES = (".vcf.gz", ".vcf", ".bcf")
INDEX_SUFFIXES = (".bgi", ".tbi", ".csi", ".idx")


def is_kinship_candidate(name: str) -> bool:
    lowered = name.lower()
    if lowered.endswith(".kin0") or lowered.endswith(".kin"):
        return True
    if "rel" in lowered and lowered.endswith(".dat"):
        return True
    return ("king" in lowered or "kinship" in lowered) and lowered.endswith(
        (".dat", ".txt", ".tsv", ".csv", ".gz")
    )


def is_genetic_map_candidate(name: str) -> bool:
    lowered = name.lower()
    if "genetic_map" in lowered or "geneticmap" in lowered:
        return True
    return lowered.endswith(".map") and ("plink" in lowered or "chr" in lowered)


def is_sample_qc_candidate(name: str) -> bool:
    lowered = name.lower()
    return any(token in lowered for token in ("sqc", "sample_qc", "covar", "pca", "pcs"))


def looks_imputed(name: str) -> bool:
    lowered = name.lower()
    return any(token in lowered for token in ("imp", "mfi", "dosage", "dose"))


def looks_haplotype(name: str) -> bool:
    lowered = name.lower()
    return any(token in lowered for token in ("hap", "phase", "shapeit"))


# -------------------------------------------------------------------- parsing


def read_bgen_header(path: Path) -> dict[str, object]:
    with path.open("rb") as handle:
        prefix = handle.read(20)
        if len(prefix) < 20:
            return {"parsed": False, "reason": "file shorter than bgen header"}
        offset, header_length, variant_count, sample_count = struct.unpack("<IIII", prefix[:16])
        magic = prefix[16:20]
        if header_length < 20 or header_length > offset:
            return {"parsed": False, "reason": "implausible bgen header length"}
        handle.seek(4 + header_length - 4)
        flags_raw = handle.read(4)
        if len(flags_raw) < 4:
            return {"parsed": False, "reason": "truncated bgen flags"}
        flags = struct.unpack("<I", flags_raw)[0]
    compression = {0: "none", 1: "zlib", 2: "zstd"}.get(flags & 0x3, "unknown")
    return {
        "parsed": True,
        "magic": magic.decode("ascii", errors="replace"),
        "variant_count": int(variant_count),
        "sample_count": int(sample_count),
        "compression": compression,
        "layout": int((flags >> 2) & 0xF),
        "sample_identifiers_present": bool((flags >> 31) & 1),
    }


def parse_bim(path: Path, max_lines: int) -> dict[str, object]:
    chroms: Counter[str] = Counter()
    extent: dict[str, list[int]] = {}
    snp_like = 0
    total = 0
    truncated = False
    with open_text(path) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 6:
                continue
            total += 1
            if total > max_lines:
                truncated = True
                break
            chrom = fields[0]
            try:
                position = int(fields[3])
            except ValueError:
                continue
            chroms[chrom] += 1
            bounds = extent.setdefault(chrom, [position, position])
            bounds[0] = min(bounds[0], position)
            bounds[1] = max(bounds[1], position)
            if len(fields[4]) == 1 and len(fields[5]) == 1:
                snp_like += 1
    return {
        "variant_rows_scanned": total if not truncated else max_lines,
        "scan_truncated": truncated,
        "chromosomes": dict(sorted(chroms.items())),
        "chromosome_extent_bp": {key: {"min": value[0], "max": value[1]} for key, value in sorted(extent.items())},
        "biallelic_single_character_allele_rows": snp_like,
    }


def parse_pvar(path: Path, max_lines: int) -> dict[str, object]:
    chroms: Counter[str] = Counter()
    extent: dict[str, list[int]] = {}
    total = 0
    truncated = False
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 2:
                continue
            total += 1
            if total > max_lines:
                truncated = True
                break
            chrom = fields[0]
            try:
                position = int(fields[1])
            except ValueError:
                continue
            chroms[chrom] += 1
            bounds = extent.setdefault(chrom, [position, position])
            bounds[0] = min(bounds[0], position)
            bounds[1] = max(bounds[1], position)
    return {
        "variant_rows_scanned": total if not truncated else max_lines,
        "scan_truncated": truncated,
        "chromosomes": dict(sorted(chroms.items())),
        "chromosome_extent_bp": {key: {"min": value[0], "max": value[1]} for key, value in sorted(extent.items())},
    }


def parse_fam(path: Path, max_lines: int) -> dict[str, object]:
    total = 0
    negative_ids = 0
    sex_counts: Counter[str] = Counter()
    duplicate_probe: set[str] = set()
    duplicates = 0
    with open_text(path) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 2:
                continue
            total += 1
            if total > max_lines:
                break
            individual = fields[1]
            if individual.startswith("-"):
                negative_ids += 1
            if individual in duplicate_probe:
                duplicates += 1
            else:
                duplicate_probe.add(individual)
            if len(fields) >= 5:
                sex_counts[fields[4]] += 1
    return {
        "sample_rows": total,
        "negative_or_withdrawn_style_ids": negative_ids,
        "duplicate_individual_ids": duplicates,
        "sex_column_counts": dict(sorted(sex_counts.items())),
    }


def parse_psam(path: Path, max_lines: int) -> dict[str, object]:
    total = 0
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                continue
            if line.strip():
                total += 1
                if total > max_lines:
                    break
    return {"sample_rows": total}


def parse_oxford_sample(path: Path, max_lines: int) -> dict[str, object]:
    rows, _ = count_lines(path, max_lines)
    return {"sample_rows": max(rows - 2, 0), "raw_line_count": rows}


# ------------------------------------------------------------------- kinship


class UnionFind:
    def __init__(self) -> None:
        self.parent: dict[str, str] = {}
        self.rank: dict[str, int] = {}

    def add(self, value: str) -> None:
        if value not in self.parent:
            self.parent[value] = value
            self.rank[value] = 0

    def find(self, value: str) -> str:
        self.add(value)
        root = value
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[value] != root:
            self.parent[value], value = root, self.parent[value]
        return root

    def union(self, left: str, right: str) -> None:
        root_left, root_right = self.find(left), self.find(right)
        if root_left == root_right:
            return
        if self.rank[root_left] < self.rank[root_right]:
            root_left, root_right = root_right, root_left
        self.parent[root_right] = root_left
        if self.rank[root_left] == self.rank[root_right]:
            self.rank[root_left] += 1


def detect_kinship_columns(header: list[str]) -> dict[str, object]:
    lowered = [token.lower() for token in header]
    id_pairs = [("id1", "id2"), ("iid1", "iid2"), ("sample1", "sample2"), ("i1", "i2")]
    id_index: tuple[int, int] | None = None
    for left, right in id_pairs:
        if left in lowered and right in lowered:
            id_index = (lowered.index(left), lowered.index(right))
            break
    value_index = None
    for candidate in ("kinship", "kin", "phi", "pi_hat", "pihat"):
        if candidate in lowered:
            value_index = lowered.index(candidate)
            break
    return {
        "header_detected": id_index is not None,
        "id_columns": list(id_index) if id_index else None,
        "kinship_column": value_index,
        "header_fields": header,
    }


def summarize_kinship(path: Path, max_lines: int) -> dict[str, object]:
    with open_text(path) as handle:
        first = handle.readline()
    header = first.split()
    detection = detect_kinship_columns(header)
    if detection["header_detected"]:
        id_left, id_right = detection["id_columns"]  # type: ignore[misc]
        value_column = detection["kinship_column"]
        skip_first = True
    else:
        id_left, id_right, value_column = 0, 1, None
        skip_first = False
        detection["fallback"] = "assumed first two columns are the related pair"

    if value_column is None:
        return {
            "parsed": False,
            "reason": "no kinship value column detected; bind columns explicitly before R2",
            "column_detection": detection,
        }

    union = UnionFind()
    degree_counts = {name: 0 for name in KINSHIP_DEGREE_THRESHOLDS}
    rows = 0
    unparsable = 0
    minimum = None
    maximum = None
    truncated = False
    with open_text(path) as handle:
        for index, line in enumerate(handle):
            if index == 0 and skip_first:
                continue
            fields = line.split()
            if len(fields) <= max(id_left, id_right, value_column):
                unparsable += 1
                continue
            try:
                value = float(fields[value_column])
            except ValueError:
                unparsable += 1
                continue
            rows += 1
            if rows > max_lines:
                truncated = True
                break
            minimum = value if minimum is None else min(minimum, value)
            maximum = value if maximum is None else max(maximum, value)
            for name, threshold in KINSHIP_DEGREE_THRESHOLDS.items():
                if value >= threshold:
                    degree_counts[name] += 1
            if value >= KINSHIP_DEGREE_THRESHOLDS["third_degree_or_closer_0.0442"]:
                union.union(fields[id_left], fields[id_right])

    component_members: Counter[str] = Counter()
    for member in list(union.parent):
        component_members[union.find(member)] += 1
    sizes = Counter(component_members.values())
    return {
        "parsed": True,
        "column_detection": detection,
        "pair_rows": rows,
        "scan_truncated": truncated,
        "unparsable_rows": unparsable,
        "kinship_min": minimum,
        "kinship_max": maximum,
        "pairs_at_or_above": degree_counts,
        "components_at_0.0442": {
            "individuals_in_any_related_pair": len(union.parent),
            "multi_member_component_count": len(component_members),
            "largest_component_size": max(component_members.values()) if component_members else 0,
            "component_size_histogram": dict(sorted(sizes.items())),
        },
        "note": (
            "Singleton individuals are not represented here; the full component set is "
            "formed in R2 by unioning this edge list with the frozen sample list."
        ),
    }


# --------------------------------------------------------------------- scanning


@dataclass
class FileRecord:
    path: str
    root: str
    name: str
    size_bytes: int
    modified_utc: str
    category: str
    digest: dict[str, object] = field(default_factory=dict)
    parsed: dict[str, object] = field(default_factory=dict)


def iter_files(root: Path, max_depth: int, max_files: int) -> Iterator[Path]:
    root = root.resolve()
    yielded = 0
    for dirpath, dirnames, filenames in os.walk(root, topdown=True, followlinks=False):
        depth = len(Path(dirpath).resolve().relative_to(root).parts)
        if depth >= max_depth:
            dirnames[:] = []
        dirnames.sort()
        for filename in sorted(filenames):
            yield Path(dirpath) / filename
            yielded += 1
            if yielded >= max_files:
                return


def categorize(path: Path) -> str:
    name = path.name
    lowered = name.lower()
    suffix = path.suffix.lower()
    if lowered.endswith(VCF_SUFFIXES):
        return "vcf"
    if suffix in INDEX_SUFFIXES:
        return "index"
    if suffix in GENOTYPE_GROUPS:
        return GENOTYPE_GROUPS[suffix][0]
    if is_kinship_candidate(name):
        return "kinship"
    if is_genetic_map_candidate(name):
        return "genetic_map"
    if is_sample_qc_candidate(name):
        return "sample_qc"
    return "other"


def fileset_key(path: Path, category: str) -> str:
    if category == "vcf":
        for suffix in VCF_SUFFIXES:
            if path.name.lower().endswith(suffix):
                return str(path)[: -len(suffix)]
    return str(path.with_suffix(""))


def build_hint(chromosome_extent: dict[str, dict[str, int]]) -> dict[str, object]:
    votes: Counter[str] = Counter()
    evidence: list[dict[str, object]] = []
    for chrom, bounds in chromosome_extent.items():
        key = chrom.replace("chr", "").replace("Chr", "").upper()
        maximum = bounds["max"]
        fits = [
            build
            for build, lengths in CONTIG_LENGTHS.items()
            if key in lengths and maximum <= lengths[key]
        ]
        if len(fits) == 1:
            votes[fits[0]] += 1
            evidence.append({"chromosome": key, "max_position": maximum, "consistent_with": fits[0]})
    if not votes:
        return {"hint": None, "basis": "no discriminating chromosome extent", "evidence": evidence}
    ranked = votes.most_common()
    if len(ranked) > 1 and ranked[0][1] == ranked[1][1]:
        return {"hint": None, "basis": "ambiguous", "votes": dict(votes), "evidence": evidence[:10]}
    return {
        "hint": ranked[0][0],
        "basis": "max variant position vs contig lengths; HINT ONLY, must be bound explicitly",
        "votes": dict(votes),
        "evidence": evidence[:10],
    }


ARRAY_VARIANT_RANGE = (100_000, 3_000_000)
IMPUTED_VARIANT_MIN = 5_000_000


def e0_eligibility(
    group: dict[str, object], array_range: tuple[int, int] = ARRAY_VARIANT_RANGE
) -> dict[str, object]:
    variants = group.get("variant_count")
    name = str(group.get("fileset", ""))
    reasons: list[str] = []
    role = "unknown"
    if looks_haplotype(name):
        role = "phased_haplotype_candidate"
        reasons.append("name suggests phased haplotypes; usable as M1b reference, not as E0 target panel")
    elif looks_imputed(name) or (
        isinstance(variants, int) and variants > max(array_range[1], IMPUTED_VARIANT_MIN)
    ):
        role = "imputed_candidate"
        reasons.append(
            "imputed dosages are model output, not observed genotype; masking them would make "
            "E0 targets and the Beagle baseline circular"
        )
    elif isinstance(variants, int) and array_range[0] <= variants <= array_range[1]:
        role = "directly_genotyped_candidate"
        reasons.append("variant count is in the array-genotyping range; eligible as E0 primary panel")
    elif isinstance(variants, int):
        role = "out_of_expected_range"
        reasons.append("variant count outside both array and imputation ranges; inspect manually")
    else:
        reasons.append("variant count not established from headers; cannot classify")
    return {
        "role": role,
        "e0_primary_target_eligible": role == "directly_genotyped_candidate",
        "reasons": reasons,
    }


# ----------------------------------------------------------------- environment


def tool_version(command: list[str]) -> dict[str, object]:
    executable = shutil.which(command[0])
    if executable is None:
        return {"available": False, "path": None, "version": None}
    try:
        completed = subprocess.run(
            command, capture_output=True, text=True, timeout=30, check=False
        )
        output = (completed.stdout or completed.stderr or "").strip().splitlines()
        return {"available": True, "path": executable, "version": output[0] if output else ""}
    except (OSError, subprocess.SubprocessError) as error:
        return {"available": True, "path": executable, "version": None, "error": str(error)}


def probe_torch() -> dict[str, object]:
    try:
        import torch  # type: ignore
    except Exception as error:  # noqa: BLE001 - torch is optional at R1
        return {"importable": False, "reason": type(error).__name__}
    info: dict[str, object] = {
        "importable": True,
        "version": getattr(torch, "__version__", None),
        "cuda_built": getattr(getattr(torch, "version", None), "cuda", None),
        "cuda_available": False,
        "devices": [],
    }
    try:
        if torch.cuda.is_available():
            info["cuda_available"] = True
            info["devices"] = [
                {
                    "name": torch.cuda.get_device_name(index),
                    "total_memory_bytes": torch.cuda.get_device_properties(index).total_memory,
                }
                for index in range(torch.cuda.device_count())
            ]
    except Exception as error:  # noqa: BLE001
        info["cuda_probe_error"] = type(error).__name__
    return info


def probe_environment(out_dir: Path) -> dict[str, object]:
    usage = shutil.disk_usage(out_dir)
    return {
        "hostname": socket.gethostname(),
        "user": getpass.getuser(),
        "platform": platform.platform(),
        "python_version": sys.version.split()[0],
        "python_executable": sys.executable,
        "cpu_count": os.cpu_count(),
        "output_disk": {
            "path": str(out_dir),
            "total_bytes": usage.total,
            "free_bytes": usage.free,
            "free_gib": round(usage.free / (1 << 30), 2),
        },
        "tools": {
            "plink2": tool_version(["plink2", "--version"]),
            "plink": tool_version(["plink", "--version"]),
            "bcftools": tool_version(["bcftools", "--version"]),
            "bgenix": tool_version(["bgenix", "-help"]),
            "qctool": tool_version(["qctool", "-help"]),
            "java": tool_version(["java", "-version"]),
            "nvidia_smi": tool_version(["nvidia-smi", "--query-gpu=name,memory.total,driver_version", "--format=csv,noheader"]),
        },
        "torch": probe_torch(),
    }


# ----------------------------------------------------------------------- main


def scan_root(
    root: Path,
    *,
    max_depth: int,
    max_files: int,
    full_hash_max_bytes: int,
    max_scan_lines: int,
    skip_hashes: bool,
) -> list[FileRecord]:
    records: list[FileRecord] = []
    for path in iter_files(root, max_depth, max_files):
        try:
            stat = path.stat()
        except OSError:
            continue
        if not path.is_file():
            continue
        category = categorize(path)
        if category == "other":
            continue
        record = FileRecord(
            path=str(path),
            root=str(root),
            name=path.name,
            size_bytes=stat.st_size,
            modified_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(stat.st_mtime)),
            category=category,
        )
        if not os.access(path, os.R_OK):
            record.parsed = {"readable": False}
            records.append(record)
            continue
        if not skip_hashes:
            try:
                record.digest = digest_file(path, stat.st_size, full_hash_max_bytes)
            except OSError as error:
                record.digest = {"error": str(error)}
        suffix = path.suffix.lower()
        try:
            if suffix == ".bim":
                record.parsed = parse_bim(path, max_scan_lines)
            elif suffix == ".pvar":
                record.parsed = parse_pvar(path, max_scan_lines)
            elif suffix == ".fam":
                record.parsed = parse_fam(path, max_scan_lines)
            elif suffix == ".psam":
                record.parsed = parse_psam(path, max_scan_lines)
            elif suffix == ".sample":
                record.parsed = parse_oxford_sample(path, max_scan_lines)
            elif suffix == ".bgen":
                record.parsed = read_bgen_header(path)
            elif suffix == ".bed":
                with path.open("rb") as handle:
                    record.parsed = {"plink1_magic_ok": handle.read(3) == b"\x6c\x1b\x01"}
            elif category == "kinship":
                record.parsed = summarize_kinship(path, max_scan_lines)
        except (OSError, ValueError, struct.error) as error:
            record.parsed = {"parse_error": f"{type(error).__name__}: {error}"}
        records.append(record)
    return records


def group_filesets(
    records: Iterable[FileRecord], array_range: tuple[int, int] = ARRAY_VARIANT_RANGE
) -> list[dict[str, object]]:
    groups: dict[str, dict[str, object]] = {}
    for record in records:
        if record.category not in {"plink1_bed", "plink2_pgen", "bgen", "vcf"}:
            continue
        key = fileset_key(Path(record.path), record.category)
        group = groups.setdefault(
            key,
            {
                "fileset": key,
                "format": record.category,
                "files": [],
                "total_bytes": 0,
                "sample_count": None,
                "variant_count": None,
                "chromosome_extent_bp": {},
            },
        )
        group["files"].append({
            "name": record.name,
            "size_bytes": record.size_bytes,
            "digest": record.digest,
        })
        group["total_bytes"] = int(group["total_bytes"]) + record.size_bytes
        parsed = record.parsed or {}
        if "sample_rows" in parsed:
            group["sample_count"] = parsed["sample_rows"]
        if parsed.get("parsed") and "sample_count" in parsed:
            group["sample_count"] = parsed["sample_count"]
            group["variant_count"] = parsed["variant_count"]
            group["bgen_header"] = {
                key_: parsed[key_]
                for key_ in ("compression", "layout", "sample_identifiers_present", "magic")
                if key_ in parsed
            }
        if "variant_rows_scanned" in parsed:
            if not parsed.get("scan_truncated"):
                group["variant_count"] = parsed["variant_rows_scanned"]
            else:
                group["variant_count_lower_bound"] = parsed["variant_rows_scanned"]
            group["chromosome_extent_bp"] = parsed.get("chromosome_extent_bp", {})
            group["chromosomes"] = parsed.get("chromosomes", {})
    result = []
    for group in groups.values():
        group["files"].sort(key=lambda item: item["name"])
        group["build_hint"] = build_hint(group.get("chromosome_extent_bp") or {})
        group["e0_eligibility"] = e0_eligibility(group, array_range)
        result.append(group)
    result.sort(key=lambda item: str(item["fileset"]))
    return result


def collect_blocking(payload: dict[str, object]) -> list[str]:
    blocking: list[str] = []
    filesets = payload["genotype_filesets"]
    if not filesets:
        blocking.append("no genotype fileset found under the scanned roots")
    if not any(item["e0_eligibility"]["e0_primary_target_eligible"] for item in filesets):
        blocking.append(
            "no directly genotyped (array-scale) fileset identified; E0 primary targets must be "
            "observed calls, not imputed dosages"
        )
    kinship = [item for item in payload["kinship_files"] if item.get("parsed", {}).get("parsed")]
    if not kinship:
        blocking.append("no parsable kinship/relatedness file found; family components cannot be frozen")
    if not payload["genetic_map_files"]:
        blocking.append("no genetic map found; the 1 cM guard rule cannot be applied")
    environment = payload["environment"]
    if not environment["torch"].get("cuda_available"):
        blocking.append("no CUDA-enabled PyTorch detected; micro-batch/precision cannot be frozen")
    if not environment["tools"]["java"]["available"]:
        blocking.append("no Java runtime detected; Beagle 5.5 M1b adapter cannot run")
    builds = {
        item["build_hint"]["hint"]
        for item in filesets
        if item["build_hint"].get("hint")
    }
    if len(builds) > 1:
        blocking.append(f"conflicting genome build hints across filesets: {sorted(builds)}")
    return blocking


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, action="append", required=True,
                        help="Directory to scan read-only; repeatable.")
    parser.add_argument("--out-dir", type=Path, required=True,
                        help="Output directory; must not live inside any scanned root.")
    parser.add_argument("--cohort-id", default=None, help="Operator-declared cohort identifier.")
    parser.add_argument("--authorization-id", default=None,
                        help="Operator-declared application/authorization identifier.")
    parser.add_argument("--access-mode", default=None,
                        choices=[None, "local_institutional_copy", "ukb_rap_dnanexus", "other"],
                        help="How this server reaches the cohort.")
    parser.add_argument("--max-depth", type=int, default=6)
    parser.add_argument("--max-files", type=int, default=200_000)
    parser.add_argument("--max-scan-lines", type=int, default=5_000_000)
    parser.add_argument("--full-hash-max-bytes", type=int, default=2 << 30,
                        help="Files larger than this get a head/tail/size fingerprint instead.")
    parser.add_argument("--array-variant-min", type=int, default=ARRAY_VARIANT_RANGE[0],
                        help="Lower variant count for an array-scale (directly genotyped) panel.")
    parser.add_argument("--array-variant-max", type=int, default=ARRAY_VARIANT_RANGE[1],
                        help="Upper variant count for an array-scale panel.")
    parser.add_argument("--skip-hashes", action="store_true",
                        help="Skip digests entirely for a first fast pass.")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    out_dir = args.out_dir.resolve()
    roots = [root.resolve() for root in args.root]

    for root in roots:
        if not root.is_dir():
            raise SystemExit(f"Root is not a directory: {root}")
        if out_dir == root or root in out_dir.parents:
            raise SystemExit(
                f"Refusing to write inventory inside a scanned root: out-dir={out_dir}, root={root}"
            )
    out_dir.mkdir(parents=True, exist_ok=True)
    inventory_path = out_dir / "INVENTORY.json"
    if inventory_path.exists() and not args.overwrite:
        raise SystemExit(f"Refusing to overwrite: {inventory_path}")

    started = time.time()
    records: list[FileRecord] = []
    for root in roots:
        records.extend(
            scan_root(
                root,
                max_depth=args.max_depth,
                max_files=args.max_files,
                full_hash_max_bytes=args.full_hash_max_bytes,
                max_scan_lines=args.max_scan_lines,
                skip_hashes=args.skip_hashes,
            )
        )

    by_category: dict[str, list[dict[str, object]]] = defaultdict(list)
    for record in records:
        by_category[record.category].append(asdict(record))

    payload: dict[str, object] = {
        "schema_version": SCHEMA_VERSION,
        "classification": "R1_READ_ONLY_INVENTORY",
        "protocol_id": PROTOCOL_ID,
        "milestone": "R1",
        "identifiers_included": False,
        "phenotype_read": False,
        "test_data_read": False,
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(started)),
        "elapsed_seconds": None,
        "operator_declared": {
            "cohort_id": args.cohort_id,
            "authorization_id": args.authorization_id,
            "access_mode": args.access_mode,
            "roots": [str(root) for root in roots],
        },
        "scan_settings": {
            "max_depth": args.max_depth,
            "max_files": args.max_files,
            "max_scan_lines": args.max_scan_lines,
            "full_hash_max_bytes": args.full_hash_max_bytes,
            "skip_hashes": args.skip_hashes,
            "array_variant_range": [args.array_variant_min, args.array_variant_max],
        },
        "environment": probe_environment(out_dir),
        "file_counts_by_category": {key: len(value) for key, value in sorted(by_category.items())},
        "genotype_filesets": group_filesets(records, (args.array_variant_min, args.array_variant_max)),
        "kinship_files": by_category.get("kinship", []),
        "genetic_map_files": by_category.get("genetic_map", []),
        "sample_qc_files": by_category.get("sample_qc", []),
        "vcf_files": by_category.get("vcf", []),
    }
    payload["blocking"] = collect_blocking(payload)
    payload["elapsed_seconds"] = round(time.time() - started, 2)

    inventory_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    hash_plan = out_dir / "HASH_PLAN.tsv"
    with hash_plan.open("w", encoding="utf-8", newline="") as handle:
        handle.write("path\tsize_bytes\tcategory\talgorithm\tcovers_whole_file\tdigest\n")
        for record in sorted(records, key=lambda item: item.path):
            digest = record.digest or {}
            handle.write(
                "\t".join(
                    [
                        record.path,
                        str(record.size_bytes),
                        record.category,
                        str(digest.get("algorithm", "none")),
                        str(digest.get("covers_whole_file", False)),
                        str(digest.get("value", "")),
                    ]
                )
                + "\n"
            )

    summary_path = out_dir / "INVENTORY_SUMMARY.md"
    summary_path.write_text(render_summary(payload), encoding="utf-8")

    print(json.dumps(
        {
            "status": "OK",
            "inventory": str(inventory_path),
            "inventory_sha256": sha256_full(inventory_path),
            "hash_plan": str(hash_plan),
            "summary": str(summary_path),
            "genotype_fileset_count": len(payload["genotype_filesets"]),
            "blocking_count": len(payload["blocking"]),
        },
        ensure_ascii=False,
        indent=2,
    ))
    return 0


def render_summary(payload: dict[str, object]) -> str:
    lines = [
        "# R1 UKB read-only inventory summary",
        "",
        f"- 协议 `{payload['protocol_id']}` · 里程碑 `R1` · 生成 `{payload['generated_utc']}`",
        f"- 声明 cohort：`{payload['operator_declared']['cohort_id']}`；授权：`{payload['operator_declared']['authorization_id']}`",
        f"- 扫描根目录：{', '.join(f'`{item}`' for item in payload['operator_declared']['roots'])}",
        "- 本文件不含任何参与者标识，只含计数、路径、大小与哈希",
        "",
        "## Genotype filesets",
        "",
        "| fileset | format | samples | variants | build hint | E0 role | primary eligible |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for group in payload["genotype_filesets"]:  # type: ignore[index]
        eligibility = group["e0_eligibility"]
        lines.append(
            "| `{name}` | {fmt} | {samples} | {variants} | {build} | {role} | {eligible} |".format(
                name=Path(str(group["fileset"])).name,
                fmt=group["format"],
                samples=group.get("sample_count"),
                variants=group.get("variant_count") or group.get("variant_count_lower_bound"),
                build=group["build_hint"].get("hint"),
                role=eligibility["role"],
                eligible="✅" if eligibility["e0_primary_target_eligible"] else "❌",
            )
        )
    lines += ["", "## Kinship", ""]
    kinship_files = payload["kinship_files"]  # type: ignore[index]
    if not kinship_files:
        lines.append("- 未发现 kinship 文件")
    for item in kinship_files:
        parsed = item.get("parsed", {})
        if not parsed.get("parsed"):
            lines.append(f"- `{item['name']}`：未解析（{parsed.get('reason', parsed.get('parse_error'))}）")
            continue
        components = parsed["components_at_0.0442"]
        lines.append(
            f"- `{item['name']}`：{parsed['pair_rows']} 对；"
            f"≥0.0442 的对 {parsed['pairs_at_or_above']['third_degree_or_closer_0.0442']}；"
            f"涉及个体 {components['individuals_in_any_related_pair']}；"
            f"多成员 component {components['multi_member_component_count']}；"
            f"最大 component {components['largest_component_size']}"
        )
    environment = payload["environment"]  # type: ignore[index]
    lines += [
        "",
        "## Environment",
        "",
        f"- host `{environment['hostname']}` · {environment['platform']} · python {environment['python_version']}",
        f"- 输出盘可用 {environment['output_disk']['free_gib']} GiB",
        f"- torch：{environment['torch'].get('version') or '未安装'}，"
        f"CUDA 可用 {bool(environment['torch'].get('cuda_available'))}",
        f"- java：{environment['tools']['java']['available']}；plink2：{environment['tools']['plink2']['available']}；"
        f"bcftools：{environment['tools']['bcftools']['available']}",
        "",
        "## Blocking",
        "",
    ]
    blocking = payload["blocking"]  # type: ignore[index]
    if not blocking:
        lines.append("- 无（仍需人工确认 build、phasing 与授权字段后才能签发 R2）")
    for item in blocking:
        lines.append(f"- {item}")
    lines += [
        "",
        "> 本 inventory 是 R1 证据，不是 `RUN_AUTHORIZED`。build、phasing、kinship 版本与阈值必须由人写入 `CONTRACT.json` 后才进入 R2。",
        "",
    ]
    return "\n".join(lines)


if __name__ == "__main__":
    raise SystemExit(main())
