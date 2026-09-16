#!/usr/bin/env python3
"""Independently validate an R1 inventory before it is used to bind CONTRACT.json.

The validator re-reads INVENTORY.json only. It checks structural completeness,
identifier hygiene, digest coverage and the read-only boundary, and it recomputes
file sizes on disk so a stale inventory cannot be promoted to R2.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

REQUIRED_TOP_LEVEL = [
    "schema_version",
    "classification",
    "protocol_id",
    "milestone",
    "identifiers_included",
    "operator_declared",
    "environment",
    "genotype_filesets",
    "kinship_files",
    "genetic_map_files",
    "blocking",
]

FORBIDDEN_KEYS = {
    "sample_ids",
    "sample_id_list",
    "eid",
    "eids",
    "iids",
    "members",
    "individual_ids",
    "participant_ids",
}

ALLOWED_DIGEST_ALGORITHMS = {"sha256_full", "sha256_head_tail_size", "none"}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def walk_keys(node: object, path: str = "") -> list[tuple[str, object]]:
    found: list[tuple[str, object]] = []
    if isinstance(node, dict):
        for key, value in node.items():
            found.append((f"{path}.{key}" if path else str(key), value))
            found.extend(walk_keys(value, f"{path}.{key}" if path else str(key)))
    elif isinstance(node, list):
        for index, value in enumerate(node):
            found.extend(walk_keys(value, f"{path}[{index}]"))
    return found


def check_identifier_hygiene(payload: dict) -> list[str]:
    errors: list[str] = []
    if payload.get("identifiers_included") is not False:
        errors.append("identifiers_included must be false")
    for key_path, value in walk_keys(payload):
        leaf = key_path.split(".")[-1].split("[")[0]
        if leaf in FORBIDDEN_KEYS:
            errors.append(f"forbidden identifier-bearing key: {key_path}")
        if isinstance(value, list) and len(value) > 200 and all(isinstance(item, str) for item in value):
            errors.append(f"suspicious long string list (possible identifier dump): {key_path}")
    return errors


def check_read_only_boundary(payload: dict, inventory_path: Path) -> list[str]:
    errors: list[str] = []
    out_dir = inventory_path.resolve().parent
    for root in payload.get("operator_declared", {}).get("roots", []):
        root_path = Path(root).resolve()
        if out_dir == root_path or root_path in out_dir.parents:
            errors.append(f"inventory output lives inside scanned root {root}")
    return errors


def check_filesets(payload: dict, verify_sizes: bool) -> tuple[list[str], list[str], dict]:
    errors: list[str] = []
    warnings: list[str] = []
    stats = {
        "fileset_count": 0,
        "e0_primary_eligible_count": 0,
        "member_files_checked": 0,
        "size_mismatches": 0,
        "full_hash_coverage": 0,
        "partial_hash_coverage": 0,
    }
    for group in payload.get("genotype_filesets", []):
        stats["fileset_count"] += 1
        if group.get("e0_eligibility", {}).get("e0_primary_target_eligible"):
            stats["e0_primary_eligible_count"] += 1
        for member in group.get("files", []):
            stats["member_files_checked"] += 1
            digest = member.get("digest") or {}
            algorithm = digest.get("algorithm", "none")
            if algorithm not in ALLOWED_DIGEST_ALGORITHMS:
                errors.append(f"unknown digest algorithm {algorithm} for {member.get('name')}")
            if digest.get("covers_whole_file"):
                stats["full_hash_coverage"] += 1
            else:
                stats["partial_hash_coverage"] += 1
    if verify_sizes:
        for record in (
            payload.get("kinship_files", [])
            + payload.get("genetic_map_files", [])
            + payload.get("sample_qc_files", [])
        ):
            path = Path(record["path"])
            if not path.exists():
                errors.append(f"inventoried file no longer exists: {path}")
                continue
            if path.stat().st_size != record["size_bytes"]:
                stats["size_mismatches"] += 1
                errors.append(f"size changed since inventory: {path}")
    if stats["e0_primary_eligible_count"] == 0:
        warnings.append(
            "no directly genotyped fileset flagged eligible; E0 primary targets must be observed "
            "calls rather than imputed dosages"
        )
    return errors, warnings, stats


def check_kinship(payload: dict) -> tuple[list[str], list[str], dict]:
    errors: list[str] = []
    warnings: list[str] = []
    parsed = [item for item in payload.get("kinship_files", []) if item.get("parsed", {}).get("parsed")]
    summary = {
        "kinship_file_count": len(payload.get("kinship_files", [])),
        "parsable_kinship_file_count": len(parsed),
        "selected": None,
    }
    if not parsed:
        errors.append("no parsable kinship file; family components cannot be frozen in R2")
        return errors, warnings, summary
    if len(parsed) > 1:
        warnings.append(
            f"{len(parsed)} parsable kinship files found; exactly one must be bound in CONTRACT.json"
        )
    best = parsed[0]["parsed"]
    summary["selected"] = {
        "name": parsed[0]["name"],
        "pair_rows": best.get("pair_rows"),
        "scan_truncated": best.get("scan_truncated"),
        "pairs_at_or_above": best.get("pairs_at_or_above"),
        "components_at_0.0442": best.get("components_at_0.0442"),
    }
    if best.get("scan_truncated"):
        errors.append("kinship scan was truncated; rerun with a larger --max-scan-lines")
    if best.get("kinship_max") is not None and best["kinship_max"] > 0.5:
        warnings.append("kinship values above 0.5 present; confirm the coefficient convention")
    return errors, warnings, summary


def check_environment(payload: dict) -> list[str]:
    warnings: list[str] = []
    environment = payload.get("environment", {})
    if not environment.get("torch", {}).get("cuda_available"):
        warnings.append("CUDA-enabled PyTorch not verified; R3 budget fields stay unbound")
    if not environment.get("tools", {}).get("java", {}).get("available"):
        warnings.append("Java runtime absent; the M1b Beagle adapter cannot run yet")
    free_gib = environment.get("output_disk", {}).get("free_gib")
    if isinstance(free_gib, (int, float)) and free_gib < 50:
        warnings.append(f"output disk has only {free_gib} GiB free")
    return warnings


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--skip-size-verification", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    payload = json.loads(args.inventory.read_text(encoding="utf-8"))

    errors: list[str] = []
    warnings: list[str] = []

    for key in REQUIRED_TOP_LEVEL:
        if key not in payload:
            errors.append(f"missing required field: {key}")
    if payload.get("classification") != "R1_READ_ONLY_INVENTORY":
        errors.append("classification is not R1_READ_ONLY_INVENTORY")
    if payload.get("phenotype_read") is not False:
        errors.append("phenotype_read must be false at R1")
    if payload.get("test_data_read") is not False:
        errors.append("test_data_read must be false at R1")

    errors.extend(check_identifier_hygiene(payload))
    errors.extend(check_read_only_boundary(payload, args.inventory))

    fileset_errors, fileset_warnings, fileset_stats = check_filesets(
        payload, verify_sizes=not args.skip_size_verification
    )
    kinship_errors, kinship_warnings, kinship_summary = check_kinship(payload)
    errors.extend(fileset_errors + kinship_errors)
    warnings.extend(fileset_warnings + kinship_warnings + check_environment(payload))

    report = {
        "classification": "R1_INVENTORY_VALIDATION",
        "status": "PASS" if not errors else "FAIL",
        "identifiers_included": False,
        "inventory_path": str(args.inventory),
        "inventory_sha256": sha256_file(args.inventory),
        "filesets": fileset_stats,
        "kinship": kinship_summary,
        "inventory_blocking": payload.get("blocking", []),
        "errors": errors,
        "warnings": warnings,
        "r2_allowed": not errors,
        "run_authorized": False,
    }
    text = json.dumps(report, ensure_ascii=False, indent=2) + "\n"
    if args.output is not None:
        if args.output.exists() and not args.overwrite:
            raise SystemExit(f"Refusing to overwrite: {args.output}")
        args.output.write_text(text, encoding="utf-8")
    print(text, end="")
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
