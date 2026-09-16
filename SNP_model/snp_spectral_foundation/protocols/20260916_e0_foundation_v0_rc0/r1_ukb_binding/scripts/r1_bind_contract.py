#!/usr/bin/env python3
"""Produce a CONTRACT draft from R1 evidence without granting authorization.

Fields that a machine can establish (hashes, formats, counts, environment) are
filled from the inventory. Fields that require a human decision (genome build,
phasing status, cohort and authorization identifiers, kinship version) stay null
unless the operator passes them explicitly. run_authorized is always false here.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def select_fileset(inventory: dict, requested: str | None) -> dict | None:
    filesets = inventory.get("genotype_filesets", [])
    if requested:
        for group in filesets:
            if group["fileset"] == requested or Path(str(group["fileset"])).name == requested:
                return group
        raise SystemExit(f"Requested fileset not present in inventory: {requested}")
    eligible = [g for g in filesets if g.get("e0_eligibility", {}).get("e0_primary_target_eligible")]
    if len(eligible) == 1:
        return eligible[0]
    return None


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", type=Path, required=True)
    parser.add_argument("--template", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--fileset", default=None, help="Fileset path or basename to bind as the E0 panel.")
    parser.add_argument("--family-summary", type=Path, default=None)
    parser.add_argument("--genome-build", default=None, help="Explicit operator binding, e.g. GRCh37.")
    parser.add_argument("--phasing-status", default=None,
                        choices=[None, "unphased_calls", "phased_haplotypes", "mixed"])
    parser.add_argument("--kinship-source", default=None)
    parser.add_argument("--kinship-version", default=None)
    parser.add_argument("--genetic-map-source", default=None)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.output.exists() and not args.overwrite:
        raise SystemExit(f"Refusing to overwrite: {args.output}")

    inventory = json.loads(args.inventory.read_text(encoding="utf-8"))
    contract = json.loads(args.template.read_text(encoding="utf-8"))

    declared = inventory.get("operator_declared", {})
    environment = inventory.get("environment", {})
    fileset = select_fileset(inventory, args.fileset)

    contract["status"] = "R1_DRAFT_NOT_AUTHORIZED"
    data = contract.setdefault("data", {})
    data["cohort_id"] = declared.get("cohort_id")
    data["authorization_id"] = declared.get("authorization_id")
    data["genome_build"] = args.genome_build
    data["phasing_status"] = args.phasing_status
    data["kinship_source"] = args.kinship_source
    data["kinship_version"] = args.kinship_version
    data["genetic_map_source"] = args.genetic_map_source

    if fileset is not None:
        data["genotype_format"] = fileset.get("format")
        contract["r1_selected_fileset"] = {
            "fileset": fileset.get("fileset"),
            "format": fileset.get("format"),
            "sample_count": fileset.get("sample_count"),
            "variant_count": fileset.get("variant_count"),
            "file_digests": {
                item["name"]: item.get("digest", {}) for item in fileset.get("files", [])
            },
            "e0_eligibility": fileset.get("e0_eligibility"),
            "build_hint": fileset.get("build_hint"),
        }

    kinship_files = [item for item in inventory.get("kinship_files", []) if item.get("parsed", {}).get("parsed")]
    if len(kinship_files) == 1:
        record = kinship_files[0]
        contract["r1_selected_kinship"] = {
            "path": record["path"],
            "digest": record.get("digest", {}),
            "summary": record.get("parsed", {}),
        }

    map_files = inventory.get("genetic_map_files", [])
    if len(map_files) == 1:
        data["genetic_map_sha256"] = (map_files[0].get("digest") or {}).get("value")

    if args.family_summary is not None:
        family = json.loads(args.family_summary.read_text(encoding="utf-8"))
        contract.setdefault("split", {})["family_manifest_sha256"] = family.get("outputs", {}).get(
            "family_manifest_sha256"
        )
        contract["split"]["family_component_cross_split_count"] = family.get("components", {}).get(
            "cross_split_count"
        )
        contract["split"]["seed"] = family.get("seed")
        data["kinship_threshold"] = family.get("kinship", {}).get("threshold", data.get("kinship_threshold"))

    resources = contract.setdefault("resources", {})
    resources["cuda_runtime_verified"] = bool(environment.get("torch", {}).get("cuda_available"))
    resources["host"] = environment.get("hostname")
    resources["output_disk_free_gib"] = environment.get("output_disk", {}).get("free_gib")

    contract["r1_evidence"] = {
        "inventory_path": str(args.inventory),
        "inventory_sha256": sha256_file(args.inventory),
        "inventory_generated_utc": inventory.get("generated_utc"),
        "inventory_blocking": inventory.get("blocking", []),
    }

    required = {
        "data.cohort_id": data.get("cohort_id"),
        "data.authorization_id": data.get("authorization_id"),
        "data.genome_build": data.get("genome_build"),
        "data.genotype_format": data.get("genotype_format"),
        "data.phasing_status": data.get("phasing_status"),
        "data.kinship_source": data.get("kinship_source"),
        "data.kinship_version": data.get("kinship_version"),
        "data.genetic_map_source": data.get("genetic_map_source"),
        "data.genetic_map_sha256": data.get("genetic_map_sha256"),
        "split.family_manifest_sha256": contract.get("split", {}).get("family_manifest_sha256"),
        "split.block_manifest_sha256": contract.get("split", {}).get("block_manifest_sha256"),
    }
    unbound = sorted(name for name, value in required.items() if value in (None, ""))
    contract["r1_unbound_fields"] = unbound

    authorization = contract.setdefault("authorization", {})
    authorization["all_required_fields_bound"] = not unbound
    authorization["run_authorized"] = False

    args.output.write_text(json.dumps(contract, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(
        {
            "status": "OK",
            "draft": str(args.output),
            "draft_sha256": sha256_file(args.output),
            "unbound_field_count": len(unbound),
            "unbound_fields": unbound,
            "run_authorized": False,
        },
        ensure_ascii=False,
        indent=2,
    ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
