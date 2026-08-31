from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_first_column(path: Path, header_names: set[str]) -> list[str]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = [line.rstrip("\r\n") for line in handle if line.strip()]
    if not rows:
        raise ValueError(f"empty manifest: {path}")
    delimiter = "\t" if "\t" in rows[0] else "," if "," in rows[0] else None
    values = []
    for line in rows:
        value = line.split(delimiter, 1)[0].strip() if delimiter else line.strip()
        values.append(value)
    if values[0].lower() in header_names:
        values = values[1:]
    if not values or any(not value for value in values):
        raise ValueError(f"missing IDs in {path}")
    if len(values) != len(set(values)):
        raise ValueError(f"duplicate IDs in {path}")
    return values


def read_population_map(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample_id", "population", "cohort"}
        if not reader.fieldnames or not required <= set(reader.fieldnames):
            raise ValueError(f"{path} requires columns {sorted(required)}")
        result: dict[str, dict[str, str]] = {}
        for row in reader:
            sample_id = row["sample_id"].strip()
            if not sample_id or sample_id in result:
                raise ValueError(f"missing/duplicate sample_id in {path}: {sample_id!r}")
            result[sample_id] = {
                "population": row["population"].strip(),
                "cohort": row["cohort"].strip(),
            }
    return result


def file_record(path: Path) -> dict[str, Any]:
    resolved = path.resolve()
    stat = resolved.stat()
    return {
        "name": resolved.name,
        "bytes": stat.st_size,
        "sha256": sha256(resolved),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Build the frozen v7.1.0 Branch-B panel manifest")
    parser.add_argument("--panel-spec", type=Path, required=True)
    parser.add_argument("--joint-chr22-panel", type=Path, required=True)
    parser.add_argument("--variants", type=Path, required=True, help="One variant ID per line or first TSV/CSV column")
    parser.add_argument("--donor-train-samples", type=Path, required=True)
    parser.add_argument("--donor-validation-samples", type=Path, required=True)
    parser.add_argument("--hgdp-primary-samples", type=Path, required=True)
    parser.add_argument("--hgdp-snpbag-calibration-samples", type=Path, required=True)
    parser.add_argument("--sample-populations", type=Path, required=True)
    parser.add_argument("--hapnest-config", type=Path, required=True)
    parser.add_argument("--neighbor-classification", type=Path, required=True)
    parser.add_argument("--source-release-record", type=Path, required=True)
    parser.add_argument("--compatibility-note", default="none")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    input_paths = {
        "panel_spec": args.panel_spec,
        "joint_chr22_panel": args.joint_chr22_panel,
        "variants": args.variants,
        "donor_train_samples": args.donor_train_samples,
        "donor_validation_samples": args.donor_validation_samples,
        "hgdp_primary_samples": args.hgdp_primary_samples,
        "hgdp_snpbag_calibration_samples": args.hgdp_snpbag_calibration_samples,
        "sample_populations": args.sample_populations,
        "hapnest_config": args.hapnest_config,
        "neighbor_classification": args.neighbor_classification,
        "source_release_record": args.source_release_record,
    }
    missing = [str(path) for path in input_paths.values() if not path.is_file()]
    if missing:
        raise SystemExit("missing T01 inputs: " + ", ".join(missing))

    spec = json.loads(args.panel_spec.read_text(encoding="utf-8"))
    if spec.get("protocol_version") != "v7.1.0":
        raise SystemExit("panel spec must declare protocol_version v7.1.0")

    sample_headers = {"sample", "sample_id", "iid", "individual_id"}
    groups = {
        "donor_train": read_first_column(args.donor_train_samples, sample_headers),
        "donor_validation": read_first_column(args.donor_validation_samples, sample_headers),
        "hgdp_primary_holdout": read_first_column(args.hgdp_primary_samples, sample_headers),
        "hgdp_snpbag_calibration": read_first_column(args.hgdp_snpbag_calibration_samples, sample_headers),
    }
    sets = {name: set(values) for name, values in groups.items()}
    overlaps: dict[str, list[str]] = {}
    names = list(groups)
    for i, left in enumerate(names):
        for right in names[i + 1 :]:
            shared = sorted(sets[left] & sets[right])
            if shared:
                overlaps[f"{left}__{right}"] = shared[:20]
    if overlaps:
        raise SystemExit("partition overlap detected: " + json.dumps(overlaps, sort_keys=True))

    population_map = read_population_map(args.sample_populations)
    all_samples = set().union(*sets.values())
    missing_metadata = sorted(all_samples - set(population_map))
    if missing_metadata:
        raise SystemExit(f"samples missing population metadata: {missing_metadata[:20]}")
    for name in ("donor_train", "donor_validation"):
        wrong = sorted(s for s in sets[name] if population_map[s]["cohort"].upper() != "1KGP")
        if wrong:
            raise SystemExit(f"non-1KGP samples in {name}: {wrong[:20]}")
    for name in ("hgdp_primary_holdout", "hgdp_snpbag_calibration"):
        wrong = sorted(s for s in sets[name] if population_map[s]["cohort"].upper() != "HGDP")
        if wrong:
            raise SystemExit(f"non-HGDP samples in {name}: {wrong[:20]}")

    calibration = groups["hgdp_snpbag_calibration"]
    cal_counts: dict[str, int] = {}
    for sample in calibration:
        pop = population_map[sample]["population"]
        cal_counts[pop] = cal_counts.get(pop, 0) + 1
    if len(calibration) != 216 or len(cal_counts) != 54 or set(cal_counts.values()) != {4}:
        raise SystemExit("SNPBag calibration must contain exactly 216 samples: 54 populations x 4")

    variants = read_first_column(args.variants, {"variant", "variant_id", "id", "site_id"})
    expected_l = int(spec["site_selection"]["snpbag_compatibility_expected_L"])
    compatibility = "EXACT_EXPECTED_L" if len(variants) == expected_l else "DATA_DERIVED_L_DIFFERS"
    records = {role: file_record(path) for role, path in input_paths.items()}
    primary_populations = sorted({population_map[s]["population"] for s in groups["hgdp_primary_holdout"]})
    payload = {
        "schema_version": "1.0",
        "protocol_version": "v7.1.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "branch": "B",
        "status": "SIGNED_INPUTS_PRESENT",
        "generator": {"primary": "HAPNEST", "sensitivity": "msprime", "hgdp_used_as_donor": False},
        "chr22": {
            "variant_count": len(variants),
            "snpbag_expected_variant_count": expected_l,
            "compatibility_status": compatibility,
            "compatibility_note": args.compatibility_note,
            "length_was_forced": False,
        },
        "partitions": {
            name: {
                "individual_count": len(values),
                "sample_manifest_sha256": records[{
                    "donor_train": "donor_train_samples",
                    "donor_validation": "donor_validation_samples",
                    "hgdp_primary_holdout": "hgdp_primary_samples",
                    "hgdp_snpbag_calibration": "hgdp_snpbag_calibration_samples",
                }[name]]["sha256"],
            }
            for name, values in groups.items()
        },
        "hgdp_primary_population_count": len(primary_populations),
        "hgdp_primary_populations": primary_populations,
        "all_partitions_disjoint": True,
        "inputs": records,
        "claim_ceiling": "chr22 masked-genotype architecture comparison under frozen public-panel and generator distributions",
    }
    args.output_dir.mkdir(parents=True, exist_ok=False)
    output = args.output_dir / "PANEL_MANIFEST.v7.1.0.json"
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": payload["status"], "chr22_L": len(variants), "primary_hgdp_n": len(groups["hgdp_primary_holdout"]), "output": str(output)}, indent=2))


if __name__ == "__main__":
    main()
