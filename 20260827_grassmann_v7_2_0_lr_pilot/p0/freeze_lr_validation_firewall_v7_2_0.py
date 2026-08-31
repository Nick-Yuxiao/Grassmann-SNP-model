from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


def stable_int(*parts: object) -> int:
    value = "\t".join(map(str, parts)).encode("utf-8")
    return int.from_bytes(hashlib.sha256(value).digest()[:8], "big")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--validation-tsv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise SystemExit(f"refusing to overwrite: {args.output_dir}")

    ids = [line for line in (args.data_dir / "VALIDATION_SAMPLES.txt").read_text(encoding="utf-8").splitlines() if line]
    if len(ids) != 249 or len(ids) != len(set(ids)):
        raise ValueError("validation sample IDs must contain 249 unique values")
    with args.validation_tsv.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    metadata = {row["sample_id"]: row for row in rows}
    if set(ids) != set(metadata):
        raise ValueError("validation TSV sample set differs from preprocessed validation IDs")

    seen: set[int] = set()
    historical_by_seed: dict[str, list[str]] = {}
    for seed in (91001, 91002):
        ranked = sorted(range(len(ids)), key=lambda index: stable_int("validation_sample", seed, ids[index]))
        selected = ranked[:32]
        seen.update(selected)
        historical_by_seed[str(seed)] = [ids[index] for index in selected]
    decision = sorted(set(range(len(ids))) - seen)
    tuning = sorted(seen)
    if set(tuning) & set(decision) or len(tuning) + len(decision) != 249:
        raise ValueError("validation firewall partition failure")

    tuning_pop = Counter(metadata[ids[index]]["population"] for index in tuning)
    decision_pop = Counter(metadata[ids[index]]["population"] for index in decision)
    all_pop = Counter(metadata[sample_id]["population"] for sample_id in ids)
    absent = sorted(pop for pop in all_pop if decision_pop[pop] == 0)
    if absent:
        raise ValueError(f"untouched decision holdout loses populations: {absent}")

    args.output_dir.mkdir(parents=True, exist_ok=False)
    outputs = {
        "LR_TUNING.indices.txt": "".join(f"{index}\n" for index in tuning),
        "LR_TUNING.samples.txt": "".join(f"{ids[index]}\n" for index in tuning),
        "DECISION_HOLDOUT.indices.txt": "".join(f"{index}\n" for index in decision),
        "DECISION_HOLDOUT.samples.txt": "".join(f"{ids[index]}\n" for index in decision),
    }
    for name, content in outputs.items():
        (args.output_dir / name).write_text(content, encoding="utf-8")
    audit = {
        "schema_version": "1.0",
        "protocol_version": "v7.2.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PASS",
        "historical_rule": "FIRST_32_BY_STABLE_INT_FOR_EACH_MASK_SEED",
        "historical_mask_seeds": [91001, 91002],
        "historical_by_seed": historical_by_seed,
        "tuning_count": len(tuning),
        "decision_holdout_count": len(decision),
        "overlap": 0,
        "decision_holdout_was_historically_viewed": False,
        "population_counts": {
            "all_validation": dict(sorted(all_pop.items())),
            "lr_tuning": dict(sorted(tuning_pop.items())),
            "decision_holdout": dict(sorted(decision_pop.items())),
        },
        "decision_populations_absent": absent,
        "inputs": {
            "validation_samples_sha256": sha256(args.data_dir / "VALIDATION_SAMPLES.txt"),
            "validation_tsv_sha256": sha256(args.validation_tsv),
        },
        "hgdp_used": False,
    }
    audit_path = args.output_dir / "VALIDATION_FIREWALL.v7.2.0.json"
    audit_path.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    manifest_inputs = sorted(path for path in args.output_dir.iterdir() if path.is_file())
    (args.output_dir / "VALIDATION_FIREWALL.v7.2.0.sha256").write_text(
        "".join(f"{sha256(path)}  {path.name}\n" for path in manifest_inputs), encoding="utf-8"
    )
    print(json.dumps({
        "status": "PASS", "tuning_count": len(tuning),
        "decision_holdout_count": len(decision), "decision_populations": len(decision_pop),
        "output_dir": str(args.output_dir),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

