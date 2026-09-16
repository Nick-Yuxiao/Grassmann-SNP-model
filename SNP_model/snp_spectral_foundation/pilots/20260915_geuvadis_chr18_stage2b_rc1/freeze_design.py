from __future__ import annotations

import csv
import hashlib
import json
import random
from collections import Counter, defaultdict
from pathlib import Path


SEED = 5_262_026
ROLE_TARGETS = {"development": 150, "task_gate": 50, "bridge_test": 230}
DEV_TARGETS = {"dev_train": 120, "dev_validation": 30}

PILOT_DIR = Path(__file__).resolve().parent
ROOT = PILOT_DIR.parents[2]
PACKAGE = ROOT / "snp_spectral_foundation"
ASSET = PACKAGE / "pilot_assets" / "geuvadis_tensorqtl_chr18"
PSAM = ASSET / "GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.chr18.psam"
PANEL = ASSET / "integrated_call_samples_v3.20130502.ALL.panel"
PED = ASSET / "integrated_call_samples_v3.20250704.ALL.ped"
RC2 = PACKAGE / "pilots" / "20260915_geuvadis_chr19_development_rc1" / "results" / "RESULTS.json"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def allocate(groups: dict[str, list[str]], metadata: dict[str, dict[str, str]], targets: dict[str, int], seed: int) -> dict[str, str]:
    roles = list(targets)
    total_target = sum(targets.values())
    if sum(len(v) for v in groups.values()) != total_target:
        raise ValueError("allocation targets do not equal eligible sample count")
    desired_by_pop: dict[str, dict[str, float]] = {}
    pop_counts = Counter(metadata[samples[0]]["pop"] for samples in groups.values())
    for pop, count in pop_counts.items():
        desired_by_pop[pop] = {role: count * targets[role] / total_target for role in roles}
    current_total = Counter()
    current_pop: dict[str, Counter[str]] = defaultdict(Counter)
    rng = random.Random(seed)
    ordered = list(groups.items())
    rng.shuffle(ordered)
    ordered.sort(key=lambda item: (-len(item[1]), metadata[item[1][0]]["pop"], item[0]))
    assignment: dict[str, str] = {}
    for group, samples in ordered:
        pop = metadata[samples[0]]["pop"]
        size = len(samples)
        candidates = []
        for role in roles:
            if current_total[role] + size > targets[role]:
                continue
            before_pop = current_pop[pop][role] - desired_by_pop[pop][role]
            after_pop = current_pop[pop][role] + size - desired_by_pop[pop][role]
            before_total = current_total[role] - targets[role]
            after_total = current_total[role] + size - targets[role]
            cost = (after_pop * after_pop - before_pop * before_pop) + 0.05 * (
                after_total * after_total - before_total * before_total
            )
            candidates.append((cost, current_total[role] / targets[role], roles.index(role), role))
        if not candidates:
            raise RuntimeError(f"cannot place family {group} of size {size} without exceeding exact targets")
        role = min(candidates)[-1]
        current_total[role] += size
        current_pop[pop][role] += size
        for sample in samples:
            assignment[sample] = role
    if dict(current_total) != targets:
        raise AssertionError((dict(current_total), targets))
    return assignment


def main() -> None:
    psam = read_rows(PSAM)
    panel = {row["sample"]: row for row in read_rows(PANEL)}
    ped = {row["Individual ID"]: row for row in read_rows(PED)}
    sample_ids = [row["#IID"] for row in psam]
    if len(sample_ids) != 445 or len(set(sample_ids)) != 445:
        raise ValueError("unexpected PSAM sample count or duplicates")
    if any(sample not in panel or sample not in ped for sample in sample_ids):
        raise ValueError("not all GEUVADIS samples map to official panel and pedigree metadata")
    rc2 = json.loads(RC2.read_text(encoding="utf-8"))
    excluded = set(rc2["split"]["test_sample_ids"])
    overlap = excluded.intersection(sample_ids)
    if overlap != excluded or len(overlap) != 15:
        raise ValueError("RC2 test exclusion set does not map exactly to the 445-sample cohort")
    eligible = [sample for sample in sample_ids if sample not in excluded]
    if len(eligible) != 430:
        raise AssertionError("expected 430 eligible samples")

    metadata: dict[str, dict[str, str]] = {}
    groups: dict[str, list[str]] = defaultdict(list)
    for sample in eligible:
        p = panel[sample]
        d = ped[sample]
        if p["pop"] != d["Population"]:
            raise ValueError(f"population mismatch for {sample}")
        sex_psam = next(row["SEX"] for row in psam if row["#IID"] == sample)
        expected_sex = "1" if p["gender"] == "male" else "2"
        if sex_psam != expected_sex or d["Gender"] != expected_sex:
            raise ValueError(f"sex mismatch for {sample}")
        family = d["Family ID"] or sample
        metadata[sample] = {"pop": p["pop"], "super_pop": p["super_pop"], "sex": sex_psam, "family": family}
        groups[family].append(sample)
    for family_samples in groups.values():
        if len({metadata[s]["pop"] for s in family_samples}) != 1:
            raise ValueError("family spans populations")

    role = allocate(groups, metadata, ROLE_TARGETS, SEED)
    dev_groups = {family: samples for family, samples in groups.items() if role[samples[0]] == "development"}
    dev_role = allocate(dev_groups, metadata, DEV_TARGETS, SEED + 1)
    rows = []
    for sample in sample_ids:
        if sample in excluded:
            rows.append({"sample_id": sample, **metadata.get(sample, {"pop": panel[sample]["pop"], "super_pop": panel[sample]["super_pop"], "sex": next(x["SEX"] for x in psam if x["#IID"] == sample), "family": ped[sample]["Family ID"] or sample}), "primary_role": "EXCLUDED_RC2_TEST", "development_role": ""})
        else:
            rows.append({"sample_id": sample, **metadata[sample], "primary_role": role[sample], "development_role": dev_role.get(sample, "")})

    PILOT_DIR.mkdir(parents=True, exist_ok=True)
    manifest = PILOT_DIR / "SAMPLE_MANIFEST.tsv"
    fields = ["sample_id", "family", "pop", "super_pop", "sex", "primary_role", "development_role"]
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "design_id": "0526-stage2b-geuvadis-chr18-rc1",
        "seed": SEED,
        "sample_counts": dict(Counter(row["primary_role"] for row in rows)),
        "development_counts": dict(Counter(row["development_role"] for row in rows if row["development_role"])),
        "population_by_role": {
            pop: dict(Counter(row["primary_role"] for row in rows if row["pop"] == pop))
            for pop in sorted({row["pop"] for row in rows})
        },
        "family_counts_by_role": {
            role_name: len({row["family"] for row in rows if row["primary_role"] == role_name})
            for role_name in ["development", "task_gate", "bridge_test", "EXCLUDED_RC2_TEST"]
        },
        "rc2_test_ids_excluded": sorted(excluded),
        "family_overlap_checks": {
            "development_task": False,
            "development_bridge": False,
            "task_bridge": False,
        },
        "source_sha256": {path.name: sha256_file(path) for path in [PSAM, PANEL, PED, RC2]},
        "manifest_sha256": sha256_file(manifest),
    }
    family_sets = {
        role_name: {row["family"] for row in rows if row["primary_role"] == role_name}
        for role_name in ["development", "task_gate", "bridge_test"]
    }
    summary["family_overlap_checks"] = {
        "development_task": bool(family_sets["development"] & family_sets["task_gate"]),
        "development_bridge": bool(family_sets["development"] & family_sets["bridge_test"]),
        "task_bridge": bool(family_sets["task_gate"] & family_sets["bridge_test"]),
    }
    if any(summary["family_overlap_checks"].values()):
        raise AssertionError("family leakage")
    (PILOT_DIR / "DESIGN_FREEZE.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
