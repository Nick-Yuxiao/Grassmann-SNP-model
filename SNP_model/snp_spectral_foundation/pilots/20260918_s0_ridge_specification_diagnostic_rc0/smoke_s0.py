"""Synthetic end-to-end rehearsal for S0. No assets, no scientific meaning.

The synthetic outcome is built so the answer is known in advance: covariates
carry real signal and the cis dosages carry some. A faithful S0 should then show
S1 losing ground as the dosage block grows while S2 and S3 stay flat.

24 synthetic genes are too few for a bootstrap interval to clear zero, so the
rehearsal ends on NOT-SUPPORTED by construction. That a preserved covariate fit
still finds real dosage signal is asserted in test_s0, not here.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

PILOT_DIR = Path(__file__).resolve().parent
if str(PILOT_DIR) not in sys.path:
    sys.path.insert(0, str(PILOT_DIR))

import run_s0  # noqa: E402
from s0_core import SPECIFICATIONS  # noqa: E402

POPULATIONS = ["CEU", "FIN", "GBR", "TSI", "YRI"]
RESULT_DIR = PILOT_DIR / "results_smoke"


def build(n_people=120, n_variants=9000, n_genes=24, seed=404):
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n_people):
        rows.append(
            {
                "sample_id": f"S{i:04d}",
                "family": f"FAM{i // 2:04d}",
                "sex": "2" if i % 2 else "1",
                "pop": POPULATIONS[i % len(POPULATIONS)],
                "development_role": "dev_train" if i < 96 else "dev_validation",
                "primary_role": "development",
            }
        )
    positions = np.sort(rng.choice(np.arange(1, 40_000_000), size=n_variants, replace=False))
    frequency = rng.uniform(0.1, 0.9, size=n_variants)
    genotype = rng.binomial(2, frequency, size=(n_people, n_variants)).astype(np.int8)

    dev_idx = np.arange(n_people)
    sex = np.asarray([1.0 if r["sex"] == "2" else 0.0 for r in rows])
    pop = np.asarray([POPULATIONS.index(r["pop"]) for r in rows], dtype=np.float64)
    tss_points = np.linspace(2_000_000, 38_000_000, n_genes).astype(int)

    records = []
    for g, tss in enumerate(tss_points):
        nearby = np.flatnonzero(np.abs(positions - tss) <= 1_000_000)
        causal = rng.choice(nearby, size=min(3, len(nearby)), replace=False)
        cis = genotype[:, causal].astype(np.float64) @ rng.normal(size=len(causal))
        # Covariates carry most of the signal and the cis dosages carry some, so a
        # specification that preserves the covariate fit should stay flat as the
        # dosage block grows while a shared penalty should visibly lose ground.
        values = (
            2.0 * sex
            + 0.8 * pop
            + 1.6 * cis
            + rng.normal(scale=1.0, size=n_people)
        )
        records.append(
            {"gene_id": f"ENSGSMOKE{g:05d}", "tss_1based": int(tss), "values": values}
        )

    config = json.loads((PILOT_DIR / "CONFIG_S0.json").read_text(encoding="utf-8"))
    config.update(
        {
            "max_genes": n_genes,
            "min_genes_required": 4,
            "snps_per_gene": 64,
            "snp_counts": [8, 32, 64],
            "pc_panel_snps": 2000,
            "cv_folds": 5,
            "bootstrap_replicates": 200,
            "ridge_alphas": [1e-4, 1e-2, 1.0, 10.0, 100.0, 1e4, 1e6],
        }
    )
    return config, rows, positions, genotype, records, dev_idx


def main() -> None:
    config, rows, positions, genotype, records, dev_idx = build()
    families = [r["family"] for r in rows]
    print(
        "synthetic chromosome: %d people x %d variants, %d genes"
        % (len(rows), len(positions), len(records))
    )

    def load_records(idx):
        if not np.array_equal(np.asarray(idx), dev_idx):
            raise SystemExit("expression may only be read for the development rows")
        return records

    results = run_s0.analyse(config, rows, positions, genotype, load_records, dev_idx, families)
    status = run_s0.verdict(results, {**config, "snps_per_gene": config["snp_counts"][-1]})

    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    (RESULT_DIR / "RESULTS_S0.json").write_text(
        json.dumps(results, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    (RESULT_DIR / "FINAL_STATUS_S0.json").write_text(
        json.dumps(status, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    print("\n--- smoke summary (synthetic; no scientific meaning) ---")
    print("%-6s %6s %10s %10s %11s" % ("spec", "snps", "R2(A)", "R2(B)", "B-A"))
    for specification in SPECIFICATIONS:
        for count in config["snp_counts"]:
            entry = results["contrasts"][f"B-A@{specification}:{count}"]
            print(
                "%-6s %6d %+10.4f %+10.4f %+11.4f"
                % (
                    specification,
                    count,
                    results["macro_r2"][f"A@{specification}"],
                    results["macro_r2"][f"B@{specification}:{count}"],
                    entry["macro_delta_r2"],
                )
            )
    print("\nstatus:", status["status"])
    print("artifacts written to  :", RESULT_DIR)


if __name__ == "__main__":
    main()
