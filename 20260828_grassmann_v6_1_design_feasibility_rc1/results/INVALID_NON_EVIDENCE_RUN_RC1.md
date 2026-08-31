# Invalid metadata count in first non-evidence run

_Applies only to `results/non_evidence_run_rc1`; recorded 2026-08-28 without deleting or overwriting that directory._

---

## ⚠️ Error

`SUMMARY.json` labels 2,880 angle-specific fits as `independent_families`. Angle cells intentionally reuse paired base seeds, so this label is incorrect. The correct accounting is 2,880 angle-specific fits, 576 unique paired base families, and 12 independent family replicates within each cell.

## 📋 Consequence

The numerical rows are unchanged, but the first run must not be used even for informal reporting because its replication metadata is misleading. The corrected script writes new output to an immutable `attempt2` directory.

