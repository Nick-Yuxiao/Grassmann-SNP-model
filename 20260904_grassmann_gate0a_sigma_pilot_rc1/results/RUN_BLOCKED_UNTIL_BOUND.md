# RUN BLOCKED until the data contract is BOUND

_No sigma-pilot result exists yet. The run is mechanically refused until `BINDING.json` (in the rc2 package) reads `status == BOUND`._

---

## Order (do not reverse)

1. Bind panel manifest hash → 2. bind LD-block version hash → 3. freeze/verify preprocessing (`BINDING.json` → `BOUND`) → 4. export the frozen panel to the `.npy` interface → 5. run this pilot on the server → 6. read `sigma_hat` per cell → 7. map to required `R` → 8. freeze `R_formal = max` over primary-relevant cells.

See `server_ops/SERVER_STEPS.sigma_pilot.md`. The run script (`scripts/run_sigma_pilot.py`) checks the binding status first and prints `RUN_BLOCKED` if the contract is not bound, before touching any data.

## What it will produce

`results/sigma_pilot_result.json`: per `(regime, k)` cell `sigma_hat` and required `R`, plus `R_formal`. Variance calibration only — it decides `R` and compute feasibility, nothing else.
