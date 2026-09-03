# Server steps: v7.7.3

Deploy this patch over validated R35 into new immutable R36. This release is
CPU-only and a contract draft: it must not inspect, reserve, allocate, or use a GPU,
must not train Grassmann, and must not fit any factorial arm.

Validate the release manifest with `validate_v7_7_3.py` and require exactly eight unit
tests in `p0/tests/test_v7_7_3.py`. Then build the architecture-pilot readiness record
from the completed v7.7.2 execution JSON:

```
python3 p0/build_architecture_pilot_readiness_v7_7_3.py \
  --source-execution <path to TASK_VALIDITY_EXECUTION.v7.7.2.json> \
  --output-dir v7/results/architecture_pilot/v7.7.3/<UTC>_architecture_pilot_v7_7_3
```

The builder accepts the source only when it is a v7.7.2 record with status
`LONG_RANGE_TASK_VALIDITY_PASS`, all four task-validity gates passed, `gpu_used`
false, and `grassmann_fitted` false. It runs only structural checks on the frozen 2x2
arm map; it fits nothing.

The expected terminal status is `LONG_RANGE_ARCHITECTURE_PILOT_CONTRACT_SIGNED_NO_GPU`.
PASS of v7.7.2 confirmed the synthetic task discriminates local-insufficient from
global-solvable; it is not a Grassmann success and not a biological claim. The only
next stage is `IMPLEMENT_V7_7_4_CPU_BLINDED_VARIANCE_PILOT_NO_GPU`, which will estimate
the primary-contrast dispersion on a reduced CPU proxy with arm means withheld and
freeze the formal truth-seed count. No GPU factorial, A3, HAPNEST, HGDP, phenotype, or
biological long-range work is authorized here.
