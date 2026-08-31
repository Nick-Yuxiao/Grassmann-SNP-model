# V7 P0 readiness verdict

- Status: **NOT_READY**
- Scope: readiness to enter A0 only; not the A1 scientific verdict.
- Blocking tasks: T00, T01, T03, T04

- T00: BLOCKED — Environment smoke, server resource record, non-interference audit and lockfile are required. A busy audit is valid evidence but no new GPU job may start until a fresh idle audit passes.
- T01: BLOCKED — A signed branch decision and hashed panel inputs are required.
- T02: PASS — Every independently named frozen constant must appear in both the decision ledger and metric dictionary.
- T03: BLOCKED — A 100-step CUDA profile at all exact P0/A1 planning lengths is required; CPU dry-runs are never accepted.
- T04: BLOCKED — The signed capacity must cover the measured A1 plan at <=80%, including the 2x engineering margin.
