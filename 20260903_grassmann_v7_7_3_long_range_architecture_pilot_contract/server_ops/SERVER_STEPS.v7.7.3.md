# v7.7.3 server steps

Deploy only after v7.7.2 result and its manifest both verify. Copy validated R35
to R36, overlay this patch, run validator, manifest verification, and unit tests.
Then build the readiness artifacts by supplying the immutable
`TASK_VALIDITY_EXECUTION.v7.7.2.json`. This stage performs no GPU inspection,
allocation, or execution.

The valid terminal status is
`LONG_RANGE_ARCHITECTURE_PILOT_CONTRACT_SIGNED_IMPLEMENTATION_ONLY`; it permits
only a future implementation-only v7.7.4 package.
