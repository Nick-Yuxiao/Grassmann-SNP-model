# Server operations: separate v6/v7 records

Server root is fixed to `/data1/home/tanyuxiao/Grassmann_model`.

The bootstrap creates, without deleting or moving existing content:

```text
Grassmann_model/
├── incoming/
├── v6/
│   ├── code/releases/
│   ├── results/
│   ├── resources/{audits,gpu_test}/
│   ├── logs/
│   ├── locks/
│   └── runbooks/
└── v7/
    ├── code/releases/
    ├── results/
    ├── resources/{audits,gpu_test}/
    ├── logs/
    ├── locks/
    └── runbooks/
```

Every result/resource directory is timestamped and newly created. Scripts do
not use `rm`, `kill`, `pkill`, scheduler cancellation or preemption.

## Server-side commands after extraction

```bash
cd /data1/home/tanyuxiao/Grassmann_model/v7/code/releases/<release_id>
bash server_ops/bootstrap_v6_v7.sh
```

Inspect the generated audit under `v7/resources/audits/`. It records all GPU
compute processes and the current user's process/Slurm state.

Create the isolated cu128 environment only when CPU/disk/network use is
acceptable:

```bash
bash p0/run_t00_nonintrusive.sh
```

Run the short GPU test with the isolated interpreter:

```bash
bash server_ops/run_gpu_test_nonintrusive.sh \
  --python /data1/home/tanyuxiao/Grassmann_model/v7/code/releases/<release_id>/.venv/bin/python
```

If every GPU has a compute process, more than 1024 MiB used memory, more than
5% utilization, or the v7 project lock is held, the test exits without starting.
It rechecks occupancy after taking the lock. A PASS is written under both
`v7/results/gpu_test/<run_id>/` and `v7/resources/gpu_test/<run_id>/`.

Only after this short test passes should T03 be started:

```bash
bash p0/run_t03_nonintrusive.sh \
  /data1/home/tanyuxiao/Grassmann_model/v7/results/t03_profile_<UTC_timestamp>
```
