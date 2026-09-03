# Grassmann v7.7.x 执行标准与代码格式对齐规范

以权威 v7.7.3 合约包（SHA-256 `b4a4aaec761a33e73882b48362860f53628b25913ab90e4c5ef4017fdd6328d1`）
为基准，固化后续我方交付（v7.7.4 起）必须遵循的执行标准、代码格式与命名细节。
本环境已复核：`validate_v7_7_3.py` → PASS；`sha256sum -c MANIFEST.v7.7.3.sha256` → 全 OK；
`unittest discover -s p0/tests -p 'test_v7_7_3.py'` → 7 tests OK；builder 端到端 →
`LONG_RANGE_ARCHITECTURE_PILOT_CONTRACT_SIGNED_IMPLEMENTATION_ONLY`，schedule 48 行。

---

## 1. 目录、包与部署

- 补丁目录名：`YYYYMMDD_grassmann_v7_7_N_<slug>`。
  v7.7.3 的 slug 带 `_contract` 后缀：`20260903_grassmann_v7_7_3_long_range_architecture_pilot_contract`。
- 上传包：`deliverables/v7_7_N_<slug>_patch_YYYYMMDD.tar.gz`，**顶层即该 dated 目录**，
  部署用 `tar -xzf "$PATCH" --strip-components=1 -C "$R36"`（strip 掉顶层目录，覆盖到 R36 根）。
- release 递进：`cp -a "$R35" "$R36"` 后覆盖补丁；若 `$R36` 已存在则 STOP，不覆盖。
- Python：`$V7_PY = .../releases/v7_20260825_p0_r1/.conda-v7/bin/python`。
  **合约/实现阶段代码必须纯 stdlib**（argparse/hashlib/json/datetime/pathlib）；
  numpy/torch 只在真正执行阶段（如 v7.7.2）出现。

## 2. 补丁文件清单（8 个文件；MANIFEST 收录其余 7 个，不含自身）

```
PROTOCOL_ADDENDUM.v7.7.N.md
<STAGE>_SPEC.v7.7.N.yaml
DECISIONS.v7.7.N.tsv
validate_v7_7_N.py
p0/<builder>.py
p0/tests/test_v7_7_N.py
server_ops/SERVER_STEPS.v7.7.N.md
MANIFEST.v7.7.N.sha256
```

- **不单独放 ARM_MAP / schedule TSV 到补丁里**。2×2 阶乘写在 SPEC yaml 内联；
  pilot schedule（`LONG_RANGE_PILOT_DRAFT.v7.7.3.tsv`）由 builder 写到**输出目录**，不进补丁。

## 3. MANIFEST 与 readiness 清单的路径约定（关键，易踩坑）

- **补丁 `MANIFEST.v7.7.N.sha256`**：每行 `<sha256><两个空格><相对路径>`，
  相对路径**不带 `./` 前缀**（如 `p0/build_architecture_pilot_readiness_v7_7_3.py`），
  按**大小写不敏感字典序**排列；共 7 行（不含 MANIFEST 自身）。
- **builder 输出的 readiness 清单**（如 `LONG_RANGE_ARCHITECTURE_PILOT_READINESS_MANIFEST.v7.7.3.sha256`）：
  每行**带 `./` 前缀**（`./LONG_RANGE_..._READINESS.v7.7.3.json`）。
- 因为是相对路径，`sha256sum -c` **必须先 `cd` 进对应目录**再校验：
  ```
  ( cd "$RUN_DIR" && sha256sum -c TASK_VALIDITY_EXECUTION_MANIFEST.v7.7.2.sha256 )
  ( cd "$OUT"     && sha256sum -c LONG_RANGE_ARCHITECTURE_PILOT_READINESS_MANIFEST.v7.7.3.sha256 )
  ```
  （上一轮我在根目录直接 `-c` 报错，就是没进目录——非补丁问题。）

## 4. validate 脚本风格

- `ROOT = Path(__file__).resolve().parent`；失败用 `def fail(msg): raise SystemExit(f"FAIL: {msg}")`。
- 断言方式 = **在 SPEC / PROTOCOL 文本里查字面子串** + 校验 `len(manifest) == 7`。
- **不在 validate 内做 manifest 的 sha256 自校验**（那一步由 server 用 `sha256sum -c` 单独跑）。
- 成功打印 `json.dumps({..., "status":"PASS", ...}, indent=2, sort_keys=True)`。
- 子串需与文档**字面一致**，注意换行：SPEC 内用无空格形式
  `NLL(R1_ROUTER)-NLL(R3_ROUTER_GRASSMANN)`；PROTOCOL 正文用带空格形式
  `NLL(R1_ROUTER) - NLL(R3_ROUTER_GRASSMANN)`；validate/test 各查各自的形式。

## 5. 单元测试风格

- 用 `importlib.util.spec_from_file_location` 的 `load()` 按路径加载 builder；
  `ROOT = Path(__file__).resolve().parents[2]`。
- 运行命令固定：`"$V7_PY" -m unittest discover -s p0/tests -p 'test_v7_7_N.py' -v`。
- 测试内容 = 查 SPEC/PROTOCOL 子串 + 校验 builder 纯函数（如 `pilot_rows()` 长度/集合）。
- **测试数量要与 server steps 里声明的一致**（v7.7.3 = 7 个）。

## 6. builder（readiness）契约

- CLI：`--task-validity <TASK_VALIDITY_EXECUTION.v7.7.2.json>`（**不是** `--source-execution`）
  + `--output-dir <OUT>`。
- 先 `if output_dir.exists(): raise FileExistsError("refuse overwrite: ...")`，再读源。
- 源校验（全部满足才继续）：
  `status == "LONG_RANGE_TASK_VALIDITY_PASS"`；
  `authorization.gpu_used is not False → raise`；
  `gates.all_pass is not True → raise`。
- 产出到 OUT：readiness JSON + schedule TSV + readiness 清单（带 `./`）。
- 所有写文件：`encoding="utf-8", newline="\n"`；JSON：`indent=2, sort_keys=True` + 末尾 `"\n"`。
- 文件头一律 `from __future__ import annotations`。

## 7. 冻结的状态字符串 / 命名（v7.7.3 权威口径）

| 项 | 值 |
| --- | --- |
| stage | `LONG_RANGE_ARCHITECTURE_PILOT_CONTRACT` |
| 源要求状态 | `LONG_RANGE_TASK_VALIDITY_PASS` |
| 终态 | `LONG_RANGE_ARCHITECTURE_PILOT_CONTRACT_SIGNED_IMPLEMENTATION_ONLY` |
| 失败态 | `LONG_RANGE_ARCHITECTURE_PILOT_NOT_READY` |
| 唯一下一步 | `IMPLEMENT_V7_7_4_LONG_RANGE_PILOT_HARNESS_NO_LAUNCH` |
| 四格 cell id | `R0_LOCAL` / `R1_ROUTER` / `R2_GRASSMANN_ONLY` / `R3_ROUTER_GRASSMANN` |
| 主对比 | `NLL(R1_ROUTER)-NLL(R3_ROUTER_GRASSMANN)`（正值利好 Grassmann） |
| 实践边界 | `0.010` nats/target |
| 未来 GO 规则 | 配对单侧 95% LCB > 0.010 且 validity+fairness+shuffled 控制全过 |
| pilot 真值 seed | `77301–77306`（6 个，**与 v7.7.2 的 77201–77205 互斥**） |
| 嵌套 init seed | `87401, 87402`（组内平均，不算独立复制） |
| effective_n / cells | `6` / `48`（6×2×4，全部 execution_authorized=false） |

## 8. DECISIONS.tsv schema（v7.7.3 起演进为 5 列）

- 列：`decision_id  state  decision  constraint  status`。
- id 用零填充 `V77N-0KK`（如 `V773-001`）；`state=FROZEN`；`status=PASS`。
- 注意：v7.7.0–v7.7.2 是旧 3 列（`decision_id decision status`，status=FROZEN）；
  **v7.7.3 起用新 5 列**，后续沿用最新权威格式。

## 9. Grassmann 臂公平性设计约束（写进合约，供 v7.7.4 实现）

- router-present 两格必须共享：source positions、token embedding、router topology、
  training split、target mask、optimizer family、LR schedule、batch construction、eval code。
- Grassmann 分支 = **在共享 router 上的加性残差**（additive residual）；
  Grassmann 臂**不得**单独引入额外原始输入 / oracle 特征 / 不同 target / 不同验证集。
- 盲化防火墙：pilot 只可释放 `paired_dispersion / parameter_counts / compute_profiles /
  control_integrity`；**arm means / rankings / p-values / GO-NO-GO 一律不释放**；
  正式样本量只能由预声明的 0.010 边界 + 独立真值 seed 精度计算得出，**不得用观测到的 pilot 均值效应选 n**。

---

## 10. 我上一版草案 → 权威版：需修正的偏差（已据此校正仓库）

| 项 | 我草案（偏差） | 权威（应对齐） |
| --- | --- | --- |
| 目录 slug | `..._architecture_pilot` | `..._architecture_pilot_contract` |
| builder CLI | `--source-execution` | `--task-validity` |
| 终态 | `..._SIGNED_NO_GPU` | `..._SIGNED_IMPLEMENTATION_ONLY` |
| 下一步 | `IMPLEMENT_V7_7_4_CPU_BLINDED_VARIANCE_PILOT_NO_GPU` | `IMPLEMENT_V7_7_4_LONG_RANGE_PILOT_HARNESS_NO_LAUNCH` |
| cell id | `LR_A00/A10/A01/A11` | `R0_LOCAL/R1_ROUTER/R2_GRASSMANN_ONLY/R3_ROUTER_GRASSMANN` |
| 主对比 | `NLL(LR_A01)-NLL(LR_A11)` | `NLL(R1_ROUTER)-NLL(R3_ROUTER_GRASSMANN)` |
| pilot seed | 未定 seed，谈 variance pilot | 固定 `77301–77306`，与 77201–77205 互斥 |
| DECISIONS | 旧 3 列 / `V773-D01` / status=FROZEN | 新 5 列 / `V773-001` / state=FROZEN,status=PASS |
| 补丁附件 | 额外放了 ARM_MAP.tsv | 不放；阶乘内联 SPEC，schedule 由 builder 出到 OUT |
| MANIFEST 路径 | 带 `./` | 补丁 MANIFEST 不带 `./`（readiness 清单才带 `./`） |
| validate | 自带 manifest sha 自校验 | 不自校验；server 用 `sha256sum -c` 单独跑 |

## 11. 阶段定位（供对齐语境，非授权）

现处于 v7.7 “长程任务有效性与实验设计”阶段，**尚未进入 Grassmann 长程架构测试**：
- v7.7.0 冻结规则；v7.7.1 验证生成器无 local 泄漏；
- v7.7.2 实测确认 local 无能、普通 global 可解、打乱标签不可学（四 gate PASS）；
- v7.7.3 冻结 Grassmann vs 同一 global baseline 的 2×2、独立重复、盲化与公平性；
- **v7.7.4** 才实现真正的 pilot harness（no launch）；之后经审计与单独授权，才启动真正的长程 A1 测试与 GPU。

这套顺序是为了不重演 A1-R：任务本身若不需要长程信息、或 baseline 不公平，
即便跑出数字也不能证明 Grassmann 的长程价值。
