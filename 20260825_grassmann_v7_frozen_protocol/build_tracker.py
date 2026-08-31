from __future__ import annotations

import argparse
import csv
from pathlib import Path


EXPECTED_IDS = [f"T{i:02d}" for i in range(16)]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    with args.source.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if [row["task_id"] for row in rows] != EXPECTED_IDS:
        raise ValueError("Unexpected task IDs or order")

    by_id = {row["task_id"]: row for row in rows}
    for row in rows:
        row["status"] = "NOT_STARTED"

    t = by_id["T00"]
    t.update(
        start_date="2026-08-25",
        action="在已授权 RTX 5090 服务器建立同时满足 Python>=3.10 / Torch>=2.7 / CUDA 12.8+ / sm_120 的环境",
        execution_standard="锁定版本；核验 GPU 型号、compute capability、驱动与显存；CUDA 张量分配、前向、反向传播均通过；环境导出为可重建 lockfile",
        deliverable="environment/requirements-cu128.lock; environment/ENV_SMOKE.json; environment/SERVER_RESOURCE.json",
        acceptance_standard="RTX 5090、compute capability>=12.0、Torch>=2.7、CUDA 可用且反向传播有限；lockfile 有 hash 并可一条命令重建",
        status="IN_PROGRESS_REMOTE_ACCESS_PENDING",
    )

    t = by_id["T02"]
    t.update(
        start_date="2026-09-07",
        end_date="2026-09-08",
        action="登记并复核已冻结判据常数（各自独立命名，禁止全表字符串替换）",
        execution_standard="delta_min=0.010; delta_NI=0.010; delta_LD=0.010; overfit_thr=0.010; pc_control_thr=0.005；Delta=NLL_comparator-NLL_candidate，正数=candidate 更好；NI 检验 UCB(NLL_candidate-NLL_comparator)<delta_NI",
        deliverable="DECISIONS.v7.0.1.tsv; METRIC_DEFINITIONS.md",
        acceptance_standard="每个常数有独立名字、数值、单位和适用行；数值不得根据 A1/A2/A3 结果调整",
    )

    t = by_id["T03"]
    t.update(
        start_date="2026-09-07",
        end_date="2026-09-09",
        action="在实际 RTX 5090 上做 profiler：单 cell 显存峰值、吞吐、端到端时长和磁盘占用",
        execution_standard="L=8192 与 32GB 显存允许的最大 L 各做一次 100-step profile；包含 CUDA synchronize、数据加载和 checkpoint；记录 OOM 边界；不得用旧 harness 或理论 FLOPs 外推",
    )

    t = by_id["T04"]
    t.update(
        execution_standard="预算单位为 T03 实测 GPU-小时；5 data_seed x 2 init_seed；含 2 倍工程余量；计划需求不超过已签署 5090 GPU-小时和存储容量的 80%",
        acceptance_standard="A1 全量在已签容量内；GPU 数量、可用时段、并发上限和存储配额明确；A2/A3 未签署时标记 UNFUNDED",
    )

    t = by_id["T05"]
    t.update(
        action="实现 harness：数据管线、五类 mask、三个模型、双公平口径、断点续跑和 GRID_SPEC 展开器",
        execution_standard="run 由 GRID_SPEC.v7.0.1.yaml 展开；run_id 重复、占位符残留或依赖缺失均报错；禁止 silent retry；resolver 只生成新的只读 resolved TSV",
        deliverable="src/; GRID_SPEC.v7.0.1.yaml; RESOLVER.md",
    )

    t = by_id["T07"]
    t.update(
        action="chr22 headroom：random/contiguous 诊断；ld_block/within_long 确认；local/local+PC/Grassmann",
        execution_standard="random/contiguous 只做 matched-parameter；ld_block/within_long 做 matched-parameter 与 matched-compute；5 data_seed x 2 init_seed；先按 init 求均值，再按 data seed 配对推断",
        acceptance_standard="所有计划 run 均有终态记录；失败 run 计入分母并列因；不得删除不利 seed；确认性 mask 两种公平口径完整",
    )

    t = by_id["T08"]
    t.update(
        action="建立参考线：四个 within-chrom mask 的 closed-form LD reference，以及 cross_chrom 的 ancestry-only reference",
        execution_standard="ancestry-only reference 不是 Bayes 上界，仅判断祖源信息后是否仍有可测增量；亲缘与远端 admixture LD 不在其覆盖内",
        acceptance_standard="closed-form LD 在 random/contiguous/ld_block/within_long 上可得；ancestry-only reference 在 cross_chrom 上可得并附局限说明",
    )

    t = by_id["T10"]
    t.update(
        execution_standard="单 mask：两种必需公平口径均 LCB>delta_min 才 GO；任一必需口径 UCB<=delta_min 则该 mask NO-GO；其余 INCONCLUSIVE。scope：任一确认性 mask GO 即 GO；全部 NO-GO 才 NO-GO",
        acceptance_standard="within-chrom 与 cross-chrom 各有三态结论；INCONCLUSIVE 最多追加一个预定 seed tranche，之后保守记为 NO-GO；完整记录量词和 CI",
    )

    t = by_id["T11"]
    t.update(
        action="resolver 对 T10 判定为 GO 的每个 mask family 展开六算子 x 双公平口径",
        execution_standard="六者仅混合原语不同；特征维、参数量和优化预算一致；每个 GO family 独立展开；resolved 依赖只包含实际生成行",
        failure_action="无 GO family 时写 A2_VERDICT.json=NOT_APPLICABLE；无签署 GPU/存储合同时保持 UNFUNDED，不排期",
    )

    t = by_id["T12"]
    t.update(
        dependencies="T10;T11（若适用）",
        execution_standard="主竞争者仅 learned generic bilinear；fixed random 是弱对照。逐 GO mask 按 T10 的双口径三态量词裁决",
        acceptance_standard="每个适用 mask 有 GO/NO-GO/INCONCLUSIVE；无 A1 GO 时有 signed NOT_APPLICABLE verdict",
    )

    t = by_id["T13"]
    t.update(
        execution_standard="A2 存活算子逐一展开为独立行；L=32768/131072/262144/524288；524288 为 held-out，不参与算子或 slope 模型选择",
    )

    t = by_id["T14"]
    t.update(
        execution_standard="四条同时成立才放行：斜率交互项 CI 有利且排除 0；优势在 held-out 最大 L 复现；预注册计算的 L_star 落在已签数据与算力内；双口径同向且当前精度在 delta_NI 内非劣",
        failure_action="任一不中：不放行更大 L；保存 L_star 诊断值但不得将其作为扩容理由或未来必然交叉的论据",
    )

    t = by_id["T15"]
    t.update(
        dependencies="T10;T12/T14 的适用 verdict 或 signed NOT_APPLICABLE verdict",
        acceptance_standard="A1 verdict 与所有适用的 A2/A3 verdict 齐备；跳过阶段必须有 signed NOT_APPLICABLE；claim 与数据分支一致",
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0].keys())
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
