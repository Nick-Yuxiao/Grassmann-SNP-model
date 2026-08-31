from __future__ import annotations
import argparse,hashlib,json,os,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def canonical_hash(p): return hashlib.sha256(("\n".join(p.read_text().splitlines())+"\n").encode()).hexdigest()
def scale_family(f,scale,SharedFamily):
 extra=np.zeros_like(f.y)
 for x,b,gamma in zip(f.regions,f.B,f.Gamma): extra+=(scale-1)*(x@b+(f.g[:,None]*x)@gamma)
 return SharedFamily(f.subject_ids,f.y+extra,f.g,f.covariates,f.regions,tuple(scale*b for b in f.B),tuple(scale*g for g in f.Gamma),f.seed)
def atomic(p,v):
 t=p.with_suffix(p.suffix+".tmp"); t.write_text(json.dumps(v,indent=2,sort_keys=True)+"\n"); t.replace(p)
def main():
 ap=argparse.ArgumentParser(); ap.add_argument("--output-dir",required=True); ap.add_argument("--parent-package"); ap.add_argument("--r1-run"); a=ap.parse_args()
 if os.environ.get("CUDA_VISIBLE_DEVICES") not in (None,"","-1"): raise RuntimeError("CPU only")
 cfg=json.loads((ROOT/"config"/"RANK_GATE_CONFIG.json").read_text()); meta=json.loads((ROOT/"PARENT_EVIDENCE.json").read_text())
 parent=Path(a.parent_package).resolve() if a.parent_package else (ROOT.parent/meta["code_parent"]).resolve()
 if canonical_hash(parent/"MANIFEST.sha256")!=meta["code_parent_manifest_canonical_sha256"]: raise RuntimeError("R1 parent mismatch")
 if a.r1_run and hashlib.sha256((Path(a.r1_run)/"RESULT_MANIFEST.sha256").read_bytes()).hexdigest()!=meta["r1_result_manifest_sha256"]: raise RuntimeError("R1 evidence mismatch")
 sys.path[:0]=[str(parent/"src"),str(parent.parent/json.loads((parent/"PARENT_EVIDENCE.json").read_text())["code_parent"]/"src"),str(ROOT/"src")]
 from shared_family_r1 import SharedFamily,generate_shared_family,reorder_candidates,reorder_subjects,run_synchronized_maxT
 from rank_gate_r1_5 import run_rank_gated_maxT
 kw=dict(resamples=cfg["resamples"],seed=cfg["bootstrap_seed"],rank=cfg["rank"],ridge_lambda=cfg["ridge_lambda"])
 strong=scale_family(generate_shared_family(seed=cfg["strong_seed"],n=cfg["n"],family_size=4,conditional_ld=True),2,SharedFamily)
 weak=scale_family(generate_shared_family(seed=cfg["weak_seed"],n=cfg["n"],family_size=4,conditional_ld=True),2,SharedFamily)
 ungated=run_synchronized_maxT(strong,**kw); zero=run_rank_gated_maxT(strong,minimum_gap=0.0,**kw); gated=run_rank_gated_maxT(weak,minimum_gap=cfg["minimum_relative_rank_gap"],**kw)
 rev=run_rank_gated_maxT(reorder_candidates(weak,np.arange(3,-1,-1)),minimum_gap=cfg["minimum_relative_rank_gap"],**kw)
 order=np.random.default_rng(109000).permutation(cfg["n"]); perm=run_rank_gated_maxT(reorder_subjects(weak,order),minimum_gap=cfg["minimum_relative_rank_gap"],**kw)
 ineligible=~gated.observed_eligible
 checks={
  "threshold_zero_observed_equals_r1":bool(np.allclose(zero.observed,ungated.observed)),
  "threshold_zero_resamples_equal_r1":bool(np.allclose(zero.resampled,ungated.resampled)),
  "threshold_zero_family_p_equal_r1":zero.family_p_value==ungated.family_p_value,
  "weak_rank_control_contains_ineligible":bool(ineligible.any()),
  "ineligible_observed_statistics_zero":bool(np.all(gated.observed[ineligible]==0)),
  "ineligible_candidate_p_values_one":bool(np.all(gated.candidate_p_values[ineligible]==1)),
  "bootstrap_gate_applied":bool(np.all(gated.resampled[~gated.resampled_eligible]==0)),
  "family_dimension_retained":len(gated.observed)==cfg["family_size"],
  "candidate_order_invariant":bool(np.allclose(gated.resampled_max,rev.resampled_max) and gated.family_p_value==rev.family_p_value),
  "subject_order_invariant":bool(np.allclose(gated.resampled_max,perm.resampled_max,rtol=1e-10,atol=1e-10) and gated.family_p_value==perm.family_p_value),
  "multiplier_fingerprints_invariant":gated.multiplier_fingerprints==rev.multiplier_fingerprints==perm.multiplier_fingerprints}
 status="R1_5_RANK_GATE_PASS" if all(checks.values()) else "R1_5_RANK_GATE_FAIL"
 gate={"gate":"R1.5-rank-gated-statistic","status":status,"minimum_relative_rank_gap":cfg["minimum_relative_rank_gap"],"checks":checks,
 "details":{"observed_raw":gated.observed_raw.tolist(),"observed_gaps":gated.observed_gaps.tolist(),"observed_eligible":gated.observed_eligible.tolist(),"candidate_p_values":gated.candidate_p_values.tolist(),"family_p":gated.family_p_value},
 "authorizes":["bounded_smoke_redesign_new_seed_block"] if status.endswith("PASS") else [],"does_not_authorize":cfg["does_not_authorize"]}
 out=Path(a.output_dir).resolve()
 if out.exists(): raise FileExistsError("duplicate output refused")
 out.mkdir(parents=True); atomic(out/"GATE_R1_5_RANK_GATE.json",gate); atomic(out/"ENVIRONMENT.json",{"python":sys.version,"numpy":np.__version__,"accelerator_used":False}); print(json.dumps(gate,indent=2,sort_keys=True)); return 0 if status.endswith("PASS") else 1
if __name__=="__main__": raise SystemExit(main())

