from __future__ import annotations
import json,sys,unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]; meta=json.loads((ROOT/"PARENT_EVIDENCE.json").read_text()); parent=ROOT.parent/meta["code_parent"]
sys.path[:0]=[str(parent/"src"),str(parent.parent/json.loads((parent/"PARENT_EVIDENCE.json").read_text())["code_parent"]/"src"),str(ROOT/"src")]
from shared_family_r1 import SharedFamily,generate_shared_family,run_synchronized_maxT
from rank_gate_r1_5 import run_rank_gated_maxT
def scale(f,s):
 e=np.zeros_like(f.y)
 for x,b,g in zip(f.regions,f.B,f.Gamma): e+=(s-1)*(x@b+(f.g[:,None]*x)@g)
 return SharedFamily(f.subject_ids,f.y+e,f.g,f.covariates,f.regions,tuple(s*b for b in f.B),tuple(s*g for g in f.Gamma),f.seed)
class RankGateTests(unittest.TestCase):
 @classmethod
 def setUpClass(c):
  c.kw=dict(resamples=19,seed=108000,rank=2,ridge_lambda=.01); c.f=scale(generate_shared_family(seed=101002,n=360,family_size=4,conditional_ld=True),2); c.g=run_rank_gated_maxT(c.f,minimum_gap=.10,**c.kw)
 def test_ineligible_retained_and_p_one(self):
  bad=~self.g.observed_eligible; self.assertTrue(bad.any()); self.assertEqual(len(self.g.observed),4); self.assertTrue(np.all(self.g.observed[bad]==0)); self.assertTrue(np.all(self.g.candidate_p_values[bad]==1))
 def test_bootstrap_gate_zeroes_ineligible(self): self.assertTrue(np.all(self.g.resampled[~self.g.resampled_eligible]==0))
 def test_threshold_zero_matches_r1(self):
  a=run_rank_gated_maxT(self.f,minimum_gap=0,**self.kw); b=run_synchronized_maxT(self.f,**self.kw); self.assertTrue(np.allclose(a.resampled,b.resampled)); self.assertEqual(a.family_p_value,b.family_p_value)
if __name__=="__main__": unittest.main()

