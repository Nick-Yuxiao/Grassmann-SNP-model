from __future__ import annotations
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def main():
 cfg=json.loads((ROOT/"config"/"RANK_GATE_CONFIG.json").read_text()); lines=(ROOT/"MANIFEST.sha256").read_text().splitlines(); dec=(ROOT/"DECISIONS_R1_5.tsv").read_text().strip().splitlines()
 checks={"gap_0_10":cfg["minimum_relative_rank_gap"]==.10,"resamples_39":cfg["resamples"]==39,"five_decisions":len(dec)==6,"manifest_nonempty":len(lines)>=8}
 checks["manifest_files_match"]=all((ROOT/r).is_file() and hashlib.sha256((ROOT/r).read_bytes()).hexdigest()==h for h,r in (x.split("  ",1) for x in lines))
 status="PASS" if all(checks.values()) else "FAIL"; print(json.dumps({"package":ROOT.name,"status":status,"checks":checks,"manifest_entries":len(lines)},indent=2,sort_keys=True)); return 0 if status=="PASS" else 1
if __name__=="__main__": raise SystemExit(main())

