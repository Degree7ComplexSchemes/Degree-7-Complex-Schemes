"""Recompute all fourteen entry-17 Seifert matrices and exact inertia counts.

This verifies the table algebra for the supplied pencil words. Exhaustion
of geometrically possible pencils remains a separate proof obligation.
"""
from pathlib import Path
import argparse,importlib.util,json,sys
from fractions import Fraction

ROOT=Path(__file__).resolve().parents[1]
path=ROOT/'obstructions/entry2/scripts/audit_table19_entry2_cover.py'
spec=importlib.util.spec_from_file_location('entry17_braid_algebra',path)
algebra=importlib.util.module_from_spec(spec);spec.loader.exec_module(algebra)
CASES=[(3,'334343'),(3,'343434'),(3,'433343'),(3,'434333'),(3,'434344'),(3,'434443'),(3,'444343'),
       (4,'333434'),(4,'343334'),(4,'343433'),(4,'343444'),(4,'344434'),(4,'434343'),(4,'443434')]

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    if not __debug__:raise SystemExit('Do not use Python -O')
    if args.output.exists():raise SystemExit('Choose a new output path')
    rows=[]
    for index,(initial,bands) in enumerate(CASES,1):
        word=[];steps=[initial,*map(int,bands),7-initial]
        for k,l in zip(steps,steps[1:]):word.extend([-k]+algebra.tau(k,l))
        word+= [j for end in range(6,0,-1) for j in range(1,end+1)]
        v=algebra.seifert_matrix(word,7)
        symmetric=[[v[i][j]+v[j][i] for j in range(len(v))] for i in range(len(v))]
        inertia,pivots=algebra.integer_inertia(symmetric)
        assert inertia==algebra.rational_inertia(symmetric)==(11,21,0),(index,inertia)
        determinant=Fraction(1)
        for pivot in pivots:determinant*=Fraction(int(pivot['numerator']),int(pivot['denominator']))
        assert determinant==-80 and len(v)==32
        assert sum(1 if g>0 else -1 for g in word)==14
        rows.append({'case':index,'initial':initial,'bands':bands,'terminal':7-initial,
            'braid_word':word,'strands':7,'seifert_matrix':v,'symmetric_matrix':symmetric,
            'determinant':int(determinant),'inertia':list(inertia),'signature':-10,'nullity_plus_one':1,
            'elimination_pivots':pivots})
    result={'all_pass':True,'cases':14,'method':'Integer fraction-free congruence and independent rational congruence; full matrices and pivots retained.',
        'scope':'Algebra for the fourteen supplied words; not a proof of geometric exhaustion or an independent derivation of the Seifert-matrix formula.',
        'records':rows}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('PASS: all 14 matrices have determinant -80 and inertia (11,21,0).',flush=True)

if __name__=='__main__':main()
