"""Verify C additions without the producer's predicate, brackets or canonicalizer.

Recompute old relations with the fresh sphere auditor, derive C from those
triple signs, expand each new 6x6 determinant by 720 permutations, and compare
raw C relations by complete backtracking. This emits a receipt, not a ledger.
"""
import hashlib,json,math,time
from fractions import Fraction as Q
from functools import reduce
from itertools import combinations
from pathlib import Path
from septic_three_obligations_v1 import POINTS,Relation,check,dump,sphere_geometry,relation_from_model,from_word,isomorphic
from septic_explicit_nonsingularity_v1 import audit_height
from protected_conic_invariant import prepared_rows,leading
ROOT=Path(__file__).resolve().parent

def C_relation(model,base):
    labels=sorted(map(int,model['colors']));n=base.n;inner={i for i,c in enumerate(base.colors) if c.startswith('i')}
    fours={}
    for four in combinations(range(n),4):
        ins=set(four)&inner
        if len(ins)!=2:continue
        pos={four[j] for j in range(4) if ((-1)**j)*base.table[four[:j]+four[j+1:]]>0}
        if pos==ins or pos==set(four)-ins:fours[four]={'four':[labels[v] for v in four],'positive_cofactor_labels':[labels[v] for v in sorted(pos)]}
    eligible={}
    for six in combinations(range(n),6):
        for four in combinations(six,4):
            if four in fours:eligible[six]=fours[four];break
    rows,scales,ed=prepared_rows(model);signs={};records=[]
    for six,certificate in eligible.items():
        rec=leading(tuple(labels[v] for v in six),rows,scales,ed)
        signs[six]=rec['sign'];records.append({**rec,**certificate})
    return Relation(base.colors,base.triples,signs),records
