"""Independent geometry/arithmetic and complete joint A/B/C separation audit.

This verifier does not call a production eligibility predicate, bracket
shortcut, or canonicalizer. All parent representatives and all previously
accepted C additions in a candidate's OLD key bucket are compared.
"""
import argparse,hashlib,json,math,time
from collections import defaultdict
from fractions import Fraction as Q
from functools import reduce
from pathlib import Path
from septic_three_obligations_v1 import POINTS,Relation,check,dump,sphere_geometry,relation_from_model,from_word,isomorphic
from septic_explicit_nonsingularity_v1 import audit_height
from septic_verify_nest_chord_additions_v1 import C_relation
ROOT=Path(__file__).resolve().parent
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def joint(m):
    b,a=relation_from_model(m);c,recs=C_relation(m,b)
    check(all(k not in b.conics or b.conics[k]==v for k,v in c.conics.items()),'overlap sign mismatch')
    return b,Relation(b.colors,b.triples,{**b.conics,**c.conics}),a,recs
