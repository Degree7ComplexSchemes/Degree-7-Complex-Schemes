"""Exact protected Veronese signs and an intrinsic colored canonical form.

Protection is the sufficient criterion in protected_sextuple_conic_lemma.txt.
This module alone does not assert any additional algebraic curve class.
"""
from collections import defaultdict
from fractions import Fraction as Q
from functools import cmp_to_key
from itertools import combinations, permutations
from math import lcm, prod
import json

PAIRS=((0,0),(0,1),(1,1),(0,2),(1,2),(2,2))
PERMS=[(p,(-1)**sum(p[i]>p[j] for i in range(6) for j in range(i+1,6)))
       for p in permutations(range(6))]

def parity(order):
    return (-1)**sum(order[i]>order[j] for i in range(len(order)) for j in range(i+1,len(order)))

def triple_data(model):
    colors={int(v):c for v,c in model['colors'].items()}
    inner=set(model['inner_regions'])
    signs={tuple(r['regions']):r['sign'] for r in model['triple_determinants']}
    def chi(a,b,c):
        t=(a,b,c)
        return signs[tuple(sorted(t))]*parity(t) if inner.intersection(t) else 0
    return colors,inner,chi

def qualifying(model):
    colors,inner,chi=triple_data(model);outer=set(colors)-inner;out={}
    for k in (4,5):
        for a in combinations(sorted(inner),k):
            for b in combinations(sorted(outer),6-k):
                six=tuple(sorted(a+b))
                for four in combinations(six,4):
                    pos=sum((-1)**j*chi(*(four[:j]+four[j+1:]))>0 for j in range(4))
                    if pos!=2:
                        out[six]={'nonconvex_four':list(four),'positive_cofactors':pos};break
    return out

def minimizing_orders(model):
    """Retain all orders minimizing the old full protected triple key."""
    colors,inner,chi=triple_data(model);vertices=sorted(colors);best=None;orders=set()
    for root in sorted(inner):
        others=[v for v in vertices if v!=root]
        if len(others)<2:circle=others
        else:
            anchor=others[0];half=[(v,chi(root,anchor,v)) for v in others[1:]]
            det=lambda a,b:a[1]*b[1]*chi(root,a[0],b[0])
            half.sort(key=cmp_to_key(lambda a,b:-det(a,b)))
            assert all(det(a,b)==1 for a,b in combinations(half,2))
            rays=[(anchor,1)]+half+[(anchor,-1)]+[(v,-s) for v,s in half]
            circle=[v for v,s in rays if s==1]
        for direction in (circle,list(reversed(circle))):
            for cut in range(max(1,len(direction))):
                order=(root,*direction[cut:],*direction[:cut])
                values=[chi(*t) for t in combinations(order,3)]
                first=next((v for v in values if v),1)
                word=''.join('+' if first*v>0 else '-' if first*v<0 else '.' for v in values)
                key=json.dumps([colors[v] for v in order],separators=(',',':'))+'|'+word
                if best is None or key<best:best=key;orders={order}
                elif key==best:orders.add(order)
    return best,sorted(orders)

def prepared_rows(model):
    points={int(v):p['vector'] for v,p in model['centers'].items() if int(v) in map(int,model['colors'])}
    ed=lcm(*(Q(e).denominator for row in points.values() for c,e in row))
    rows={};scales={}
    for v,row in points.items():
        cs=[Q(c) for c,e in row];es=[int(Q(e)*ed) for c,e in row]
        scale=lcm(*(c.denominator for c in cs));cs=[int(c*scale) for c in cs]
        rows[v]=[(cs[i]*cs[j],es[i]+es[j]) for i,j in PAIRS];scales[v]=scale
    return rows,scales,ed

def leading(six,rows,scales,ed):
    data=[rows[v] for v in six];terms=defaultdict(int)
    for p,sg in PERMS:
        exponent=0;coefficient=sg
        for i,j in enumerate(p):
            c,e=data[i][j];exponent+=e;coefficient*=c
        terms[exponent]+=coefficient
    nonzero=[e for e,c in terms.items() if c]
    assert nonzero,('Protected conic determinant vanished identically',six)
    exponent=min(nonzero);coefficient=Q(terms[exponent],prod(scales[v]**2 for v in six))
    return {'regions':list(six),'sign':1 if coefficient>0 else -1,
            'leading_exponent':str(Q(exponent,ed)),'leading_coefficient':str(coefficient)}

def canonical(model,records):
    base,orders=minimizing_orders(model)
    signs={tuple(r['regions']):r['sign'] for r in records};best=None
    for order in orders:
        word=''.join('.' if tuple(sorted(t)) not in signs else
            ('+' if parity(t)*signs[tuple(sorted(t))]>0 else '-') for t in combinations(order,6))
        best=word if best is None else min(best,word)
    return base+'|'+best

def evaluate(model):
    eligible=qualifying(model);rows,scales,ed=prepared_rows(model)
    records=[{**leading(six,rows,scales,ed),**criterion} for six,criterion in sorted(eligible.items())]
    return {'protected_conic_determinants':records,'conic_augmented_signature':canonical(model,records)}
