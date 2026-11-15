"""Regional crossing blocks of Cimasoni--Ferretti--Liu, Figure 2."""
from fractions import Fraction as F
import itertools as it
import math

def tau(k, l):
    if k < l:
        return list(range(-(k + 1), -l - 1, -1)) + list(range(k, l))
    if k > l:
        return list(range(-(k - 1), -l + 1)) + list(range(k, l, -1))
    return []

def pencil_braid(bands, surgery=False):
    result = [-1] if surgery else []
    path = (2,) + tuple(bands) + (5,)
    for k, l in zip(path, path[1:]):
        result.extend([-k] + tau(k, l))
    for last in range(6, 0, -1):
        result.extend(range(1, last + 1))
    return result

def integer_inertia(symmetric):
    a = [list(row) for row in symmetric]
    denominator = 1
    positive = negative = 0
    pivots = []
    while a:
        k = next((j for j in range(len(a)) if a[j][j]), None)
        shear = None
        if k is None:
            edge = next(((i, j) for i in range(len(a)) for j in range(i + 1, len(a)) if a[i][j]), None)
            if edge is None:
                return (positive, negative, len(a)), pivots
            i, j = edge
            a[i] = [x + y for x, y in zip(a[i], a[j])]
            for row in a:
                row[i] += row[j]
            k, shear = i, edge
        pivot = a[k][k]
        positive += pivot * denominator > 0
        negative += pivot * denominator < 0
        pivots.append({"index": k, "shear": shear, "numerator": str(pivot), "denominator": str(denominator)})
        indices = [j for j in range(len(a)) if j != k]
        result = [[0] * len(indices) for _ in indices]
        for p, i in enumerate(indices):
            for q in range(p, len(indices)):
                j = indices[q]
                value, remainder = divmod(pivot * a[i][j] - a[i][k] * a[k][j], denominator)
                if remainder:
                    raise ArithmeticError("Nonintegral fraction-free elimination")
                result[p][q] = result[q][p] = value
        a, denominator = result, pivot
    return (positive, negative, 0), pivots

def model(word,palette,reverse=False):
    n=len(palette); col=list(palette); strands=list(range(n)); gaps=list(range(n+1))
    parent=list(range(n+1+len(word)))
    def find(i):
        while parent[i]!=i: i=parent[i]
        return i
    events=[]; wm=0; lk=[[0]*n for _ in range(n)]
    for t,g in enumerate(word):
        i=abs(g); sign=1 if g>0 else -1
        j,k=col[i-1:i+1]; a,d=strands[i-1:i+1]
        lk[a][d]+=sign;lk[d][a]+=sign
        fresh=n+1+t
        # Incoming gap, right side, outgoing gap, left side.
        corners=(gaps[i],gaps[i+1],fresh,gaps[i-1])
        if reverse:corners=(gaps[i],gaps[i-1],fresh,gaps[i+1])
        events.append((corners,sign,j,k)); gaps[i]=fresh
        col[i-1:i+1]=[k,j];strands[i-1:i+1]=[d,a]
        wm+=sign*(j==k)
    assert col==list(palette),(word,palette,col)
    for i,end in enumerate(gaps):parent[find(end)]=find(i)
    region_ids=sorted({find(i) for i in range(len(parent))})
    assert len(region_ids)==len(word)+2
    mark=palette.index(0)
    region_ids=[i for i in region_ids if i not in (find(mark),find(mark+1))]
    ix={j:i for i,j in enumerate(region_ids)}
    events=[(tuple(ix.get(find(i),-1) for i in r),s,j,k) for r,s,j,k in events]
    cycles=[]; todo=set(range(n))
    while todo:
        a=min(todo); c=[]
        while a in todo:
            c.append(a);todo.remove(a);a=strands[a]
        cycles.append(c)
    lk=[[sum(lk[i][j] for i in c for j in d)//2 if c!=d else 0 for d in cycles] for c in cycles]
    return len(region_ids),wm,events,[[i+1 for i in c] for c in cycles],lk,[palette[c[0]] for c in cycles]

def blocks(half):
    blocks={}
    for j,k in it.product(range(len(half)),repeat=2):
        x,y=half[j]; z,v=half[k]
        c=x*z-y*v; d=x*z+y*v
        blocks[j,k]=[[a/(y*v) for a in row] for row in
          [[c,x,1,z],[x,d,z,1],[1,z,c,x],[z,1,x,d]]]
    return blocks

def matrix(m,half):
    n,_,events,*_=m
    H=[[F(0) for j in range(n)] for i in range(n)]
    block=blocks(half)
    for r,s,j,k in events:
        A=block[j,k]
        for a,ra in enumerate(r):
            if ra<0:continue
            for b,rb in enumerate(r):
                if rb>=0:H[ra][rb]+=s*A[a][b]
    return H

def half(q):
    return ((1-q*q)/(1+q*q),2*q/(1+q*q))

def exact(m,q,r):
    A=matrix(m,[half(q),half(r)])
    d=math.lcm(*(x.denominator for row in A for x in row))
    integer=[[int(x*d) for x in row] for row in A]
    inertia,pivots=integer_inertia(integer)
    plus,minus,zero=inertia
    assert (plus-minus-m[1])%2==0 and zero%2==0
    return dict(q=str(q),r=str(r),inertia=inertia,writhe=m[1],sigma=(plus-minus-m[1])//2,eta=zero//2,pivots=pivots)

