"""Fresh sphere-cover, bracket-identity, and exhaustive-isomorphism audit.

No old geometry, determinant, cyclic-order, or canonicalization module is
imported. Published geometric theorems are supplied in the companion proof.
"""
import argparse
import gzip
import hashlib
import json
import math
import time
from collections import Counter, defaultdict, deque
from fractions import Fraction as Q
from itertools import combinations, permutations, product
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parent
POINTS = [(i, j) for i in range(8) for j in range(8-i)]
PERMS3 = [(p, (-1)**sum(p[i] > p[j] for i in range(3) for j in range(i+1, 3)))
          for p in permutations(range(3))]


def check(b, message):
    if not b:
        raise ValueError(message)


def dump(x):
    return json.dumps(x, separators=(",", ":"), sort_keys=True)


def neg(v):
    return tuple(-x for x in v)


def add(a, b):
    return tuple(x+y for x, y in zip(a, b))


def det3(a, b, c):
    return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])


class UF:
    def __init__(self):
        self.p = []
    def new(self):
        n = len(self.p); self.p.append(n); return n
    def find(self, a):
        while self.p[a] != a:
            self.p[a] = self.p[self.p[a]]; a = self.p[a]
        return a
    def join(self, a, b):
        self.p[self.find(a)] = self.find(b)


def sphere_geometry(w, model, typ, scheme):
    """Eight octants on S^2; no orientation-transition bits are propagated."""
    tris = [tuple(map(tuple, t)) for t in w["triangles"]]
    phase = dict(zip(POINTS, map(int, w["phase_mask"])))
    owners = defaultdict(list)
    for k, t in enumerate(tris):
        for e in combinations(t, 2): owners[tuple(sorted(e))].append(k)
    # Compute Harnack compatibility by enumerating its defining distributions.
    harnacks = [{p: mu ^ int(p != special) for p in product((0, 1), repeat=2)}
                for special in product((0, 1), repeat=2) for mu in (0, 1)]
    adj = defaultdict(list)
    for e, ks in owners.items():
        if len(ks) != 2: continue
        a, b = ks; union = set(tris[a]) | set(tris[b])
        compatible = any(all(h[x % 2, y % 2] == phase[x, y] for x, y in union) for h in harnacks)
        adj[a].append((b, int(not compatible))); adj[b].append((a, int(not compatible)))
    theta = {0: 0}; queue = deque([0]); consistent = True
    while queue:
        a = queue.popleft()
        for b, bit in adj[a]:
            if b not in theta: theta[b] = theta[a] ^ bit; queue.append(b)
            elif theta[b] != theta[a] ^ bit: consistent = False
    check(len(theta) == 49 and (typ == "I") == consistent, "Klein type mismatch")
    faces = UF(); edge_faces = defaultdict(list); vertex_faces = defaultdict(list)
    arcs = []; endpoint_arcs = defaultdict(list); face_vertices = []
    for ti, original in enumerate(tris):
        t = list(original)
        if (t[1][0]-t[0][0])*(t[2][1]-t[0][1])-(t[1][1]-t[0][1])*(t[2][0]-t[0][0]) < 0:
            t[1], t[2] = t[2], t[1]
        if consistent and theta[ti]: t[1], t[2] = t[2], t[1]
        for sx, sy, sz in product((1, -1), repeat=3):
            vs = [(sx*x, sy*y, sz*(7-x-y)) for x, y in t]
            bits = [phase[x, y] ^ ((sx < 0)*(x % 2)) ^ ((sy < 0)*(y % 2))
                    ^ ((sz < 0)*((7-x-y) % 2)) for x, y in t]
            part = {}
            for bit in set(bits):
                f = faces.new(); part[bit] = f
                fv = [v for v, b in zip(vs, bits) if b == bit]
                face_vertices.append(fv)
                for v in fv: vertex_faces[v].append(f)
            for i, j in combinations(range(3), 2):
                edge_faces[tuple(sorted((vs[i], vs[j])))].append({vs[i]: part[bits[i]], vs[j]: part[bits[j]]})
            if len(part) == 2:
                i = next(i for i in range(3) if bits[i] == bits[(i+1) % 3])
                x, y, z = vs[i], vs[(i+1) % 3], vs[(i+2) % 3]
                u, v = add(x, z), add(y, z)
                k = len(arcs); arcs.append((u, v, part[bits[(i+2) % 3]], part[bits[i]]))
                endpoint_arcs[u].append(k); endpoint_arcs[v].append(k)
    for edge, entries in edge_faces.items():
        check(len(entries) == 2, "sphere edge incidence")
        for v in edge: faces.join(entries[0][v], entries[1][v])
    check(all(len(ks) == 2 for ks in endpoint_arcs.values()), "curve is not a closed one-manifold")
    curve_uf = UF()
    for _ in arcs: curve_uf.new()
    for ks in endpoint_arcs.values(): curve_uf.join(*ks)
    components = defaultdict(list)
    for k in range(len(arcs)): components[curve_uf.find(k)].append(k)
    endpoint_component = {p: curve_uf.find(ks[0]) for p, ks in endpoint_arcs.items()}
    one_sided = [c for c, ks in components.items() if endpoint_component[neg(arcs[ks[0]][0])] == c]
    check(len(one_sided) == 1, "unique lifted pseudoline")
    J = one_sided[0]
    # Quotient complementary regions by antipodes to obtain RP^2 regions.
    region_anti = {}
    for v, fs in vertex_faces.items():
        r = faces.find(fs[0]); ar = faces.find(vertex_faces[neg(v)][0])
        check(all(faces.find(f) == r for f in fs), "vertex region ambiguity")
        check(r != ar and (r not in region_anti or region_anti[r] == ar), "region antipode")
        region_anti[r] = ar
    def rp(r): return min(r, region_anti[r])
    tree = defaultdict(set); root_candidates = set()
    for u, v, a, b in arcs:
        a, b = rp(faces.find(a)), rp(faces.find(b))
        if a == b: root_candidates.add(a)
        else: tree[a].add(b); tree[b].add(a)
    check(len(root_candidates) == 1, "root complementary region")
    root = next(iter(root_candidates)); parent = {root: None}; queue = deque([root])
    while queue:
        a = queue.popleft()
        for b in tree[a]:
            if b not in parent: parent[b] = a; queue.append(b)
            else: check(parent[a] == b or parent[b] == a, "region adjacency not a tree")
    check(set(parent) == {rp(r) for r in region_anti}, "missing complementary region")
    # Join across all oval arcs, leaving ONLY the lifted pseudoline cut.
    sheet = UF()
    sheet.p = faces.p.copy()
    for k, (_, _, a, b) in enumerate(arcs):
        if curve_uf.find(k) != J: sheet.join(a, b)
    check(len({sheet.find(f) for f in range(len(sheet.p))}) == 2, "pseudoline complement on sphere")
    # Select the sheet in which oriented J has its interior on the left.
    k = components[J][0]; u, v, a, b = arcs[k]
    left = a if det3(u, v, face_vertices[a][0]) > 0 else b
    preferred = sheet.find(left)
    signs = {}
    if consistent:
        outgoing = Counter(u for u, v, a, b in arcs); incoming = Counter(v for u, v, a, b in arcs)
        check(outgoing == incoming and all(n == 1 for n in outgoing.values()), "symmetric arc orientation")
        for u, v, a, b in arcs:
            if sheet.find(a) != preferred: continue
            ra, rb = rp(faces.find(a)), rp(faces.find(b))
            if ra == rb: continue
            f, r = (a, ra) if parent[ra] == rb else (b, rb)
            value = "-" if det3(u, v, face_vertices[f][0]) > 0 else "+"
            check(r not in signs or signs[r] == value, "complex orientation inconsistency")
            signs[r] = value
    children = {r: [s for s in parent if parent[s] == r] for r in parent}
    containers = [r for r in children[root] if children[r]]
    signed_counts = lambda rs: [sum(signs[r] == s for r in rs) for s in ("+", "-")]
    if not containers:
        actual_scheme = ["unnested", signed_counts(children[root]) if consistent else len(children[root])]
    elif len(containers) == 1 and not any(children[r] for r in children[containers[0]]):
        container = containers[0]; exterior = [r for r in children[root] if r != container]
        actual_scheme = (["nest", signed_counts(exterior), signs[container], signed_counts(children[container])]
                         if consistent else ["nest", len(exterior), len(children[container])])
    else:
        check(len(children[root]) == 1, "unsupported multiple-nest scheme")
        chain = []; r = children[root][0]
        while True:
            chain.append(r)
            if not children[r]: break
            check(len(children[r]) == 1, "unsupported branching nest")
            r = children[r][0]
        check(len(chain) == 3, "unsupported nest depth")
        actual_scheme = ["deep", [1 if signs[r] == "+" else -1 for r in chain]] if consistent else ["deep"]
    check(actual_scheme == scheme, "full intrinsic scheme mismatch")
    if model is None:
        return {"sphere_curve_components": len(components), "regions": len(parent), "marked_points": 0}
    regions = {}; chosen_sheets = set(); height = {tuple(p): Q(h) for p, h in w["heights"]}
    for label, center in model["centers"].items():
        i, j = center["vertex"]; sx, sy = center["quadrant"]; ls = center["lift_sign"]
        vertex = (sx*i, sy*j, 7-i-j)
        f = vertex_faces[vertex][0]; regions[int(label)] = rp(faces.find(f))
        lifted = vertex if ls == 1 else neg(vertex)
        chosen_sheets.add(sheet.find(vertex_faces[lifted][0]))
        ax, ay = map(Q, center["amplitudes"]); vector = center["vector"]
        check(ax > 0 and ay > 0 and [Q(q[0]) for q in vector] == [ls*sx*ax, ls*sy*ay, ls], "point coefficients")
        alpha, beta = Q(vector[0][1]), Q(vector[1][1])
        check(Q(vector[2][1]) == 0, "point Z exponent")
        gap = min(height[p]-height[i, j]+alpha*(p[0]-i)+beta*(p[1]-j) for p in POINTS if p != (i, j))
        check(gap > 0 and gap == Q(center["dominance_gap"]), "dominant point gap")
    check(len(chosen_sheets) == 1, "incoherent marked-point lifts")
    check(len(regions) == len(set(regions.values())) == len(parent)-1 and root not in regions.values(), "marked region bijection")
    inner = set(model["inner_regions"]); outer = set(model["outer_regions"])
    check(len(containers) == 1, "single-nest scope")
    container = containers[0]
    check({regions[v] for v in inner} == {r for r in parent if parent[r] == container}, "inner point regions")
    check({regions[v] for v in outer} == {r for r in parent if parent[r] == root and r != container}, "exterior point regions")
    check(set(map(int, model["colors"])) == inner | outer and not inner & outer, "colored oval set")
    for label in inner | outer:
        color = ("i" if label in inner else "o") + (signs[regions[label]] if consistent else "")
        check(model["colors"][str(label)] == color, "intrinsic complex color")
    actual = ["nest", len(outer), len(inner)] if not consistent else ["nest", [sum(signs[regions[v]] == s for v in outer) for s in ("+", "-")], signs[container], [sum(signs[regions[v]] == s for v in inner) for s in ("+", "-")]]
    check(actual == scheme, "scheme mismatch")
    return {"sphere_curve_components": len(components), "regions": len(parent), "marked_points": len(regions), "coherent_sphere_sheet": True, "dominance_and_colors": True}


def parity(t):
    return (-1)**sum(a > b for i, a in enumerate(t) for b in t[i+1:])


def poly_mul(a, b):
    out = defaultdict(int)
    for e, c in a.items():
        for f, d in b.items(): out[e+f] += c*d
    return {e: c for e, c in out.items() if c}


def relation_from_model(model):
    labels = sorted(map(int, model["colors"])); colors = [model["colors"][str(v)] for v in labels]
    inner = {i for i, c in enumerate(colors) if c.startswith("i")}; n = len(labels)
    vectors = [model["centers"][str(v)]["vector"] for v in labels]
    D = math.lcm(*(Q(e).denominator for v in vectors for c, e in v))
    rows = []
    for v in vectors:
        scale = math.lcm(*(Q(c).denominator for c, e in v))
        rows.append([(int(Q(c)*scale), int(Q(e)*D)) for c, e in v])
    brackets = {}
    chi = {}
    for abc in combinations(range(n), 3):
        pol = defaultdict(int)
        for perm, s in PERMS3:
            c = s; e = 0
            for i, j in enumerate(perm):
                a, b = rows[abc[i]][j]; c *= a; e += b
            pol[e] += c
        pol = {e: c for e, c in pol.items() if c}; brackets[abc] = pol
        if inner.intersection(abc):
            check(bool(pol), "zero protected triple")
            chi[abc] = 1 if pol[min(pol)] > 0 else -1
    def bracket(t):
        p = brackets[tuple(sorted(t))]
        return p if parity(t) == 1 else {e: -c for e, c in p.items()}
    def oriented(t): return chi.get(tuple(sorted(t)), 0)*parity(t)
    conics = {}; slow = 0
    for six in combinations(range(n), 6):
        if not 0 < len(inner.intersection(six)) < 6: continue
        eligible = False
        for four in combinations(six, 4):
            if len(inner.intersection(four)) >= 2 and sum(((-1)**i)*oriented(four[:i]+four[i+1:]) > 0 for i in range(4)) != 2:
                eligible = True; break
            for r in inner.intersection(four):
                a, b, c = [v for v in four if v != r]
                values = [oriented((r, a, b)), oriented((r, b, c)), oriented((r, c, a))]
                if values[0] == values[1] == values[2]: eligible = True; break
            if eligible: break
        if not eligible: continue
        a, b, c, d, e, f = six
        # Exact Veronese determinant identity in columns x2,xy,y2,xz,yz,z2.
        plus = [bracket(t) for t in ((a,b,d),(a,c,e),(b,c,f),(d,e,f))]
        minus = [bracket(t) for t in ((a,b,c),(a,d,e),(b,d,f),(c,e,f))]
        def leading(ps):
            if any(not p for p in ps): return None
            degrees = [min(p) for p in ps]
            return sum(degrees), math.prod(p[k] for p, k in zip(ps, degrees))
        lp, lm = leading(plus), leading(minus)
        if lp is not None and (lm is None or lp[0] < lm[0]): coeff = lp[1]
        elif lm is not None and (lp is None or lm[0] < lp[0]): coeff = -lm[1]
        elif lp is not None and lp[1] != lm[1]: coeff = lp[1]-lm[1]
        else:
            slow += 1; p = {0: 1}; m = {0: 1}
            for q in plus: p = poly_mul(p, q)
            for q in minus: m = poly_mul(m, q)
            for exponent, coef in m.items(): p[exponent] = p.get(exponent, 0)-coef
            p = {ex: coef for ex, coef in p.items() if coef}
            check(bool(p), "zero protected conic"); coeff = p[min(p)]
        conics[six] = 1 if coeff > 0 else -1
    return Relation(colors, chi, conics), {"triples": len(chi), "conics": len(conics), "full_bracket_expansions": slow}


class Relation:
    def __init__(self, colors, triples, conics):
        self.colors = colors; self.n = len(colors); self.triples = triples; self.conics = conics
        self.table = {}
        for t, value in triples.items():
            for p in permutations(t): self.table[p] = value*parity(p)
        self.palette = sorted(set(colors)); self.pair = {}; self.vertex = {}
        for g in (1, -1):
            for a in range(self.n):
                observations = []
                for b in range(self.n):
                    if a == b: continue
                    count = Counter((colors[c], g*self.table.get((a,b,c), 0)) for c in range(self.n) if c not in (a,b))
                    fp = tuple(count[color, s] for color in self.palette for s in (-1,0,1))
                    self.pair[g,a,b] = fp; observations.append((colors[b], fp))
                self.vertex[g,a] = tuple(sorted(observations))


def from_word(word):
    a, b, c = word.split("|"); colors = json.loads(a); n = len(colors)
    sign = {"+": 1, "-": -1}
    ts = list(combinations(range(n), 3)); ss = list(combinations(range(n), 6))
    check(len(b) == len(ts) and len(c) == len(ss), "word dimensions")
    return Relation(colors, {t: sign[v] for t,v in zip(ts,b) if v != "."}, {t: sign[v] for t,v in zip(ss,c) if v != "."})


def isomorphic(A, B, stats):
    """Complete finite backtracking, not a cyclic order or canonical minimum."""
    if A.n != B.n or Counter(A.colors) != Counter(B.colors): return None
    n = A.n
    for g in (1, -1):
        base = {a: [b for b in range(n) if A.colors[a] == B.colors[b] and A.vertex[g,a] == B.vertex[1,b]] for a in range(n)}
        if any(not bs for bs in base.values()): continue
        mapping = {}; used = set()
        def search():
            stats["nodes"] += 1
            if len(mapping) == n:
                if len(A.conics) != len(B.conics): return None
                for six, value in A.conics.items():
                    target = tuple(mapping[v] for v in six)
                    if value != B.conics.get(tuple(sorted(target)), 0)*parity(target): return None
                return dict(mapping)
            choices = []
            for a in range(n):
                if a in mapping: continue
                domain = []
                for b in base[a]:
                    if b in used: continue
                    if any(A.pair[g,a,x] != B.pair[1,b,y] for x,y in mapping.items()): continue
                    if any(g*A.table.get((a,x,z),0) != B.table.get((b,mapping[x],mapping[z]),0) for x,z in combinations(mapping,2)): continue
                    domain.append(b)
                if not domain: return None
                choices.append((len(domain), a, domain))
            _, a, domain = min(choices)
            for b in domain:
                mapping[a] = b; used.add(b)
                result = search()
                if result is not None: return result
                used.remove(b); del mapping[a]
            return None
        result = search()
        if result is not None: return {"triple_multiplier": g, "mapping": result}
    return None
