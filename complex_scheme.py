"""
complex_scheme.py  --  complex scheme detector for T-curves of odd degree.

Given a primitive triangulation Gamma of Delta_m (m odd) and a sign
distribution eps, decide the Klein type and, when type I, return the complex
scheme in Viro's recursive notation with every oval signed relative to J.

Handles arbitrary nesting depth.

---------------------------------------------------------------------------
CONVENTIONS, DERIVED RATHER THAN GUESSED
---------------------------------------------------------------------------

(1) Symmetric orientation.  Parenti Def. 2.7: adjacent triangles have
    compatible orientations (= they induce opposite orientations on the common
    edge) iff eps restricted to their union is a Harnack distribution.  With
    every reference orientation counterclockwise, two adjacent triangles in the
    plane ALWAYS induce opposite orientations on the shared edge, so r == 1 and
    the system is

        x_tau + x_tau' = h(tau,tau') + 1   (mod 2).

    h is computed by the closed form in `is_harnack_pair` below.

(2) Arc orientation inside a triangle.  Parenti Fig. 2.  Given the cyclic order
    of a nonempty triangle copy, let (X -> Y) be the empty edge directed
    consistently with that cyclic order and Z the remaining (isolated) vertex.
    The arc is directed  mid(XZ) -> mid(YZ).

    Well defined: cyclic rotation leaves the directed empty edge unchanged;
    reversing the cyclic order reverses it, hence reverses the arc.  Choosing
    the opposite global convention reverses EVERY arc, including J, which
    leaves all oval signs unchanged (see (4)).  So this choice is free.

(3) Orientation of RP^2 \\ J.  J is one-sided, so RP^2 \\ J is an open disc and
    carries exactly two orientations.  Model: the diamond |i|+|j| <= m with
    antipodal boundary identification.  The diamond is an oriented disc; the
    identification is orientation-REVERSING, which is what makes RP^2
    non-orientable.  So a local orientation is (plane orientation) x (-1)^(number
    of boundary crossings), and this is consistent precisely after cutting along
    J.  Implemented as a two-colouring o of the faces, propagated across face
    adjacencies with a flip at the boundary identification, merging across
    oval-arcs but NOT across J-arcs.

    INVARIANT A: propagation must be consistent (no odd cycle).
    INVARIANT B: o must REVERSE across every J-arc and agree across every
    oval-arc.  This is forced, and is a strong check: going from one side of a
    one-sided curve to the other reverses the ambient orientation.

    Global normalisation: fix the sign of o so that the induced boundary
    orientation of the disc projects onto the complex orientation of J.  This is
    the only place J's orientation enters, and it is what makes the answer
    canonical rather than defined up to a flip.

(4) Sign of an oval.  For odd degree, an oval O is NEGATIVE iff [O] = 2[J] in
    H_1 of the Moebius band RP^2 \\ Int(Delta_O).

    Equivalent computational form, used here: O is negative iff its complex
    orientation agrees with the boundary orientation of Delta_O induced by the
    canonical orientation of RP^2 \\ J.

    Derivation.  Let A = D_J \\ Delta_O, an annulus, oriented by D_J.  Then
    [O]_A + [J-hat] = 0 and [J-hat] = 2[J], so [O]_A = -2[J].  The A-induced
    orientation on O is the reverse of the Delta_O-induced one, hence
    "complex orientation = Delta_O-induced" gives [O] = +2[J], i.e. negative.

    Reversing the whole semiorientation reverses J too, and both sides of the
    comparison flip, so the sign is an invariant of the curve, not of the choice.

    INVARIANT C: every arc of a given oval must yield the same sign.

(5) Injective pairs.  For nested ovals o1 > o2, the annulus between them lies in
    D_J; orienting it by D_J gives o1 its Delta-induced orientation and o2 the
    reverse of its own.  Hence

        the pair is POSITIVE  <=>  o1 and o2 have OPPOSITE signs,

    for ANY nesting depth, not just immediate containment.  Used by the Rokhlin
    self-check.

(6) Rokhlin-Mishachev, odd degree m = 2k+1:

        2(Pi+ - Pi-) + Lambda+ - Lambda- = l - k(k+1)

    where l is the number of ovals, Lambda+/- count positive/negative ovals, and
    Pi+/- count injective pairs over ALL nested pairs.

    INVARIANT D: this must hold for every computed complex scheme.  A violation
    is a bug, never a discovery.
---------------------------------------------------------------------------
"""

import itertools
from collections import defaultdict, Counter


class NotTypeI(Exception):
    """The curve admits no symmetric orientation."""


class InvariantViolation(Exception):
    """An internal consistency check failed -- indicates a bug, not a result."""


# --------------------------------------------------------------------------
# 1.  Harnack predicate and the symmetric orientation system
# --------------------------------------------------------------------------

def is_harnack_pair(sig, v1, v2, v3, v4):
    """eps restricted to tau u tau' is a Harnack distribution?

    tau, tau' share the edge {v1,v2}; v3, v4 are the opposite vertices.
    sig maps a lattice point to 0/1 with eps = (-1)^sig.

    Primitivity makes v1,v2,v3 pairwise distinct in parity, likewise v1,v2,v4.

      * v3, v4 of different parity -- all four parity classes occur.  A Harnack
        distribution gives one class mu and the other three -mu, so the product
        of the four signs is mu*(-mu)^3 = -1.
      * v3, v4 of the same parity -- only three classes occur, and any
        assignment extends to a Harnack distribution iff eps(v3) = eps(v4).

    Verified equivalent to the 8-case search over (c,a,b) in F2^3 on all 5040
    adjacent-pair configurations arising in cen, spl, fra, hon, 7H2.
    """
    p3 = (v3[0] % 2, v3[1] % 2)
    p4 = (v4[0] % 2, v4[1] % 2)
    if p3 != p4:
        return (sig[v1] + sig[v2] + sig[v3] + sig[v4]) % 2 == 1
    return sig[v3] == sig[v4]


def dual_adjacency(tris):
    """[(i, j, v1, v2, v3, v4)] for each interior edge of Gamma."""
    edge = defaultdict(list)
    for k, t in enumerate(tris):
        for e in itertools.combinations(t, 2):
            edge[tuple(sorted(e))].append(k)
    out = []
    for e, ts in edge.items():
        if len(ts) != 2:
            continue
        i, j = ts
        v3 = next(v for v in tris[i] if v not in e)
        v4 = next(v for v in tris[j] if v not in e)
        out.append((i, j, e[0], e[1], v3, v4))
    return out


def symmetric_orientation(tris, sig):
    """Solve Parenti's F2 system.  Returns [x_tau] or raises NotTypeI.

    x_tau = 0 keeps the counterclockwise reference orientation, 1 reverses it.
    Never returns a fallback: inconsistency is the only route to type II.
    """
    adj = defaultdict(list)
    for (i, j, v1, v2, v3, v4) in dual_adjacency(tris):
        h = 1 if is_harnack_pair(sig, v1, v2, v3, v4) else 0
        w = (h + 1) % 2                      # r == 1 for ccw references
        adj[i].append((j, w))
        adj[j].append((i, w))
    x = {0: 0}
    stack = [0]
    while stack:
        u = stack.pop()
        for v, w in adj[u]:
            if v not in x:
                x[v] = x[u] ^ w
                stack.append(v)
            elif x[v] != x[u] ^ w:
                raise NotTypeI("inconsistent cycle in the Definition 2.7 system")
    if len(x) != len(tris):
        raise InvariantViolation("dual graph of Gamma is disconnected")
    return [x[k] for k in range(len(tris))]


# --------------------------------------------------------------------------
# 2.  The chart
# --------------------------------------------------------------------------

QUADRANTS = [(1, 1), (1, -1), (-1, 1), (-1, -1)]


def _ccw(t):
    a, b, c = t
    cr = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    return (a, b, c) if cr > 0 else (a, c, b)


def extended_sign(sig, p):
    """eps on the diamond: eps(+-i,+-j) = eps(i,j)*(-1)^(i*[x<0] + j*[y<0])."""
    x, y = p
    v = sig[(abs(x), abs(y))]
    if x < 0:
        v += abs(x)
    if y < 0:
        v += abs(y)
    return v & 1


def build_chart(tris, sig, degree, orientations=None):
    """Return (copies, cells).

    copies[ci] = ordered vertex triple of a triangle copy, carrying the
                 reflected symmetric orientation.
    cells[ci]  = (m1, m2, Z, X, Y): one arc directed m1 -> m2 cutting off Z,
                 with m1 on edge (X,Z) and m2 on edge (Y,Z).
    """
    x = symmetric_orientation(tris, sig) if orientations is None else orientations
    copies = []
    for k, t in enumerate(tris):
        base = _ccw(t)
        if x[k]:
            base = (base[0], base[2], base[1])
        for g in QUADRANTS:
            copies.append(tuple((g[0] * v[0], g[1] * v[1]) for v in base))

    cells = {}
    for ci, vs in enumerate(copies):
        s = [extended_sign(sig, v) for v in vs]
        if s[0] == s[1] == s[2]:
            continue
        # empty edge = consecutive pair, in the cyclic order, with equal signs
        for a in range(3):
            b, c = (a + 1) % 3, (a + 2) % 3
            if s[a] == s[b]:
                X, Y, Z = vs[a], vs[b], vs[c]
                break
        m1 = ((X[0] + Z[0]) / 2.0, (X[1] + Z[1]) / 2.0)
        m2 = ((Y[0] + Z[0]) / 2.0, (Y[1] + Z[1]) / 2.0)
        cells[ci] = (m1, m2, Z, X, Y)
    return copies, cells


# --------------------------------------------------------------------------
# 3.  Reading the complex scheme off a chart
# --------------------------------------------------------------------------

def _left_of(m1, m2, q):
    return (m2[0] - m1[0]) * (q[1] - m1[1]) - (m2[1] - m1[1]) * (q[0] - m1[0]) > 0


def read_chart(copies, cells, degree, signed=True):
    """-> dict(tree=..., signs=..., n_ovals=...).  Raises InvariantViolation."""
    D = degree

    def edge_class(p, q):
        e = tuple(sorted((p, q)))
        if abs(p[0]) + abs(p[1]) == D and abs(q[0]) + abs(q[1]) == D:
            f = tuple(sorted(((-p[0], -p[1]), (-q[0], -q[1]))))
            return min(e, f)
        return e

    # --- faces: a cell with an arc splits into the Z-side and the XY-side ---
    face = {}

    def fid(ci, part):
        k = (ci, part)
        if k not in face:
            face[k] = len(face)
        return face[k]

    for ci in range(len(copies)):
        if ci in cells:
            fid(ci, 1); fid(ci, 2)
        else:
            fid(ci, 0)

    def face_at(ci, v):
        if ci not in cells:
            return fid(ci, 0)
        return fid(ci, 1) if v == cells[ci][2] else fid(ci, 2)

    # --- INVARIANT: the symmetric orientation must glue head to tail ---
    indeg, outdeg = Counter(), Counter()
    for ci, (m1, m2, Z, X, Y) in cells.items():
        outdeg[edge_class(Y, Z)] += 1
        indeg[edge_class(X, Z)] += 1
    for e in (set(indeg) | set(outdeg)) if signed else ():
        if indeg[e] != 1 or outdeg[e] != 1:
            raise InvariantViolation(
                "arcs do not glue head-to-tail at an edge midpoint")

    # --- adjacency of faces across diamond edges ---
    owner = defaultdict(list)
    for ci, vs in enumerate(copies):
        for a in range(3):
            owner[edge_class(vs[a], vs[(a + 1) % 3])].append((ci, vs[a], vs[(a + 1) % 3]))

    adjacency = []
    for e, ow in owner.items():
        if len(ow) != 2:
            continue
        (c1, p1, q1), (c2, p2, q2) = ow
        boundary = (abs(p1[0]) + abs(p1[1]) == D and abs(q1[0]) + abs(q1[1]) == D
                    and {p1, q1} != {p2, q2})
        mate = (lambda v: (-v[0], -v[1])) if boundary else (lambda v: v)
        flip = 1 if boundary else 0
        crossed = any(ci in cells and ((pp == cells[ci][2]) != (qq == cells[ci][2]))
                      for (ci, pp, qq) in ow)
        if crossed:
            for v in (p1, q1):
                adjacency.append((face_at(c1, v), face_at(c2, mate(v)), flip))
        else:
            adjacency.append((face_at(c1, p1), face_at(c2, mate(p1)), flip))

    # --- regions of RP^2 minus the curve ---
    n = len(face)
    parent_uf = list(range(n))

    def find(a):
        while parent_uf[a] != a:
            parent_uf[a] = parent_uf[parent_uf[a]]
            a = parent_uf[a]
        return a

    for a, b, _ in adjacency:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent_uf[ra] = rb
    region = {}
    for k in range(n):
        r = find(k)
        if r not in region:
            region[r] = len(region)
    R = lambda p: region[find(p)]

    # --- J: the component whose arcs have the same region on both sides ---
    j_cells = {ci for ci in cells if R(fid(ci, 1)) == R(fid(ci, 2))}
    if not j_cells:
        raise InvariantViolation("no one-sided component found (odd degree?)")
    root = R(fid(next(iter(j_cells)), 1))

    # --- nesting tree (ovals only) ---
    radj = defaultdict(set)
    for ci in cells:
        if ci in j_cells:
            continue
        a, b = R(fid(ci, 1)), R(fid(ci, 2))
        radj[a].add(b); radj[b].add(a)
    parent = {root: None}
    stack, seen = [root], {root}
    while stack:
        u = stack.pop()
        for v in radj[u]:
            if v not in seen:
                seen.add(v); parent[v] = u; stack.append(v)
    if len(seen) != len(region):
        raise InvariantViolation("region adjacency graph is not a spanning tree")

    if not signed:
        children = defaultdict(list)
        for v, u in parent.items():
            if u is not None:
                children[u].append(v)
        return dict(root=root, children=children, signs={},
                    n_ovals=len(region) - 1)

    # --- orientation of RP^2 \ J ---
    o = {}
    graph = defaultdict(list)
    for a, b, f in adjacency:
        graph[a].append((b, f)); graph[b].append((a, f))
    for ci in cells:                       # merge across OVAL arcs only
        if ci in j_cells:
            continue
        graph[fid(ci, 1)].append((fid(ci, 2), 0))
        graph[fid(ci, 2)].append((fid(ci, 1), 0))
    start = fid(next(iter(j_cells)), 1)
    o[start] = 0
    stack = [start]
    while stack:
        u = stack.pop()
        for v, f in graph[u]:
            if v not in o:
                o[v] = o[u] ^ f
                stack.append(v)
            elif o[v] != o[u] ^ f:
                raise InvariantViolation(
                    "orientation of RP^2 minus J is inconsistent (INVARIANT A)")
    if len(o) != n:
        raise InvariantViolation("RP^2 minus J is disconnected")
    for ci in j_cells:
        if o[fid(ci, 1)] == o[fid(ci, 2)]:
            raise InvariantViolation(
                "ambient orientation fails to reverse across J (INVARIANT B)")

    # --- normalise o against J's own complex orientation ---
    ci0 = next(iter(j_cells))
    m1, m2, Z, X, Y = cells[ci0]
    shift = 0 if (_left_of(m1, m2, Z) == (o[fid(ci0, 1)] == 0)) else 1

    # --- oval signs ---
    signs = {}
    for ci in cells:
        if ci in j_cells:
            continue
        m1, m2, Z, X, Y = cells[ci]
        a, b = fid(ci, 1), fid(ci, 2)
        ra, rb = R(a), R(b)
        if parent.get(ra) == rb:
            inner, iface, rep = ra, a, Z
        else:
            inner, iface, rep = rb, b, ((X[0] + Y[0]) / 2.0, (Y[1] + X[1]) / 2.0)
        s = '-' if (_left_of(m1, m2, rep) == ((o[iface] ^ shift) == 0)) else '+'
        if inner in signs and signs[inner] != s:
            raise InvariantViolation(
                "an oval receives two different signs (INVARIANT C)")
        signs[inner] = s

    children = defaultdict(list)
    for v, u in parent.items():
        if u is not None:
            children[u].append(v)
    return dict(root=root, children=children, signs=signs,
                n_ovals=len(region) - 1)


# --------------------------------------------------------------------------
# 4.  Notation and the Rokhlin self-check
# --------------------------------------------------------------------------

def format_scheme(tree, signed=True):
    """Viro recursive notation; signed=False gives the underlying real scheme."""
    root, children, signs = tree['root'], tree['children'], tree['signs']

    def word(u):
        subs = sorted(word(v) for v in children[u])
        if not subs:
            return ""
        c = Counter(subs)
        out = []
        for w in sorted(c, key=lambda z: (len(z), z)):
            k = c[w]
            out.append(str(k) if not w else (f"{k}<{w}>" if k > 1 else f"1<{w}>"))
        return " u ".join(out)

    def node(u):
        inner = word_signed(u) if signed else word(u)
        tag = signs[u] if signed else ""
        return f"1{tag}<{inner}>" if inner else f"1{tag}"

    def word_signed(u):
        subs = sorted(node(v) for v in children[u])
        if not subs:
            return ""
        c = Counter(subs)
        out = []
        for w in sorted(c, key=lambda z: (len(z), z)):
            k = c[w]
            if k == 1:
                out.append(w)
            elif w.startswith("1") and "<" not in w:
                out.append(f"{k}{w[1:]}")
            else:
                out.append(" u ".join([w] * k))
        return " u ".join(out)

    body = word_signed(root) if signed else word(root)
    return "J" + (" u " + body if body else "")


def rokhlin_residual(tree, degree):
    """LHS - RHS of the Rokhlin-Mishachev formula.  Must be 0."""
    root, children, signs = tree['root'], tree['children'], tree['signs']
    ovals = [v for v in signs]
    l = len(ovals)
    Lp = sum(1 for v in ovals if signs[v] == '+')
    Lm = l - Lp
    # all injective pairs: ancestor/descendant in the nesting tree
    depth_stack = [(root, [])]
    Pip = Pim = 0
    while depth_stack:
        u, anc = depth_stack.pop()
        for v in children[u]:
            for a in anc:
                if signs[a] != signs[v]:
                    Pip += 1
                else:
                    Pim += 1
            depth_stack.append((v, anc + [v]))
    k = (degree - 1) // 2
    return 2 * (Pip - Pim) + Lp - Lm - (l - k * (k + 1))


# --------------------------------------------------------------------------
# 5.  Entry point
# --------------------------------------------------------------------------

def complex_scheme(tris, sig, degree=7, check=True):
    """-> (klein_type, real_scheme, complex_scheme).

    klein_type is 'I' or 'II'.  For type II the complex scheme is None.
    Raises InvariantViolation on any internal inconsistency; never silently
    downgrades to type II.
    """
    try:
        copies, cells = build_chart(tris, sig, degree)
        dividing = True
    except NotTypeI:
        # type II: the real scheme is still well defined -- it does not depend
        # on orientations -- so build the chart with the reference orientations.
        copies, cells = build_chart(tris, sig, degree, orientations=[0]*len(tris))
        dividing = False
    tree = read_chart(copies, cells, degree, signed=dividing)
    if not dividing:
        return 'II', format_scheme(tree, signed=False), None
    real = format_scheme(tree, signed=False)
    comp = format_scheme(tree, signed=True)
    if check:
        r = rokhlin_residual(tree, degree)
        if r != 0:
            raise InvariantViolation(
                f"Rokhlin-Mishachev residual {r} != 0 for {comp} (INVARIANT D)")
    return 'I', real, comp
