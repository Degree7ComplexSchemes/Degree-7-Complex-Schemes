"""Exact data audit for an explicit small-parameter nonsingularity argument.

This is NOT a pairwise rigid-isotopy separation verifier. The accompanying
proof note supplies the mathematical argument using GKZ and Sombra.
Only the Python standard library is used; no prior verifier is imported.
"""

import hashlib
import json
import math
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "extension_round4_lower_bound_families.json"
OUT = ROOT / "septic_explicit_nonsingularity_v1.json"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def audit_height(points, entry):
    heights = entry["integer_heights"]
    require(len(heights) == len(points), "height count")
    require(all(type(h) is int for h in heights), "noninteger height")
    require(min(heights) == 0, "height minimum not zero")
    hmap = dict(zip(points, heights))
    triangles = [tuple(tuple(v) for v in tri) for tri in entry["triangles"]]
    require(len(triangles) == 49, "triangle count")
    require(len({tuple(sorted(t)) for t in triangles}) == 49, "repeated triangle")
    valence = Counter()
    edges = Counter()
    facets = []
    gap_count = 0
    min_gap = None
    for tri in triangles:
        require(len(tri) == 3 and len(set(tri)) == 3, "invalid triangle")
        require(all(v in hmap for v in tri), "unknown vertex")
        p, q, r = tri
        dxq, dyq = q[0] - p[0], q[1] - p[1]
        dxr, dyr = r[0] - p[0], r[1] - p[1]
        determinant = dxq * dyr - dyq * dxr
        require(abs(determinant) == 1, "nonprimitive triangle")
        dhq, dhr = hmap[q] - hmap[p], hmap[r] - hmap[p]
        a = (dhq * dyr - dhr * dyq) // determinant
        b = (dxq * dhr - dxr * dhq) // determinant
        c = hmap[p] - a * p[0] - b * p[1]
        require(all(a*x + b*y + c == hmap[x, y] for x, y in tri), "plane error")
        local_min = None
        for v in points:
            if v in tri:
                continue
            gap = hmap[v] - (a * v[0] + b * v[1] + c)
            require(gap > 0, "not a strict lower facet")
            gap_count += 1
            local_min = gap if local_min is None else min(local_min, gap)
            min_gap = gap if min_gap is None else min(min_gap, gap)
        facets.append({"triangle": tri, "plane": [a, b, c], "minimum_gap": local_min})
        valence.update(tri)
        edges.update(tuple(sorted((tri[i], tri[(i + 1) % 3]))) for i in range(3))
    require(set(valence) == set(points), "unused lattice point")
    require(all(n in (1, 2) for n in edges.values()), "edge incidence")
    boundary = [edge for edge, n in edges.items() if n == 1]
    require(len(boundary) == 21 and len(edges) == 84, "triangulated disk counts")
    for p, q in boundary:
        require((p[0] == q[0] == 0) or (p[1] == q[1] == 0)
                or (sum(p) == sum(q) == 7), "interior boundary edge")
    # All triangles lie in Delta_7, are distinct strict lower faces, and
    # have total normalized area 49. Thus they cover Delta_7 without overlap.
    kappa = []
    for x, y in points:
        corner = (x, y) in ((0, 0), (7, 0), (0, 7))
        on_boundary = x == 0 or y == 0 or x + y == 7
        subtract = 1 if corner else (2 if on_boundary else 0)
        kappa.append(valence[x, y] - subtract)
    require(all(k >= 0 for k in kappa), "negative discriminant exponent")
    require(sum(kappa) == 108, "discriminant degree")
    return {"strict_lower_facets": facets, "strict_gap_checks": gap_count,
            "minimum_gap": min_gap, "discriminant_initial_exponents": kappa,
            "initial_t_exponent": sum(k*h for k, h in zip(kappa, heights)),
            "maximum_height": max(heights)}
