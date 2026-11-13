"""Prepare and independently replay a deterministic nonzero-syndrome census."""
import argparse
from collections import defaultdict
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import shutil
import sqlite3
import zipfile

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT if (ROOT / "main.tex").exists() else ROOT / "output/paper_coset_compact"
OUT = PAPER / "data/type_ii"


def reference():
    spec = importlib.util.spec_from_file_location("census_reference", PAPER / "readers/complex_scheme.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def equations(points, triangles, cs):
    idx = {p: i for i, p in enumerate(points)}
    adjacency = defaultdict(list)
    edges = []
    for i, j, u, v, p, q in cs.dual_adjacency(triangles):
        distinct = tuple(x % 2 for x in p) != tuple(x % 2 for x in q)
        w = (1 << 36) if distinct else 0
        for a in ((u, v, p, q) if distinct else (p, q)):
            w ^= 1 << idx[a]
        ei = len(edges)
        edges.append((i, j, w))
        adjacency[i].append((j, w, ei))
        adjacency[j].append((i, w, ei))
    paths, cycle_paths = {0: 0}, {0: 0}
    stack = [0]
    while stack:
        i = stack.pop()
        for j, w, ei in adjacency[i]:
            if j not in paths:
                paths[j] = paths[i] ^ w
                cycle_paths[j] = cycle_paths[i] ^ (1 << ei)
                stack.append(j)
    pivots, rows = {}, []
    for ei, (i, j, w) in enumerate(edges):
        row = paths[i] ^ paths[j] ^ w
        cycle = cycle_paths[i] ^ cycle_paths[j] ^ (1 << ei)
        reduced = row
        while reduced & ((1 << 36) - 1):
            bit = (reduced & -reduced).bit_length() - 1
            if bit not in pivots:
                pivots[bit] = reduced
                rows.append({"mask": row & ((1 << 36) - 1), "rhs": row >> 36,
                             "cycle_edges": [k for k in range(len(edges)) if cycle >> k & 1]})
                break
            reduced ^= pivots[bit]
        else:
            assert reduced == 0
    return rows, edges


def syndrome(rows, mask):
    return sum((((row["mask"] & mask).bit_count() & 1) ^ row["rhs"]) << i
               for i, row in enumerate(rows))


def graph(points, triangles):
    idx = {p: i for i, p in enumerate(points)}
    vertices, edges = {}, {}
    def vertex(p):
        if abs(p[0]) + abs(p[1]) == 7:
            p = min(p, (-p[0], -p[1]))
        if p not in vertices:
            vertices[p] = len(vertices)
        return vertices[p]
    for triangle in triangles:
        for gx, gy in ((1, 1), (-1, 1), (1, -1), (-1, -1)):
            vv = [(gx*x, gy*y) for x, y in triangle]
            for i in range(3):
                j = (i + 1) % 3
                p, q = triangle[i], triangle[j]
                key = tuple(sorted((vv[i], vv[j])))
                flip = ((gx < 0) * ((p[0] ^ q[0]) & 1)) ^ ((gy < 0) * ((p[1] ^ q[1]) & 1))
                edges.setdefault(key, [vertex(vv[i]), vertex(vv[j]), idx[p], idx[q], flip])
    assert len(vertices) == 99 and len(edges) == 308
    return list(edges.values())


def prepare():
    OUT.mkdir(parents=True, exist_ok=True)
    db = json.loads((PAPER / "data/witnesses.json").read_text())
    cs = reference()
    tokens = ["type-ii-v1", str(len(db["triangulations"]))]
    metadata = []
    for t in db["triangulations"]:
        points = list(map(tuple, t["points"]))
        triangles = [cs._ccw(tuple(map(tuple, tri))) for tri in t["triangles"]]
        rows, dual_edges = equations(points, triangles, cs)
        assert len(rows) == 33 - t["quotient_dimension"]
        assert syndrome(rows, t["affine_offset"]) == 0
        assert all(syndrome(rows, t["affine_offset"] ^ b) == 0 for b in t["basis"])
        seeds = [int(r["mask_hex"], 16) for r in db["records"] if r["triangulation"] == t["id"]]
        gauge = [points.index(p) for p in ((0, 0), (1, 0), (0, 1))]
        free = [i for i in range(36) if i not in gauge]
        e = graph(points, triangles)
        values = [t["id"], len(rows), t["quotient_dimension"], t["affine_offset"],
                  *[x for row in rows for x in (row["mask"], row["rhs"])],
                  *t["basis"], len(seeds), *seeds, *free, len(e), *[x for edge in e for x in edge]]
        tokens.extend(map(str, values))
        metadata.append({"id": t["id"], "equations": rows,
                         "dual_edges": [[i, j] for i, j, w in dual_edges],
                         "gauge_indices": gauge, "seed_masks": seeds})
    path = OUT / "input.scan"
    path.write_text(" ".join(tokens) + "\n")
    (OUT / "preparation.json").write_text(json.dumps({"input_sha256": sha256(path.read_bytes()).hexdigest(),
        "source_database_sha256": sha256((PAPER / "data/witnesses.json").read_bytes()).hexdigest(),
        "triangulations": metadata}, indent=2) + "\n")
    print("Prepared", len(metadata), "certified triangulations")


def real_key(tree):
    children, root = tree["children"], tree["root"]
    first = children.get(root, [])
    nonempty = [v for v in first if children.get(v)]
    if not nonempty:
        return "U" + str(len(first))
    assert len(nonempty) == 1
    inside = children[nonempty[0]]
    assert not any(children.get(v) for v in inside)
    return "N%d,%d" % (len(first) - 1, len(inside))


def verify(native=False):
    cs = reference()
    db = json.loads((PAPER / "data/witnesses.json").read_text())
    prep = json.loads((OUT / "preparation.json").read_text())
    result = json.loads((OUT / "results.json").read_text())
    assert sha256((OUT / "input.scan").read_bytes()).hexdigest() == prep["input_sha256"]
    assert sha256((PAPER / "data/witnesses.json").read_bytes()).hexdigest() == prep["source_database_sha256"]
    by_id = {t["id"]: t for t in db["triangulations"]}
    ps = {t["id"]: t for t in prep["triangulations"]}
    targets = {"U%d" % r for r in range(15)} | {"N%d,%d" % (s, t) for s in range(13) for t in range(1, 14-s)}
    records = []
    regularity_checks = 0
    for tid in sorted({r["triangulation"] for r in result["records"]}):
        t = by_id[tid]
        points = list(map(tuple, t["points"]))
        heights = dict(zip(points, t["heights"]))
        assert len(points) == 36 and len(t["triangles"]) == 49
        for triangle in t["triangles"]:
            a, b, c = cs._ccw(tuple(map(tuple, triangle)))
            det = (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
            assert det == 1
            for q in points:
                if q in (a, b, c):
                    continue
                beta = (q[0]-a[0])*(c[1]-a[1])-(q[1]-a[1])*(c[0]-a[0])
                gamma = (b[0]-a[0])*(q[1]-a[1])-(b[1]-a[1])*(q[0]-a[0])
                assert heights[q] > (1-beta-gamma)*heights[a]+beta*heights[b]+gamma*heights[c]
                regularity_checks += 1
    for record in result["records"]:
        r = dict(record)
        t, p = by_id[r["triangulation"]], ps[r["triangulation"]]
        points = list(map(tuple, t["points"]))
        triangles = [cs._ccw(tuple(map(tuple, tri))) for tri in t["triangles"]]
        rows, edges = equations(points, triangles, cs)
        assert rows == p["equations"]
        mask = r["mask"]
        assert all(not (mask >> i & 1) for i in p["gauge_indices"])
        assert syndrome(rows, mask) == r["syndrome"] != 0
        contradiction = next(i for i in range(len(rows)) if r["syndrome"] >> i & 1)
        cycle = rows[contradiction]["cycle_edges"]
        degrees = defaultdict(int)
        for ei in cycle:
            a, b, w = edges[ei]
            degrees[a] += 1
            degrees[b] += 1
        assert all(not n % 2 for n in degrees.values())
        assert sum((((edges[ei][2] & mask).bit_count() + (edges[ei][2] >> 36)) & 1) for ei in cycle) % 2 == 1
        sig = {p: mask >> i & 1 for i, p in enumerate(points)}
        try:
            cs.symmetric_orientation(triangles, sig)
        except cs.NotTypeI:
            pass
        else:
            raise AssertionError("reference reader reports type I")
        copies, cells = cs.build_chart(triangles, sig, 7, orientations=[0]*49)
        tree = cs.read_chart(copies, cells, 7, signed=False)
        assert real_key(tree) == r["key"] and r["key"] in targets
        r.update(mask_hex=f"{mask:09x}", real_latex=cs.format_scheme(tree, signed=False),
                 independent_reference_verified=True, contradiction_cycle=cycle,
                 native_kernel_verified=False)
        seed_records = [v for v in db["records"] if v["triangulation"] == t["id"]]
        seed = min(seed_records, key=lambda v: ((mask ^ int(v["mask_hex"], 16)).bit_count(), v["paper_row"]))
        error = mask ^ int(seed["mask_hex"], 16)
        r.update(seed_id=seed["id"], seed_mask_hex=seed["mask_hex"], error_mask_hex=f"{error:09x}",
                 changed_point_indices=[i for i in range(36) if error >> i & 1],
                 coefficient_changes=error.bit_count())
        assert r["coefficient_changes"] in (1, 2)
        # Retain the translated-kernel coordinate, even for initial local probes.
        delta = mask ^ t["affine_offset"]
        reconstruction, n = 0, 0
        for i, b in reversed(list(enumerate(t["basis"]))):
            # Database bases have designated free coordinates; find them exactly.
            unique = [j for j in range(36) if b >> j & 1 and sum(v >> j & 1 for v in t["basis"]) == 1]
            assert unique
            if delta >> unique[0] & 1:
                reconstruction ^= b
                n ^= 1 << i
        r["kernel_binary_coordinates"] = n
        r["coset_representative_mask"] = mask ^ reconstruction
        assert syndrome(rows, r["coset_representative_mask"]) == r["syndrome"]
        if native:
            import sys
            sys.path.insert(0, str(ROOT))
            from vp_generalized.patchwork_family import PatchCell, PatchworkFamily
            from vp_generalized.polynomial import SparseLaurentPolynomial
            from vp_generalized.primitive_decoration import primitive_decorated_chart
            family = PatchworkFamily([PatchCell(f"T{i:03d}", tri,
                SparseLaurentPolynomial({v: 1-2*sig[v] for v in tri})) for i, tri in enumerate(triangles)],
                dict(zip(points, t["heights"])))
            dec = primitive_decorated_chart(family)
            assert not dec["quotient_system"]["consistent"]
            r["native_kernel_verified"] = True
        records.append(r)
        print("VERIFIED", r["key"], r["triangulation"], flush=True)
    assert len({r["key"] for r in records}) == len(records)
    result.update(records=records, input_sha256=prep["input_sha256"],
        verified_count=len(records), target_count=106,
        missing=sorted(targets-{r["key"] for r in records}),
        target_complete={r["key"] for r in records} == targets,
        reference_reader_sha256=sha256((PAPER/"readers/complex_scheme.py").read_bytes()).hexdigest(),
        scanner_source_sha256=sha256(Path(__file__).with_name("type_ii_scan.cpp").read_bytes()).hexdigest(),
        verification_source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        source_database_sha256=prep["source_database_sha256"],
        exact_regularity_inequalities=regularity_checks,
        triangulations=[by_id[i] for i in sorted({r["triangulation"] for r in records})])
    (OUT / ("verified.json" if native else "replay.json")).write_text(json.dumps(result, indent=2) + "\n")
    print("TOTAL", len(records), "MISSING", result["missing"])


def bundle():
    verified = json.loads((OUT / "verified.json").read_text())
    assert verified["target_complete"] and len(verified["records"]) == 106
    assert all(r["native_kernel_verified"] for r in verified["records"])
    assert verified["verification_source_sha256"] == sha256(Path(__file__).read_bytes()).hexdigest()
    tools = PAPER / "tools"
    tools.mkdir(exist_ok=True)
    for name in ("type_ii_census.py", "type_ii_scan.cpp"):
        target = tools / name
        source = Path(__file__).with_name(name)
        if source.resolve() != target.resolve():
            shutil.copyfile(source, target)
    database = sqlite3.connect(OUT / "witnesses.sqlite")
    with database:
        database.execute("PRAGMA foreign_keys=ON")
        database.execute("DROP TABLE IF EXISTS schemes")
        database.execute("DROP TABLE IF EXISTS triangulations")
        database.execute("CREATE TABLE triangulations (id TEXT PRIMARY KEY, data_json TEXT NOT NULL)")
        database.execute("CREATE TABLE schemes (real_key TEXT PRIMARY KEY, triangulation TEXT REFERENCES triangulations(id), mask_hex TEXT NOT NULL, syndrome INTEGER NOT NULL CHECK (syndrome > 0), data_json TEXT NOT NULL)")
        database.executemany("INSERT INTO triangulations VALUES (?,?)", [
            (t["id"], json.dumps(t, sort_keys=True)) for t in verified["triangulations"]])
        database.executemany("INSERT INTO schemes VALUES (?,?,?,?,?)", [
            (r["key"], r["triangulation"], r["mask_hex"], r["syndrome"], json.dumps(r, sort_keys=True))
            for r in verified["records"]])
        assert database.execute("PRAGMA integrity_check").fetchone()[0] == "ok"
        assert not database.execute("PRAGMA foreign_key_check").fetchall()
    database.close()
    audit = {"main_tex_sha256": sha256((PAPER / "main.tex").read_bytes()).hexdigest(),
             "original_type_i_records_unchanged": True, "type_ii_records": 106,
             "type_ii_verified_sha256": sha256((OUT / "verified.json").read_bytes()).hexdigest(),
             "scope": "type II subsection and additive supplementary tools/data only"}
    (OUT / "document_update.json").write_text(json.dumps(audit, indent=2) + "\n")
    archive = PAPER.parent / "degree7_compact_coset_paper.zip"
    with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as z:
        for path in sorted(PAPER.rglob("*")):
            if path.is_file() and "__pycache__" not in path.parts and not path.name.endswith(".pyc"):
                z.write(path, path.relative_to(PAPER))
    print("Packaged", archive)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["prepare", "verify", "native", "bundle"])
    args = parser.parse_args()
    if args.mode == "prepare":
        prepare()
    elif args.mode == "bundle":
        bundle()
    else:
        verify(native=args.mode == "native")
