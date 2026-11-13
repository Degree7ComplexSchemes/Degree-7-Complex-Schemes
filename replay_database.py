from fractions import Fraction
from hashlib import sha256
import importlib.util
import json
from pathlib import Path


def det(a,b,c):
    return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])


def canonical(nodes):
    return tuple(sorted((node["sign"],canonical(node["children"])) for node in nodes))


def tree_canonical(tree):
    def visit(v):
        return tuple(sorted((tree["signs"][w],visit(w)) for w in tree["children"][v]))
    return visit(tree["root"])


def main():
    root=Path(__file__).resolve().parent
    db=json.loads((root/"data/witnesses.json").read_text(encoding="utf-8"))
    reader=root/"readers/complex_scheme.py"
    assert sha256(reader.read_bytes()).hexdigest()==db["reference_reader_sha256"]
    spec=importlib.util.spec_from_file_location("paper_complex_reader",reader)
    cs=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cs)
    triangles={t["id"]:t for t in db["triangulations"]}
    total_checks=0
    for t in triangles.values():
        points=list(map(tuple,t["points"]))
        assert points==[(i,j) for i in range(8) for j in range(8-i)]
        assert len(t["triangles"])==49
        heights=dict(zip(points,t["heights"]))
        geometry={k:t[k] for k in ("points","triangles","heights")}
        assert sha256(json.dumps(geometry,sort_keys=True,separators=(",",":")).encode()).hexdigest()==t["sha256"]
        for tri in t["triangles"]:
            a,b,c=cs._ccw(tuple(map(tuple,tri)))
            assert det(a,b,c)==1
            for q in points:
                if q in (a,b,c):
                    continue
                beta,gamma=det(a,q,c),det(a,b,q)
                assert heights[q]>(1-beta-gamma)*heights[a]+beta*heights[b]+gamma*heights[c]
                total_checks+=1
    seen=set()
    receipts={}
    archived=json.loads((root/"data/type_i_receipts.json").read_bytes())
    for r in db["records"]:
        t=triangles[r["triangulation"]]
        mask=int(r["mask_hex"],16)
        gray=r["gray_index"]^(r["gray_index"]>>1)
        reproduced=t["affine_offset"]
        for i,b in enumerate(t["basis"]):
            if (gray>>i)&1:
                reproduced^=b
        assert reproduced==mask and r["gray_index"]<t["normalized_classes"]
        points=list(map(tuple,t["points"]))
        assert all(not((mask>>points.index(p))&1) for p in ((0,0),(1,0),(0,1)))
        copies,raw=cs.build_chart([tuple(map(tuple,x)) for x in t["triangles"]],
                                  {p:(mask>>i)&1 for i,p in enumerate(points)},7)
        arcs={i:tuple(tuple(Fraction(x) for x in p) for p in a) for i,a in raw.items()}
        tree=cs.read_chart(copies,arcs,7)
        actual=tree_canonical(tree)
        assert actual==canonical(r["complex_scheme"]["ovals"])
        assert actual not in seen
        seen.add(actual)
        assert cs.rokhlin_residual(tree,7)==0
        assert r["native_kernel_replay_verified"] and r["exhaustive_space_membership"]
        path=r["portable_receipt"]
        receipt_bytes=archived[path].encode("utf-8")
        assert sha256(receipt_bytes).hexdigest()==r["receipt_sha256"]
        if str(path) not in receipts:
            receipt=json.loads(receipt_bytes)
            assert receipt["examined"]==1<<receipt["dimension"]
            assert receipt.get("complete") is not False
            assert sum(receipt["component_histogram"].values())==receipt["examined"]
            receipts[str(path)]=receipt
        receipt=receipts[str(path)]
        compact={k:t[k] for k in ("name","affine_offset","basis","triangles")}
        space_hash=sha256(json.dumps(compact,separators=(",",":")).encode()).hexdigest()
        assert space_hash==receipt["space_sha256"]
        assert receipt["affine_offset"]==t["affine_offset"]
        assert receipt["dimension"]==t["quotient_dimension"]
        print(r["id"],cs.format_scheme(tree))
    print(f"PASS: {len(seen)} distinct signed schemes; {len(triangles)} strict regular primitive triangulations;")
    print(f"{total_checks} exact plane inequalities; {len(receipts)} completed-space receipts.")
    print("Historical discovery and later exhaustive recovery remain separately recorded.")


if __name__=="__main__":
    main()
