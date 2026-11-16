"""Replay the finite separation of all 1,549 patchworked families in Section 9."""

from collections import Counter, defaultdict
import csv
from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path
import sys

from septic_explicit_nonsingularity_v1 import audit_height
from septic_three_obligations_v1 import isomorphic, sphere_geometry
from septic_verify_joint_nest_chord_v1 import joint

ROOT = Path(__file__).resolve().parent


def read_inputs():
    data = ROOT
    with (data / "lattice_points.csv").open(newline="", encoding="utf-8") as handle:
        points = [(int(row["i"]), int(row["j"])) for row in csv.DictReader(handle)]
    with (data / "heights.csv").open(newline="", encoding="utf-8") as handle:
        heights = {
            str(row["height_index"]): [int(row[f"h_{i:02}"]) for i in range(36)]
            for row in csv.DictReader(handle)
        }
    with (data / "families.csv").open(newline="", encoding="utf-8") as handle:
        families = list(csv.DictReader(handle))
    with (data / "marked_models.json").open(encoding="utf-8") as handle:
        models = json.load(handle)
    assert points == list(map(tuple, models["points"]))
    assert len(heights) == len(models["triangulations"]) == 117
    return points, heights, families, models


def main():
    assert __debug__, "Run with assertions enabled."
    points, heights, families, models = read_inputs()
    report = json.loads((ROOT / "whole_bank_audit.json").read_bytes())
    assert len(families) == len(models["models"]) == len(report["receipts"]) == 1549
    assert [int(row["family_id"]) for row in families] == list(range(1, 1550))
    for index, values in heights.items():
        audit_height(points, {
            "triangles": models["triangulations"][index], "integer_heights": values,
        })
    groups = defaultdict(list)
    counts = Counter()
    stats = {"nodes": 0}
    for family, data, saved in zip(families, models["models"], report["receipts"], strict=True):
        ident, kind, scheme = data["id"], data["type"], data["scheme"]
        assert ident == int(family["family_id"]) == saved["id"]
        assert kind == family["klein_type"] == saved["type"]
        assert scheme == saved["scheme"]
        index = family["height_index"]
        witness = data["witness"]
        witness["phase_mask"] = family["phase_bits_p0_first"]
        assert len(witness["phase_mask"]) == 36
        assert int(witness["phase_mask"], 2) == int(family["phase_hex_msb_p0"], 16)
        if witness["triangles"] is None:
            witness["triangles"] = models["triangulations"][index]
        assert sorted(map(sorted, witness["triangles"])) == sorted(map(sorted, models["triangulations"][index]))
        old = dict((tuple(p), Fraction(h)) for p, h in witness["heights"])
        values = heights[index]
        j = next(k for k, h in enumerate(values) if h != values[0])
        scale = (old[points[j]] - old[points[0]]) / (values[j] - values[0])
        assert scale > 0
        assert all(old[p] - old[points[0]] == scale * (h - values[0]) for p, h in zip(points, values, strict=True))
        model = data["model"]
        geometry = sphere_geometry(witness, model, kind, scheme)
        assert geometry == saved["sphere_geometry"]
        relation = None
        counts[f"type_{kind}"] += 1
        if model is not None:
            _, relation, signs, conics_c = joint(model)
            assert len(relation.triples) == saved["protected_triples"]
            assert len(relation.conics) == saved["joint_conics"]
            counts.update(signs)
            counts["C_determinants"] += len(conics_c)
            counts["point_bearing_models"] += 1
        groups[json.dumps([kind, scheme])].append((ident, relation))
        if ident % 100 == 0:
            print(f"Checked geometry and signs for {ident} families", flush=True)
    tested = set()
    for records in groups.values():
        for (i, a), (j, b) in combinations(records, 2):
            assert a is not None and b is not None, "Repeated scheme-only record"
            assert isomorphic(a, b, stats) is None, (i, j)
            tested.add((i, j))
    assert tested == {(r["a"], r["b"]) for r in report["pairwise_checks"] if not r["isomorphic"]}
    assert (len(groups), len(tested)) == (199, 24905)
    assert all(counts[name] == report["summary"][name] for name in
               ("type_I", "type_II", "point_bearing_models", "triples", "conics", "C_determinants"))
    print("PASS: 1,549 families, 117 height presentations, 199 groups, 24,905 pairs, zero collisions.")


if __name__ == "__main__":
    if not __debug__:
        sys.exit("Do not use Python -O.")
    main()
