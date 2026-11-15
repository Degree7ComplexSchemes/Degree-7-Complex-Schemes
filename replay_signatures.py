"""Verify the exact ordinary and colored signature calculations in Section 7.6."""

from fractions import Fraction
import json
from pathlib import Path
import sys

from regional_form import exact, integer_inertia, model, pencil_braid, tau
from seifert import realification, seifert_matrix


def exponent(word):
    return sum(1 if g > 0 else -1 for g in word)


def exterior_cut():
    cases = []
    for phase in (-1, 1):
        for positive_slot in range(7):
            caps = [3 if (1 if j == positive_slot else -1) == phase * (-1)**j
                    else 4 for j in range(7)]
            path = caps + [7 - caps[0]]
            braid = [g for k, l in zip(path, path[1:]) for g in [-k] + tau(k, l)]
            for last in range(6, 0, -1):
                braid.extend(range(1, last + 1))
            assert exponent(braid) == 14
            v = seifert_matrix(braid, 7)
            form = [[v[i][j] + v[j][i] for j in range(len(v))] for i in range(len(v))]
            inertia, pivots = integer_inertia(form)
            assert inertia == (11, 21, 0)
            cases.append(dict(phase=phase, positive_slot=positive_slot, caps=path,
                              braid=braid, seifert_matrix=v, inertia=inertia, pivots=pivots,
                              signature=-10, adjusted_nullity=1))
    return cases


def colored(bands, palette, q, r, cycles_color, links, expected_inertia):
    braid = pencil_braid(list(map(int, bands)))
    assert exponent(braid) == 12
    regions = model(braid, palette)
    assert regions[3] == [[1, 7], [2], [3], [4, 6], [5]]
    assert regions[5] == cycles_color
    linking = regions[4]
    assert all(linking[i][j] == value for i, j, value in links)
    result = exact(regions, q, r)
    assert result["inertia"] == expected_inertia and result["writhe"] == 12
    assert (result["sigma"], result["eta"]) == (-8, 0)
    return dict(bands=bands, palette=palette, braid=braid, closure_cycles=regions[3],
                linking=linking, **result)


def three_inner():
    bands = "243434"
    braid = pencil_braid(list(map(int, bands)))
    assert exponent(braid) == 14
    v = seifert_matrix(braid, 7)
    inertia, pivots = integer_inertia(realification(v, 8))
    assert inertia == (26, 46, 0)
    return dict(bands=bands, braid=braid, seifert_matrix=v, realified_inertia=inertia,
                pivots=pivots, signature=-10, adjusted_nullity=1)


def main():
    assert __debug__, "Run with assertions enabled."
    data = {
        "7.6 N_I+(8;7,6)": exterior_cut(),
        "7.7 N_I-(10;3,3)": colored("34322222", [1, 1, 0, 1, 1, 1, 1],
                                   Fraction(1, 8), Fraction(1), [1, 1, 0, 1, 1],
                                   [(0, 1, 1), (0, 3, 2), (0, 4, 1)], (20, 24, 0)),
        "7.8 N_I+(8;5,5)": three_inner(),
        "7.9 N_I+(10;3,3)": colored("22224342", [0, 0, 1, 0, 1, 0, 0],
                                   Fraction(1, 20), Fraction(19, 20), [0, 0, 1, 0, 1],
                                   [(0, 1, 1), (0, 3, 2), (2, 4, 3)], (22, 26, 0)),
    }
    if len(sys.argv) == 2:
        with Path(sys.argv[1]).open("x", encoding="utf-8", newline="\n") as handle:
            json.dump(data, handle, indent=2)
            handle.write("\n")
    else:
        expected = Path(__file__).with_name("signature_certificates.json")
        assert json.loads(json.dumps(data)) == json.loads(expected.read_bytes())
    print("PASS: the exact signature calculations of Theorems 7.6–7.9.")


if __name__ == "__main__":
    main()
