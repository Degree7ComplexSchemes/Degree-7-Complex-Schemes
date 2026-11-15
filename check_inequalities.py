"""Replay and check the signature bounds in Theorems 7.6–7.9."""

import json
from pathlib import Path
import sys

from replay_signatures import exponent, main as replay


def check_violation(label, left, right):
    if left <= right:
        raise SystemExit(f"{label}: {left} <= {right}; no obstruction")
    print(f"{label}: {left} > {right}")


def check_ordinary(theorem, cases):
    # Equation (9). The pencil has seven strands and n_zeta is adjusted nullity.
    for case in cases:
        left = abs(case["signature"]) - case["adjusted_nullity"]
        right = exponent(case["braid"]) - 7
        if "caps" in case:
            word = "".join(map(str, case["caps"]))
            label = (f"{theorem} caps {word}, phase {case['phase']:+d}, "
                     f"positive slot {case['positive_slot']}")
        else:
            label = f"{theorem} word {case['bands']}"
        check_violation(label, left, right)


def check_colored(theorem, case):
    # Theorem 7.2 of Cimasoni–Florens. The manuscript proves that the two
    # color surfaces are connected and disjoint, with total chi = -5.
    colors = len(set(case["palette"]))
    total_chi = -5
    beta_1 = colors - total_chi
    left = abs(case["sigma"]) + abs(case["eta"] - colors + 1)
    check_violation(f"{theorem} word {case['bands']}", left, beta_1)


def main():
    if len(sys.argv) != 1:
        raise SystemExit("Usage: python -B check_inequalities.py")
    if not __debug__:
        raise SystemExit("Run with Python assertions enabled.")
    replay()
    path = Path(__file__).with_name("signature_certificates.json")
    data = json.loads(path.read_bytes())
    check_ordinary("7.6", data["7.6 N_I+(8;7,6)"])
    check_colored("7.7", data["7.7 N_I-(10;3,3)"])
    check_ordinary("7.8", [data["7.8 N_I+(8;5,5)"]])
    check_colored("7.9", data["7.9 N_I+(10;3,3)"])


if __name__ == "__main__":
    main()
