# Degree 7 Complex Schemes

Computational companion to Brennan Williams, *Degree 7 Complex Schemes*.
These files accompany the manuscript with Theorem 9.1, which gives the lower
bound of **1,557** rigid isotopy classes. Upload the files to
[Degree-7-Complex-Schemes](https://github.com/Degree7ComplexSchemes/Degree-7-Complex-Schemes)
using the relative paths shown here. [File hashes](SHA256SUMS.txt)
identify the distributed inputs, code, and reports.

## Paper tables

- [Representatives by scheme](data/representatives.csv) follows Tables 2–4
  in print order. It gives the patchworked and published contributions
  separately. Table 2 contributes 258 representatives across 93 schemes,
  Table 3 contributes 1,295 across 106 schemes, and Table 4 contributes 4
  across 3 schemes. The total is 1,549 patchworked families and 8 published
  additions.
- [Type I framework census](data/type_i_framework.csv) follows Table 8.
  It includes 96 algebraically realized codes, 168 excluded codes, and one
  algebraically unresolved code, `N_I+(12;3,2)`.
- The [exclusion crosswalk](data/excluded_candidates.csv) and its
  [LaTeX table](tables/excluded_candidates.tex) retain the old row numbers
  for the 167 codes still present, add `N_I-(6;4,1)` as row 168, and
  point to the revised case-specific theorems. Row 166 now points to
  Theorem 7.9.

## Exact replays

Run from the repository root with Python 3.10 or later and assertions enabled.
These checks use only the standard library and require no network access.

```sh
python -B census/replay_database.py
python -B scripts/replay_type_ii.py --output-dir results/type-ii
python -B signatures/replay_signatures.py
python -B signatures/verify_two_inner_one_band.py
python -B scripts/export_paper_bindings.py
python -B scripts/replay_separation.py
```

The witness scripts read the retained [93 dividing witnesses](census/data/witnesses.json)
and [106 nondividing witnesses](census/data/type_ii/verified.json). They
reconstruct the signed schemes or check the nonzero syndromes and nesting.
The type II replay uses saved scan results.

The [signature certificates](signatures/signature_certificates.json)
cover the fourteen exterior-cut words in Theorem 7.6 and the words
`34322222`, `243434`, and `22224342` in Theorems 7.7–7.9.
The band-diagram script verifies Theorem 7.10's eleven planar moves
and the resulting 18 positive crossings and six Seifert circles.
These are finite algebra and diagram checks. The geometric reductions
and surface arguments are given in the paper.

## Section 9 family data

The [lattice order](data/lattice_points.csv), [height presentations](data/heights.csv),
and [family masks](data/families.csv) give all 1,549 polynomials. The original
1,430-family bank from Table 7 uses IDs 1–1430 and height indices 1–103.
The expanded bank uses 117 height presentations in the audit's convention,
which includes the ordering of the listed triangles. Two presentations
have the same ordered height vector and list the same triangles in different
orders. The replay checks strict lower faces and uses the recorded triangle
order when reading the marked geometry.

The 36-bit family string gives the coefficients in lattice order, with
0 meaning +1 and 1 meaning −1. The first point is the most significant
hexadecimal bit in `data/families.csv` and the least significant bit in
the separate census masks.

The [marked models](data/marked_models.json) give exact Laurent-monomial
coordinates, region labels, lift signs, and dominance gaps for the expanded
bank. The
[expanded audit](expected/whole_bank_audit.json) is the original record
from the companion archive cited in the paper. Its source hashes document
that archive's provenance. The portable replay here reads the consolidated
inputs and checks the reported geometry, signs, and all 24,905 pairs in the
199 same-scheme groups. It takes several minutes. The eight published
additions are accounted for in the representative table. The
[jump comparison](data/jump_comparison.csv) identifies the patchworked
families tested against the published jump constructions.

The two marked-point files cited in Section 9 apply to the original
Table 7 bank. The export command regenerates and compares its
[13,835 marked points](expected/paper-bindings/marked_points.csv) and
[21 scheme-only records](expected/paper-bindings/coarse_only.json).
The geometry that turns these finite computations into a rigid-isotopy
lower bound is supplied by the paper and its cited theorems.

## Curve plates

- [93 dividing curves with triangulations](plates/pdf/degree7-93-dividing-T-curves-with-triangulations.pdf)
- [106 nondividing curves with triangulations](plates/pdf/degree7-106-nondividing-T-curves-with-triangulations.pdf)
- [Plate index](plates/plate_index.json) and [rendering provenance](plates/PROVENANCE.json)

The plate PDFs are ready to view. For regeneration with SageMath and the
bundled Viro.sage snapshot, use the following commands. Assembly uses the
optional Python packages in `plates/requirements.txt`.

```sh
sage -python plates/render_sage.py
python -m pip install -r plates/requirements.txt
python -B plates/build_plates.py --panel-dir plates/generated/panels
```
