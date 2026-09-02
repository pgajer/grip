# Supplement S1: Search-sensitive complexity of weighted GRIP

Supplement to *grip: Multiscale Graph Layout and Geodesic Embedding Tools in R*,
by Pawel Gajer and Jacques Ravel.

The canonical source is `S1-weighted-grip-complexity.tex`, accompanied by its
bibliography and citation-verification evidence. It analyzes the weighted
2D/3D core of grip 0.2.0. The five implementation files identified in the
supplement are unchanged between release tag `v0.2.0` and the pinned commit
`cebcc627179a09564f446982d783fb84ea836ea0` (verified with `git diff`).

The conditional near-linear bound requires bounded degree, bounded overlap
of filtration search balls, the stated coarse-vertex coverage condition,
bounded adjacent level-size ratios, logarithmic depth, and fixed tuning
parameters. It is not an unconditional
guarantee or an empirical scaling result. The derivations do not cover optional
GMDS refinement.

The insertion-search term is sharpened by a telescoping level-size sum to
`O((b + 1 + a gamma) n log(n/q) log n)`, without a separate depth factor.
A compressed-block weighted path shows why bounded degree and doubling do
not imply the coarse-vertex coverage assumption: uncapped searches for three
higher-level anchors can require quadratic total insertion work. The
expected-repulsion calculation explicitly assumes independent uniform draws
and counts rejection of the vertex itself and duplicate partners.

`check-complexity-arguments.py` provides independent finite checks using only
the Python standard library. It checks the telescoping inequality on 4,095
nested size sequences, verifies 36 expected draw counts using exact rational
finite-state recurrences, and simulates the compressed-path stopping rules
for three graph sizes and three greedy center orders. It uses integer-scaled
edge lengths to avoid floating-point ambiguity, includes the extra eligible
vertex required for cache completion, and requires three higher-level anchors.
It models the first coarse level and level-zero insertions only, not the full
package algorithm. These are mathematical checks, not timing benchmarks.

Run `make argument-check` from this directory (or invoke
`python3 check-complexity-arguments.py` directly). The standard supplement
build also runs these checks.

From the repository root, run `make paper-supplement`. The PDF and intermediates
are written to `papers/grip-software-paper/build/supplement/`. To build directly
from this directory, run `make`; LaTeX, latexmk, BibTeX, and Python 3 are needed.
There is no visible internal build stamp in the publication PDF.

In the submission archive, the PDF is supplied alongside the sources and a
portable copy of the citation checker. To rebuild after extraction, run
`make BUILD_DIR=build` from the archive's `supplement/` directory.

The related `reproducibility/` directory contains data, benchmark artifacts,
and computational reproduction scripts; this supplement is a mathematical
analysis and does not replace those materials.
