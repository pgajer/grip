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
logarithmic depth, and fixed tuning parameters. It is not an unconditional
guarantee or an empirical scaling result. The derivations do not cover optional
GMDS refinement.

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
