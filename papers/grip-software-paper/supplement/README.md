# Supplementary materials

Supplement to *grip: Multiscale Graph Layout and Geodesic Embedding Tools in R*,
by Pawel Gajer and Jacques Ravel.

## Software versions and provenance

The article describes [grip 0.2.0](https://CRAN.R-project.org/package=grip)
and uses [dgraphs 0.2.0](https://CRAN.R-project.org/package=dgraphs) in its
examples. Extended worked examples are available in the
[grip vignettes](https://pgajer.github.io/grip/).

The sampled-saddle configurations were fitted with grip 0.2.0. Additional
reference diagnostics evaluated those saved coordinates with grip development
commit [`b72f61d`](https://github.com/pgajer/grip/tree/b72f61d9b5f20a822d3e87dacc1b45de025aabc7)
without refitting them. The triangulated interactive visualizations use ivue
commit [`872f9d4`](https://github.com/pgajer/ivue/tree/872f9d45827c7617005e7938f429a56d58b3e8b7).
The `reproducibility/README.md` file gives the full dependency record,
checksums, and commands for repeating the computational results.

## Supplement S1: Search-sensitive complexity of weighted GRIP

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

## Supplements S2 and S3

S2 provides the higher-dimensional Möbius-strip comparison, with canonical
source `S2-mobius-comparison.tex` and separate bibliography and citation evidence.

S3 documents the five independent, 1,000-observation sampled-saddle experiments:
graph calibration, numerical surface-reference checks, primary and additional
iteration results, identity-embedding controls, and end-to-end errors. It also
retains the executable variable-density circle example, full fixed-grid
saddle comparison, and edge-KK mesh illustration moved from the main text.
The mesh example uses `../reproducibility/scripts/edge-kk-workflow.R`;
the manuscript instead shows the sampled-saddle workflow. S3's canonical source is
`S3-controlled-examples.Rmd`. S3 also contains the reference table of method
families and principal functions, including their standard or experimental
status. `render-S3.R` renders it using R, rmarkdown,
bookdown, grip, dgraphs, igraph, Pandoc, and LaTeX. Python is needed to repeat
the surface-distance experiment, not to render the supplement.

The regular `make` target renders S3 from the supplied compact RDS and runs its
citation gate. From the repository root, its narrow build is:

```sh
Rscript papers/grip-software-paper/supplement/render-S3.R \
  papers/grip-software-paper/build/supplement
```

All supplementary PDFs and generated figures remain under `build/`. The
submission archive includes the source and input files needed to rebuild them.

## Supplement S4

`S4-interactive-saddle.html` provides self-contained interactive views of the
sampled-saddle configurations. Figure S4.1 carries a fixed triangulation from
the original parameter coordinates to the fitted layouts, Figure S4.2 overlays
the generating saddle surface, and Figure S4.3 compares metric-MDS with its
edge-KK refinement. The source is `S4-interactive-saddle.Rmd`; the reusable
widget builder is `../reproducibility/scripts/saddle-widgets.R`.
