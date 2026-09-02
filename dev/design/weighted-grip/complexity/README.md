# Weighted GRIP search-sensitive complexity

This directory retains the original standalone technical note. The maintained
publication-facing version is now
`papers/grip-software-paper/supplement/S1-weighted-grip-complexity.tex`
relative to the repository root; build it with `make paper-supplement`.
Future manuscript-related revisions should be made there, not in this snapshot.

This user-requested technical note derives workload-sensitive estimates for the
existing weighted 2D/3D GRIP core. It is a companion derivation, not a change to
the R Journal manuscript, a package change, or a new benchmark.

Canonical source: `weighted-grip-complexity.tex`, with its bibliography and
`citation_verification.html`. Source behavior is pinned to commit `cebcc62`.
The ball-overlap proof, cache accounting, and conditional coverage estimates
are derived in the note, not attributed to the original GRIP publications.

From this directory, run `make`. The PDF and intermediates are generated in
`output/pdf/weighted-grip-complexity/` relative to the repository root. The build
refreshes the Eastern-time metadata and runs the citation verification gate.
For an independent directory, use `make pdf BUILD_DIR=build`; the PDF needs
LaTeX with latexmk, BibTeX, and the packages listed in the source.

The source bundle includes the citation checker and supports
`make BUILD_DIR=build CHECKER=check-citation-verification.py` after extraction.
No private files or benchmark outputs are required to compile the document.
