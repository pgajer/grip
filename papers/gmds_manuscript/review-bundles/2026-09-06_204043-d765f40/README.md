# Graph-geodesic methods paper: author-review package

Start with `geodesic_mds.pdf`. The focused paper combines the completed MDS
initializer and radius studies, with their experimental designs kept separate.

- `typesetting-source.zip`: minimal active TeX/BibTeX/bibliography, six figure
  PDFs, five generated table/macro files, and a Makefile. Extract and run
  `make pdf` with latexmk on PATH.
- `evidence-reproduction.zip`: active source, figure builders, validation
  scripts, citation-support HTML, frozen numerical exports and provenance.
  Requires Python with NumPy, pandas and Matplotlib, zsh, and latexmk. Extract,
  run `make verify`, then `make pdf`. To regenerate figures explicitly, run
  `make figures` first. Set PYTHON and LATEXMK if needed.
- `artifact-manifest.json`: paper commit, branch, input counts and SHA-256 hashes.
  Each source archive also has `bundle-files.json` with per-input hashes.

The reproduction archive supplies explicit figure regeneration and checks selected
claims from saved evidence. Ordinary PDF builds use the frozen figure assets. It does not rerun the original MDS/geodesic/edge-KK fits.
Their pinned source paths and commits are in the evidence provenance records.
The original study documents are preserved as source records; their relative
links refer to the source study trees and are not a complete mirror here.

This is an author-review draft. Internal checks are complete; independent
scientific review informed this revision; approval by both authors, disclosure review, overlap review
with the software article, and submission/category/license choices remain.
The canonical source is grip/papers/gmds_manuscript. No submission has been performed.

Distribute this directory as one matched set. Its PDF, source archives and
manifest describe the same revision; older dated review folders are historical
artifacts and must not be substituted for any member of this set.
