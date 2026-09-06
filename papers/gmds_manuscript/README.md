# Graph-geodesic embedding methods paper

The active source is `geodesic_mds.tex`. It is being refocused as
**What Graph-Geodesic Embeddings Preserve: Path Fidelity, Scale, and Geometric
Recovery**, using the completed grip initializer and radius studies.

The focused version is an author-review draft. Current citation support,
algebraic diagnostics and 58 numerical claims are checked against frozen data,
and a focused primary-source related-work review is recorded. An independent
mathematical and numerical audit has informed the revision.
Final author approval remains a separate requirement.
Do not describe it as submission-ready.

## Canonical sources and provenance

- One active manuscript: `geodesic_mds.tex` and `geodesic_mds.bib`.
- The pre-refocus manuscript is preserved in Git at commit `740df96`.
- Historical snapshots under `archive/` remain historical, not active versions.
- `evidence/source-manifest.json` pins the grip source commit and file hashes.
- `evidence/claims.json` records the retained claims and their assumptions.
- `evidence/study-initializers/` and `evidence/study-radius/` keep study identities
  separate. `evidence/surface-theory/` preserves the analytical source report.
- Generated focused figures and tables are in `figures/focused/` and
  `tables/focused/`; older figures are retained for history but not included
  automatically in the focused paper or its eventual submission bundle.
- `citation_verification.html` maps every current citation to source evidence.
- `evidence/numerical-claims.json` maps generated numerical macros to the
  selection rules and source tables, distinguishing recomputation from recorded
  validation summaries. `evidence/spatial-display-checks.json` records the
  orthogonal transformations and uniform scales used in the spatial figure.

## Build and check

From this directory:

```sh
make pdf
make verify
```

Figure generation requires Python with NumPy, pandas and Matplotlib. The PDF
build requires zsh and a LaTeX installation with latexmk. The build uses latexmk on PATH, with a macOS TeX fallback; set `LATEXMK`
to select another executable. Builds work from a
clean worktree and do not require a neighboring grip checkout or private files.
Build timestamps use DST-aware America/New_York wall time.

The PDF is written to `build/geodesic_mds.pdf`. Building verifies LaTeX syntax;
visual inspection of the resulting PDF remains necessary. The numerical scripts
check exported byte identities, declared counts, reported paired comparisons,
and explicit algebraic identities/counterexamples. Numerical diagnostics do not
replace mathematical proofs or an independent scientific review.

See `evidence/README.md` for explicit refresh commands. Refreshing source exports
requires a specified grip checkout and rejects selected files with uncommitted
changes. Ordinary rendering does not refit experiments.

## Scope

The manuscript analyzes endpoint, fixed-path, graph-reference and coordinate
fidelity, the scale behavior of edge stress, and the expanding paraboloid/saddle
limits. It does not introduce a new SMACOF algorithm or promise general manifold
recovery. The earlier unsupported path-augmented update is excluded from the
focused draft; its counterexample is preserved by the mathematics checker.

The broader project in `../data_to_geodesic_embedding_paper/` and the grip
software article remain separate. The methods manuscript must keep its own
contribution clear while documenting shared experimental material.

## Author-review packaging

`make review-bundle` writes the PDF, a minimal typesetting source ZIP, a
frozen-evidence reproduction ZIP, and an artifact hash manifest under
`build/review-bundle/`. The packager uses an explicit active-source list and
excludes historical manuscripts, unrelated figures, and private review notes.
Each archive includes its own file-hash manifest. The reproduction archive
regenerates figures and verifies frozen numerical evidence; it does not refit
the original experiments. The source provenance identifies those fitting
scripts separately. Independent author review is still required before any
submission, even when these computational checks pass.

The September 6 audit revision adds component-wise coordinate diagnostics, local
rigidity checks, contraction prevalence, measured neighborhood sensitivity and
paired stiffness comparisons. Distribute one complete `review-bundle` directory;
do not mix a PDF from one dated revision with archives from another.
