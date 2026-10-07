# Sparse metric MDS: implementation and validation

[Methods note PDF](methods.pdf) · [LaTeX source](methods.tex) · [Citation verification](citation_verification.html) · [Build manifest](build_manifest.json)

Implemented October 7, 2026 in the local `grip` and `dgraphs` source packages.
The methods note has no abstract. It covers uniform landmarks, Euclidean regional
sparse stress and graph-geodesic regional sparse stress, including internally
selected pivots and regional weights. Component-MST repair streams distances
without retaining a full distance matrix. Its worst-case distance work remains
quadratic. Dense calculations occur only in the deliberately small references.

The final validation completed 288/288 fits on four 120-point quadratic forms:
latent dimensions 2 and 3, index 0 and 1, outputs 2 and 3, three seeds,
16 pivots, q=5/10/20 and 100 SGD epochs with the published weight-based schedule.
The report shows all cases and distinguishes Euclidean from geodesic targets.
It does not establish convergence, large-data performance or V4 parameter adequacy.

## Evidence

- [Final results and saved preparations](../../../output/sparse-mds-methods/validation_2026-10-07_final/validation.rds), [288 scores](../../../output/sparse-mds-methods/validation_2026-10-07_final/fits.csv), [source fingerprints](../../../output/sparse-mds-methods/validation_2026-10-07_final/source_manifest.csv), [session information](../../../output/sparse-mds-methods/validation_2026-10-07_final/session-info.txt).
- [Figures 1–4: primary comparisons](../../../output/sparse-mds-methods/error-comparison.pdf) and [Figures 5–8: paired neighborhood sensitivity](../../../output/sparse-mds-methods/neighborhood-sensitivity.pdf). Seed ranges are descriptive, not confidence intervals.
- [Full grip tests: 6,665 passed](../../../output/sparse-mds-methods/grip_full_tests_final.log); [dgraphs focused tests: 112 passed](../../../output/sparse-mds-methods/dgraphs_tests_final_corrected.log).
- [Final execution log](../../../output/sparse-mds-methods/validation_dispatch_final.log): 41.12 s serial process wall time, 443,318,272 bytes maximum resident memory, including small dense references and saved outputs.
- [Package check blocker](../../../output/sparse-mds-methods/check_fast.log): the existing `inst/extdata/hmp_u01_gc_coarse/weighted_layouts.rds` contains a personal path; repository hygiene stopped `make check-fast` before package build. That unrelated object was preserved. No CRAN-readiness claim is made.

The first benchmark dispatch failed before fitting because exports were stale.
An initial successful benchmark remains in `validation_2026-10-07_exported`.
A projected-coincidence unit test later exposed an overly strict search-coordinate
validator; after correction the final benchmark was repeated with source hashes.
Test-harness failures and export-inventory updates are preserved in the logs.
The final source/library hashes were checked both before and after fitting and
are checked again by the note builder. The production data were not analyzed.

## Reproduce

Use local source packages, with optimized native libraries compiled by
`pkgbuild::compile_dll(debug=FALSE, force=TRUE)` in each package, then
`pkgload::load_all(..., recompile=FALSE)`. Generate documentation with
`make document` before loading the new exported interface. Installed older
versions do not contain these extensions.

From `~/current_projects/vaginal_microbiome`, the benchmark dispatch was:

```sh
python3 scripts/storage_retention_v4.py run --local-peak-gib 2 --external-peak-gib 0 -- \
  /usr/bin/time -l Rscript ../grip/inst/scripts/sparse-mds-methods-validation.R \
  ../grip ../grip/output/sparse-mds-methods/validation_2026-10-07_final
```

For a new scientific rerun, choose a new output directory; the runner refuses
to overwrite an existing completed score table. Data seeds are
`20261020 + 10*d + index`, preparation seeds `20261030 + 1:3`, fitting seeds
`1:3`; initialization is matched across methods within each case, dimension and
seed. All settings and input coordinates are saved.

From the `grip` root, regenerate figures and render the note without rerunning fits:

```sh
Rscript --vanilla tools/pkg/summarize-sparse-methods.R \
  output/sparse-mds-methods/validation_2026-10-07_final output/sparse-mds-methods
python3 vignettes/articles/sparse-mds-methods/build_report.py --pdf
```

The builder regenerates result tables, refreshes the Eastern build timestamp,
checks saved implementation fingerprints, compiles twice, runs the citation gate
and writes SHA-256 fingerprints. Without `--pdf`, it refreshes the source for the
built-in editor. The citation gate is also available through `make citation-check`
in this note directory. Use the built-in `.tex` editor for source and live preview.
The linked figure supplements are separate PDFs so the main source remains
standalone for the native compiler.

The later schedule-documentation clarification changed comments in
`R/metric_mds.R` without changing its parsed R expressions. The exact benchmark-era
source is preserved in `output/sparse-mds-methods/documentation_revision_2026-10-07/`.
The builder accepts this documented difference only after checking the archived
file against the original run hash and independently verifying identical parsed
expressions; the build manifest records both hashes. Other implementation
changes still fail the source check.
