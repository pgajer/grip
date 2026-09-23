# Package tools

`tools/pkg/` contains package documentation and release utilities, including
self-contained inputs for the public README previews. Repository checks live
at this directory's top level.

Manuscript experiments, figure builders, reports, and submission utilities
belong inside their owning private manuscript workspace outside this repository.
The public package must build and pass checks without those workspaces.

`make repo-hygiene` checks public source placement, development-file hygiene,
Rd documentation for exported functions, and private-data/reference exclusion.
`make test-release-content` exercises the release-content check in isolated
fixtures, including absent documentation and private-path failures. Historical
manuscript catalogs are not inputs to either target.

## Metric-MDS comparison vignette

`inst/scripts/metric-mds-comparison.R` contains the public-API synthetic-data,
fitting, scoring and plotting recipes. The installed vignette renders reviewed
results from `inst/extdata/metric-mds-comparison/`; its full protocol and all cases
are retained in the installed `metric-mds-results` appendix.

To run a new, isolated timing experiment, install the candidate package into a
separate library, then call:

```
Rscript tools/pkg/build-metric-mds-comparison.R <fresh-output-dir> <candidate-library> <install-log>
Rscript tools/pkg/render-metric-mds-comparison.R <fresh-output-dir>/benchmark.rds <render-dir>
node tools/pkg/capture-metric-mds-comparison.cjs <render-dir> <image-dir>
```

For the separate 60-fit weighting comparison on the 64-vertex saddle,
paraboloid, and helix graphs, add `weighting` as the runner's fourth argument.
Both backends fit both objectives from five paired starts; every result is
scored under both objectives. The same renderer detects that bundle and writes
four-panel comparisons plus 3D overlays. Reviewed weighting artifacts belong
in `inst/extdata/metric-mds-weighting/`; the original uniform experiment stays
in `inst/extdata/metric-mds-comparison/`.

The explicit benchmark needs `callr`, `smacof` and the candidate package; the
renderer needs `ivue`, `rgl`, `htmlwidgets` and Pandoc; image capture needs
Playwright and Chromium. The script preserves each attempt and records failures.
Review results before replacing the public bundle or images. Raw worker logs,
intermediate attempts and browser audit records belong outside the repository.
The reviewed bundle, images and recipes are public vignette artifacts.

## Graph examples gallery

The six SuiteSparse graphs are installed as `zheng.graphs`. Rebuild them
without downloading anything using `Rscript data-raw/zheng_graphs.R`; the
archived matrices and attribution are in `inst/extdata/zheng-graphs/`.

The gallery includes six generated graphs, karate club, and these six datasets.
To recreate all 3D layouts, install this checkout into a separate library and run:

```
Rscript tools/pkg/build-graph-examples.R <fresh-output-dir> <candidate-library>
Rscript tools/pkg/render-graph-examples.R <fresh-output-dir>/layouts.rds <render-dir>
node tools/pkg/capture-metric-mds-comparison.cjs <render-dir> <image-dir>
```

The fitter needs `callr`; rendering needs `ivue`, `rgl`, `htmltools`,
`htmlwidgets`, and Pandoc. The shared capture script needs Playwright and Chromium;
`PLAYWRIGHT_MODULE` and `CHROMIUM_EXECUTABLE` can locate them. Review the layouts
and scores before replacing `inst/extdata/graph-examples/layouts.rds` and its
images. Rendering uses the saved layouts without refitting. Each interactive view
has a selector for metric-MDS, metric-MDS + edge-KK, or both; captures show both.
The capture script checks visibility, fixed camera and coordinates, rotation,
and the absence of legend tables. Use a fresh render directory to avoid
capturing obsolete figures. Each graph runs in a separate process with a
180-second limit, and failures are saved while later graphs continue. Edge-KK failures retain the
completed MDS layout. Settings and all MDS start summaries are retained.


## Sparse SGD comparison vignette

After a clean optimized installation as above, run:

```sh
Rscript tools/pkg/build-sparse-mds-comparison.R <fresh-dir> <candidate-library> <install-log>
Rscript tools/pkg/render-sparse-mds-comparison.R <fresh-dir>/benchmark.rds <render-dir>
node tools/pkg/capture-metric-mds-comparison.cjs <render-dir> <image-dir>
```

The fixed 60-fit experiment uses three seeds, 100 SGD epochs, several pivot
counts, and weighted GRIP defaults. A worker gets 90 seconds; the run gets 20
minutes. It measures worker peak resident memory using `/usr/bin/time` (macOS
or Linux), preserving failed attempts and logs. Quality evaluation uses igraph
outside the timed workers. Source-sampled distances for the large graph are
independent of fitting pivots. Reviewed coordinates, scores and figure assets
belong in `inst/extdata/sparse-mds-comparison/`. Full worker results and logs
remain private; ordinary vignette builds only read the reviewed bundle.

## Documentation collection

All teaching material is maintained as installed `vignettes/*.Rmd` sources.
Pkgdown publishes those same sources. Ten primary guides are followed by four
appendices: further graph examples, the full MDS experiment, the full sparse
experiment, and the historical grip 0.2.0 engine comparison. `_pkgdown.yml`
controls both the site navigation and local preview index.

Run `make vignette-previews` to render all pages, local help, an index, a build
manifest, and redirects for former article locations. This also checks the
complete API catalog, selector behavior, archived checksums, and local links.
The output is `output/vignette-previews/index.html`. For another destination:

```
Rscript tools/pkg/render-vignette-previews.R /absolute/output/directory
```

Routine rendering uses saved experiment results and surface registrations;
it runs only the modest teaching examples. `make documentation-assets` is the
explicit, more expensive step for rebuilding surface registrations after
changed inputs or alignment code. It does not rerun optimizers or benchmarks.
The saved transforms have input/algorithm fingerprints checked at render time.
Static fallbacks use the current display coordinates rather than old screenshots.

The immutable snapshot in `vignettes/archives/` records the sources before this
restructuring. Its underscored Rmd names prevent pkgdown discovery and the whole
directory is excluded from R builds. Generated archive HTML remains local.
Do not regenerate or edit archived sources when updating current documentation.

## README and showcase assets

`make readme-assets` generates the recursive triangle/carpet construction
animations and their provenance. Their current teaching counterpart is the
recorded-stage slider in the tracing vignette. The old saddle showcase is a
historical display fixture, with its known and missing metadata retained in
`tools/pkg/fixtures/saddle-preview.json`; it is not evidence about current SGD.
The former showcase URL redirects to the tracing guide. Rendering commands for
the historical saddle remain the explicit `readme-saddle-animation` and
`saddle-overlay-animation` Makefile targets.
