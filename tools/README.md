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
results from `inst/extdata/metric-mds-comparison/` and runs only its small example.

To run a new, isolated timing experiment, install the candidate package into a
separate library, then call:

```
Rscript tools/pkg/build-metric-mds-comparison.R <fresh-output-dir> <candidate-library> <install-log>
Rscript tools/pkg/render-metric-mds-comparison.R <fresh-output-dir>/benchmark.rds <render-dir>
node tools/pkg/capture-metric-mds-comparison.cjs <render-dir> <image-dir>
```

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
images. Each graph runs in a separate process with a 180-second limit, and
failures are saved while later graphs continue. Edge-KK failures retain the
completed MDS layout. Settings and all MDS start summaries are retained.
