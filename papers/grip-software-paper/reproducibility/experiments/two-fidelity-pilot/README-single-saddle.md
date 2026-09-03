# Figure 7: one-sample interactive experiment

`single-saddle-ivue.R` repeats the Figure 7 experiment for **one** cloud, with a
configurable sample size, and creates interactive `ivue` views. It does not
edit the manuscript or overwrite the five-cloud pilot inputs. Sourcing it only
defines functions; it does not launch a computation or a browser.

## Where Figure 7 comes from

- `../../scripts/plot-two-fidelity.R`: `plot_pilot_calibration_layouts()` draws
  the publication figure. The Rmd supplies its caption.
- `export-paper-data.R`: creates `../../precomputed/two-fidelity-saddle.rds`,
  the compact input read by the manuscript. It selects cloud 5 as the median
  cloud by minimum surface-to-graph error, not by appearance.
- `sample-clouds.R`, `surface-reference.py`, `graph-sweep.R`, `calibrate.R`,
  and `fit-layouts.R`: generate the original samples, numerical surface
  distances, graphs, calibration results, and fitted coordinates.
- `path-lengths.cpp`: independent sums along the retained graph paths.

The new R script follows these calculations, replacing the original fixed
sample count and size with parameters. `single-saddle-reference.py` calls the
existing Python reference functions and adapts their validation to one cloud
of arbitrary size. The original experiment scripts are unchanged.

## Dependencies

Use `grip >= 0.2.0`, `dgraphs >= 0.2.0`, `igraph`, `Rcpp`, `htmlwidgets`,
`htmltools`, `rgl`, and `pkgload`. By default the script loads `ivue` directly
from `~/current_projects/ivue`; set `ivue.source = NULL` to use an installed
version instead. No package source is edited and no native rgl window is opened.

If the correct grip/dgraphs versions are in a separate library, set
`GRIP_RJOURNAL_PACKAGE_LIBRARY` before launching R, or pass
`package.library = "/path/to/library"` to `single.saddle.config()`. The script
stops early when it finds an older version; it does not silently change methods.

Python requires the packages in the adjacent `requirements.txt`. The script
uses `GRIP_SADDLE_PYTHON`, then the existing pilot environment at
`../../../build/two-fidelity-pilot/venv/bin/python`, then `python3` on PATH.
Override with `python = "/path/to/venv/bin/python"`. No automatic installation
occurs. Rcpp route checks require a working C++ toolchain.

## Quick trial in R

```r
source("/Users/pgajer/current_projects/grip/papers/grip-software-paper/reproducibility/experiments/two-fidelity-pilot/single-saddle-ivue.R")

cfg <- single.saddle.config(
  n = 200,
  cloud = 5,
  ks = 3:40,
  reference.sources = 16,
  reference.grids = c(21, 41),
  fine.grid = 61,
  fine.sources = 8,
  max.iter = 30,
  audit.iter = 50,
  plot.ks = c(3, 10, 20, 40),
  out = "~/current_projects/grip/papers/grip-software-paper/build/saddle-trial-n200"
)
result <- run.single.saddle(cfg)
result$calibration$selected.k
result$calibration$curve
result$reference$validation
result$fits[[1]]$scores
result$fits[[1]]$edge.metadata
result$fits[[1]]$audit$scores
```

These reduced grid sizes and iteration budgets are for testing, not the
published numerical protocol. The configuration can be modified directly
before the run, but calling `single.saddle.config()` again provides input
validation. Choose a **new output directory** whenever scientific settings,
software versions, or source code change. Existing numerical results are
reused only when their saved protocol matches. View-only settings may change
without discarding the experiment.

## Original Figure 7 cloud and primary settings

```r
cfg <- single.saddle.config(n = 1000, cloud = 5)
result <- run.single.saddle(cfg)
```

Defaults: `ks = 3:80`, 128 reference sources, grids 41/81, a 161-grid check on
16 sources, three-dimensional metric MDS, and the five-stage density mixing
schedule `c(0, .25, .5, .75, 1)` with at most 200 edge-KK iterations per stage.
A separate 1,000-iteration continuation is saved as an audit, not substituted
for the primary fit. The default sample coordinates exactly match the
published cloud 5. `cloud` selects one reproducible seed family; **no other
cloud is generated**. For another sample use, for example, `cloud = 6`.

The original n=1,000 pilot used about 150 seconds per selected-graph fitting
run on its recorded environment, in addition to reference and graph
construction. Timings depend on hardware and software. Dense scaling and
all-pair route preparation become expensive as n grows. This exploratory
script does not provide a large-n approximation or a memory guarantee.

Command-line equivalent, from any working directory:

```sh
Rscript /Users/pgajer/current_projects/grip/papers/grip-software-paper/reproducibility/experiments/two-fidelity-pilot/single-saddle-ivue.R 1000 5 /path/to/new-output
```

## Run each part independently

```r
cfg <- setup.single.saddle(single.saddle.config(n = 1000, cloud = 5))
cloud <- sample.single.saddle(cfg)
reference <- reference.single.saddle(cfg)
calibration <- graphs.single.saddle(cfg, cloud, reference)
fits <- fit.single.saddle(cfg, cloud, reference, calibration)
views <- view.single.saddle(cfg, cloud, calibration, fits)
browseURL(views$index)

# Returned widgets can also be printed directly in an RStudio viewer.
views$widgets[["original-saddle"]]
names(views$widgets)
```

The built-in checks cover connectedness and MST counts, Euclidean input edge
lengths, reference-pair bookkeeping, mesh/BVP/plane controls, agreement of
prepared and strict graph distances, independent C++ and R fixed-route sums,
agreement with `score.gmds()`, and zero edge/path error in the original
coordinates. Scores use all `choose(n, 2)` unordered pairs and all graph edges.

`reference$validation` reports numerical resolution sensitivity; it does not
certify exact distances on the smooth surface. `calibration$boundary` flags
a minimum at the edge of the tested k range. The script still fits that graph
so an exploratory run can be inspected, but warns that the range should be
expanded before interpreting the minimum. With `ks` restricted to a sparse
set, selection is only among the supplied values.

## Compare graphs and change the views

- `ks`: graphs to construct and compare against surface distances.
- `plot.ks`: graphs to show in their original saddle coordinates. NULL chooses
  a small set including the selected k. Set `plot.ks = cfg$ks` for every graph.
- `fit.ks`: graphs to fit with MDS and MDS→edge-KK. NULL fits the calibrated
  graph only. To compare embedding methods at several k values, supply, for
  example, `fit.ks = c(10, 20, 40)` and include them all in `ks`.
- `edge.alpha`: background-edge opacity (default 0.03). Every fitted
  configuration also gets a point-only view. The highlighted route and its
  endpoint chord remain visible in these views.
- `align.display`: similarity-align fitted coordinates to the generating
  sample for display (default TRUE); all scores use unaligned coordinates.

```r
# Redraw from a saved run, with no reference, graph, or layout recomputation.
result <- readRDS(file.path(cfg$out, "experiment.rds"))
result$config$edge.alpha <- 0.01
result$config$plot.ks <- c(3, 10, 20, 40, 73)
result$config <- setup.single.saddle(result$config)
views <- with(result, view.single.saddle(config, cloud, calibration, fits))
browseURL(views$index)
```

`ivue::plot3D.cont()` renders points and `ivue::plot3D.graph()` renders graphs.
`color.scale.cont()` preserves original x-coordinate colors across all views.
`layer3D.path()` / `layer3D.edges()` show a route, its chord, and MST bridges;
the `layer3D.callback()` escape hatch explicitly sets rgl line opacity.
Dragging rotates the scene and scrolling zooms. Axes have equal data-unit
scales, but each widget auto-fits its own bounds: apparent sizes between
separate windows are not an additional fidelity measure.

The original-coordinate graph views retain the same endpoints but find a route
on each respective k graph. Within a fitted-graph comparison, the route is
held fixed across all three configurations, never recomputed after embedding.

## Saved outputs

- `views/index.html`: links, settings, scores, and interaction instructions.
- `views/*.html`: self-contained interactive widgets; no server is required.
- `views/calibration.pdf`: one-cloud error-versus-k curve.
- `calibration.csv`, `reference-checks.csv`, `layout-scores.csv`: numeric results.
- `graphs/k-*.rds`: graph edges, Euclidean weights, and MST information.
- `fit-k-*.rds`: candidate coordinates, scores, highlighted route, iteration
  traces/metadata, timings, and the separate extra-iteration audit.
- `experiment.rds`: combined numerical result, excluding the duplicated widget
  scenes. Load it and rerun only the viewing stage to change presentation.
- `config.rds`, `session-info.txt`, `reference-environment.json`: settings,
  source checksums, and software versions.

For a served index, serve the **experiment output directory**, not just its
`views` subdirectory, so links to neighboring CSV files work. For example:

```sh
python3 -m http.server 8769 --bind 127.0.0.1 --directory /path/to/new-output
```

Then open `http://127.0.0.1:8769/views/index.html`.

## Automated checks

```sh
Rscript /Users/pgajer/current_projects/grip/papers/grip-software-paper/reproducibility/experiments/two-fidelity-pilot/check-single-saddle.R
Rscript /Users/pgajer/current_projects/grip/papers/grip-software-paper/reproducibility/experiments/two-fidelity-pilot/check-single-saddle.R --integration
```

The fast test verifies parameter validation, RNG-state restoration, an exact
match to the published cloud-5 coordinates, and pair accounting at several n.
The integration test runs every stage on n=80, with both sparse and denser
graphs, two MDS/edge-KK fits, an extra-iteration check, and widget generation.
It checks the captured line opacity as well. Its temporary files are only
test outputs, not a replacement for the n=1,000 experiment.
