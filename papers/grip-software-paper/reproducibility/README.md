# Reproducing the `grip` software-paper results

This directory makes the precomputed components of the paper inspectable and
reproduces all reported results, including reconstruction of the UMB-HMP-only graph
from the supplied minimal upstream tables. The manuscript evaluates the smaller
examples directly.

## Requirements

### Single-graph example for Figure 6E

Subsection 3.2 uses `scripts/panel-e-workflow.R` to reproduce the representative
metric-MDS + edge-KK configuration (cloud 5, n = 1,000, k = 73). From the
manuscript directory, or the root of the extracted submission, run:

```r
source("reproducibility/scripts/panel-e-workflow.R")
```

This prepares all 499,500 vertex pairs, fits metric-MDS, runs the original five
200-iteration edge-KK stages, scores the unaligned coordinates, and creates the
interactive `panel.e` widget. It neither repeats the graph-calibration sweep
nor overwrites the supplied fits or figure assets. Allow roughly two minutes
on a machine comparable to the original pilot host; runtime is machine-dependent.
The package versions used for the original fits and visualizations are recorded
below. The graph, seed, stage schedule, profiled scales, and zero edge-length
epsilon match `experiments/two-fidelity-pilot/fit-layouts.R`.

Normal article builds print this code but use the supplied scores and figure
assets. To redraw the saved configuration without fitting, run the input block,
set `fit <- list(coords = example$candidates[["MDS + edge-KK"]])`, and run the
view block. `scripts/panel-e-display.R` implements display-only similarity
alignment, the saved parameter-plane triangulation, original-x colors, and the
common bounds of panels C--E. It does not change graph lengths or scores.
The ivue view uses the same mesh, retained path, endpoint chord, and camera as
Figure S4.1C; the paper's static panel is a cropped capture of that scene.
Optional HTML export is `htmlwidgets::saveWidget(panel.e, "panel-e.html")`.

Fitting and scoring require grip; visualization additionally requires ivue and
rgl, with htmlwidgets for export. The view uses the saved Delaunay triangles,
so geometry is required to regenerate those triangles, not to load them.
CRAN availability of the required ivue functionality remains a manuscript
submission prerequisite unless an editorial exception is obtained; a pinned
development source is not represented here as a CRAN release.

### Additional reference diagnostics and interactive views

The original five-cloud fits remain grip 0.2.0 results. New coordinate and
surface scores evaluate those saved configurations with grip commit `b72f61d`.
The compact `precomputed/saddle-reference-diagnostics.rds` contains all five
clouds and coordinates, fixed triangles, scores for rigid and similarity
alignment, source-fit checksums, and software identifiers. Its companion CSV
contains all scalar scores. `precomputed/single-saddle-reference.csv` is the
separate n=500, k=3:20 illustration, not part of the five-cloud pilot.

From the repository root, regenerate these additions with:

```sh
Rscript papers/grip-software-paper/reproducibility/experiments/two-fidelity-pilot/export-reference-diagnostics.R
Rscript papers/grip-software-paper/reproducibility/experiments/two-fidelity-pilot/check-reference-diagnostics.R
IVUE_SOURCE=/path/to/ivue Rscript papers/grip-software-paper/supplement/render-S4.R
node papers/grip-software-paper/reproducibility/scripts/capture-saddle-widgets.cjs
Rscript papers/grip-software-paper/reproducibility/scripts/reference-manifest.R
make paper-all
make readme
```

The exporter reads existing `build/two-fidelity-pilot/fit-*.rds` files; use the
pilot protocol below to recreate them. It also reads the separately generated
500-point scores from `run-saddle-reference.R`. Ordinary document rendering
uses the supplied compact inputs and screenshots and does not rerun fitting,
reference scoring, or a browser. S4 rebuilding requires ivue commit `872f9d4`,
rgl, htmlwidgets, geometry, and rmarkdown; `IVUE_SOURCE` optionally selects a
source checkout. The screenshot script requires Playwright and Chromium;
`PLAYWRIGHT_MODULE` and `CHROMIUM_EXECUTABLE` optionally specify their locations.
It uses an identical crop for every view and verifies offline JavaScript startup.
The exporter optionally accepts the fit directory and the separate single-cloud
output directory as its first and second arguments. `reference-assets-md5.tsv`
records the distributed reference-score and visualization assets.

The additional APIs are defined in
[grip's pinned reference-scoring source](https://github.com/pgajer/grip/blob/b72f61d9b5f20a822d3e87dacc1b45de025aabc7/R/reference_scores.R).
The renderer uses the
[pinned ivue source](https://github.com/pgajer/ivue/tree/872f9d45827c7617005e7938f429a56d58b3e8b7).

`supplement/S4-interactive-saddle.html` is a deliberately distributed,
self-contained HTML supplement, generated from its Rmd source. The Pages
workflow copies it to `supplements/` without requiring ivue on the site builder.
The README image is regenerated from the same reviewed screenshots.

Reference diagnostics use 2,000 and 8,000 samples per direction, with two
reference-mesh subdivisions; scores retain Monte Carlo standard errors.
These measure integration variability, not reference triangulation bias or
between-cloud uncertainty. Fixed connectivity is a display/scoring mesh, not
the calibrated skNN graph. All original graph-score values and saved method
keys remain unchanged. See S3 for formulas and baseline interpretation.

### Original experiment requirements

- R and `grip` 0.2.0 or later;
- `dgraphs` 0.2.0 or later for the manuscript examples and to rebuild the
  UMB-HMP graph from the upstream count and metadata tables;
- `igraph` and `graphlayouts` for the cross-package benchmark;
- the packages listed in the manuscript-level `_Rpackages.txt`.

The article describes `grip` 0.2.0 and uses `dgraphs` 0.2.0 in its examples.
The corresponding source-release tags resolve to the following commits:

- [`grip` v0.2.0](https://github.com/pgajer/grip/tree/v0.2.0):
  [`520902ad3f1b2aeabd287a379f5c08729c7b2c5d`](https://github.com/pgajer/grip/tree/520902ad3f1b2aeabd287a379f5c08729c7b2c5d);
- [`dgraphs` v0.2.0](https://github.com/pgajer/dgraphs/tree/v0.2.0):
  [`8733d2a74dc09d57b453b88eff3119610c6440f3`](https://github.com/pgajer/dgraphs/tree/8733d2a74dc09d57b453b88eff3119610c6440f3).

The full commit identifiers pin the source trees independently of the tag
names. The environment used to generate the precomputed benchmark results is
recorded separately in `BENCHMARK_PROVENANCE.md` and the RDS metadata.

The `grip` commit identifier above reflects the September 3, 2026 repository
history cleanup, which removed internal working records outside the R package.
The cleanup did not change the package files in the release snapshot or the
published CRAN archive; see the
[source-history notice](https://github.com/pgajer/grip/blob/main/dev/release/source-history.md).
The SHA-256 of the published `grip_0.2.0.tar.gz` archive is
`5fbaecee890cd1f402e6776c661f5934325d37ac93a7d39120b54cf7acaf31a8`.

## Included artifacts

- `precomputed/two-fidelity-saddle.rds`: compact manuscript and S3 inputs for
  the sampled-saddle experiment, including representative coordinates, graphs,
  all calibration curves, and primary/additional-budget scores and checks;
- `precomputed/vs_alternatives/benchmark_results.rds`: cross-package layouts,
  scores, timings, and generation-session information;
- `BENCHMARK_PROVENANCE.md`: timing boundary, repeat policy, hardware,
  software versions, and the benchmark medians reported in the paper;
- `hmp_gc/`: the exact UMB-HMP-only Illumina 16S edge list, vertex metadata,
  feature-screening manifest, graph summary, and provenance record;
- `hmp_gc/upstream/`: the minimal UMB-HMP-only feature-count and technical
  metadata tables needed to reconstruct that graph, together with an archive
  manifest, checksums, provenance and licensing information;
- `SHA256SUMS`: checksums for every supplied binary or tabular input and
  precomputed result.

The manuscript uses the full 4,391-vertex giant component of the UMB-HMP-only
graph. UMB-HMP denotes the Ravel-led longitudinal vaginal microbiome HMP
demonstration project, not the NIH HMP healthy-reference cohort.

## Paper-to-artifact map

| Manuscript component | Supplied input | Regeneration or validation path |
|---|---|---|
| Sampled-saddle calibration, layouts, and paired-score figures | `precomputed/two-fidelity-saddle.rds` | `scripts/plot-two-fidelity.R` renders; `experiments/two-fidelity-pilot/` regenerates and validates the experiment |
| Fixed weighted graph comparison (Supplement S3) | `precomputed/vs_alternatives/benchmark_results.rds`, components `weighted_saddle` and `weighted_saddle_resolutions` | `scripts/compare-saddle-resolutions.R` and its numerical/plotting helpers, also used by the full benchmark |
| Executable variable-density circle example (Supplement S3) | Generated directly by `supplement/S3-controlled-examples.Rmd` | Render S3; no precomputed input is used for this example |
| UMB-HMP-only graph construction summary | `hmp_gc/` graph files and `hmp_gc/upstream/` input tables | `scripts/build-hmp-only-graph.R` |
| Repeated HMP runtime and shared-scoring benchmark | `precomputed/vs_alternatives/benchmark_results.rds`, component `hmp` | `scripts/precompute-vs-alternatives.R` |
| Two-panel UMB-HMP layout figure | `precomputed/vs_alternatives/benchmark_results.rds`, component `hmp$layouts`, and `hmp.gc$vertex_data$cst` from the package data | `scripts/precompute-vs-alternatives.R` supplies coordinates; the manuscript plots them with the supplied CST labels |
| Benchmark hardware, software, timing boundary, and repeat policy | `BENCHMARK_PROVENANCE.md` and the RDS `benchmark_metadata` component | Recorded by `scripts/precompute-vs-alternatives.R` |
| Small-graph figures and tables | Evaluated directly from `grip-software-paper.Rmd` | Render the manuscript; no precomputed supplement input is used |
| Two-edge angle and coincident-endpoint diagram | Analytic three-vertex, unit-edge path | The `path-angle-freedom` chunk checks all three fixed-path distances at 1,001 deterministic angles (including zero), then draws the acute-angle and coincident-endpoint configurations; no random sampling or fitted layout is involved |
| Higher-dimensional Möbius-strip comparison | Unit-edge grid with a reversed seam, generated from dimensions in the script | The `mobius-data` chunk regenerates the 300-vertex pair using `scripts/mobius-dimension-comparison.R`; `mobius-dimension` draws it and computes its caption's retained-variance percentage from the same coordinates. Supplement S2 uses the same functions for 150, 300, and 1,500 vertices |

The table identifies the manuscript result that consumes each supplied
artifact. Figure and table numbering may change during editing, so the mapping
uses stable section descriptions and RDS component names.

## Sampled-saddle experiment (main text and Supplement S3)

The manuscript renders both sampled-saddle figures from the compact RDS using
`scripts/plot-two-fidelity.R`. S3 reads the same results and evaluates the
smaller circle example directly. Neither render repeats the five-cloud fits
or requires Python or the multi-gigabyte pilot working directory.

`experiments/two-fidelity-pilot/README.md` documents full regeneration, including
surface-area sampling, exact symmetric-neighbor search, MST repair, the
adaptive k range, numerical surface-reference checks, fixed-path scoring,
refinement budgets, and environment records. After that sequence, export the
publication input from the manuscript directory with:

```sh
Rscript reproducibility/experiments/two-fidelity-pilot/export-paper-data.R \
  build/two-fidelity-pilot reproducibility/precomputed/two-fidelity-saddle.rds
```

The export is checked against the validated scores. The representative sample
is selected by median minimum surface-to-graph loss, not appearance. Full
per-sample results distinguish numerical reference error, source-pair sampling,
graph selection, and finite-budget optimization. Five samples illustrate
behavior; they do not establish precise population estimates or a general
method ranking.

## Higher-dimensional layout illustration (Supplement S2)

From the manuscript directory, regenerate the full six-panel comparison with:

```sh
Rscript reproducibility/scripts/mobius-dimension-comparison.R build/supplement
make -C supplement
```

The script writes a vector PDF, a PNG preview, and an RDS containing the graphs,
3D and 4D coordinates, PCA projections, display coordinates, explicit settings,
and generation-session information. The submission archive includes the figure
PDF and RDS beside `supplement/S2-mobius-comparison.pdf`. The manuscript computes
its 300-vertex pair directly using the same script; it does not read this RDS.

This is a modern illustration inspired by Gajer et al. (2004), not a replication
of the historical executable, graphs, or projection procedure. Both dimensions
use the current `weighted.grip.nd()` backend with unit edge lengths and seed 1;
the other shared settings and graph construction are explicit in the script.
PCA retains three components without scaling the four input axes. Display
normalization and orthogonal alignment occur only afterward; the analytic
Möbius strip is never supplied as layout initialization. The saved
`projection.edge.ratio` values compare projected edge lengths with their
original 4D lengths before display normalization. High retained coordinate
variance and smoother-looking drawings do not establish preservation of
individual edges or improved graph-distance fidelity. No seed search or
multi-seed performance comparison was performed.

## Software environment and validation

`BENCHMARK_PROVENANCE.md` records the released R version, operating system,
hardware, BLAS, package versions, timing boundary, and random-number policy used
for the supplied benchmark artifact. The RDS file also retains its generation
session metadata. The complete package list needed to render the article is in
the manuscript-level `_Rpackages.txt` file.

From this directory, verify the distributed inputs and artifacts with:

```sh
shasum -a 256 -c SHA256SUMS
```

The graph builder validates dimensions, identifiers, edge endpoints, edge
weights, and the feature manifest against the supplied HMP graph. The benchmark
script records the software session with its regenerated artifact.

## UMB-HMP-only graph

`scripts/build-hmp-only-graph.R` documents and executes the cohort filter,
feature screening, and graph construction from the supplied minimal UMB-HMP-only
tables. It retains only rows explicitly labeled as HMP and Illumina, then
constructs the symmetric 3-nearest-neighbor graph described in
`hmp_gc/PROVENANCE.md`:

```sh
Rscript scripts/build-hmp-only-graph.R
```

The minimal inputs exclude U01 records, clinical fields, participant
identifiers, and host sequence. They were extracted from the historical
combined working tables by `scripts/prepare-hmp-upstream-deposit.R`; the
combined tables are neither required nor distributed. A clean rebuild from the
minimal inputs reproduces the exact vertex identifiers, edge endpoints, edge
weights, and feature manifest used in the paper. Paths can still be overridden
with `GRIP_HMP_METADATA_TSV`, `GRIP_HMP_FEATURE_MATRIX_TSV`, and
`GRIP_HMP_OUTPUT_DIR`.

## Cross-package benchmark

From this directory, run:

```sh
GRIP_VS_ALTERNATIVES_OUTPUT=generated/benchmark_results.rds \
  Rscript scripts/precompute-vs-alternatives.R
```

The benchmark resets the layout seed before each of five calls and reports the
median elapsed time and interquartile range. Timing covers layout generation
only; scoring and garbage collection are excluded. Every method receives the
same unweighted HMP topology, and the common score is sampled hop-distance
stress. Elapsed times remain machine-dependent; `BENCHMARK_PROVENANCE.md`
records the machine used for the supplied artifact.

To regenerate both weighted-saddle resolutions, export all six-panel figures,
and retain the other benchmark results and their timing metadata, run:

```sh
mkdir -p generated
Rscript scripts/compare-saddle-resolutions.R generated/saddle-resolutions \
  precomputed/vs_alternatives/benchmark_results.rds \
  generated/benchmark_results.rds
```

The command generates 10 by 10 and 15 by 15 saddle meshes from scratch, exporting
a PDF and PNG for each, a combined score CSV, and the complete RDS objects. It
selects 10 by 10 for Supplement S3 and stores both cases and their shared
plotting limits in `weighted_saddle_resolutions`. To export new comparisons
without updating a benchmark RDS, omit the final two arguments.

The helper computes three-dimensional metric MDS and both edge-KK workflows,
reuses each weighted-GRIP matrix for both its refinements, and independently
checks all four score columns by summing the stored fixed paths.
Both edge-KK calls use identical package-default settings, including five
density-continuation stages with at most 50 iterations per stage. The supplied
artifact retains these settings, refinement traces, fixed paths, and check
results within each case. `plot-weighted-saddle.R` is shared by the standalone
exports and the Supplement S3 figure, with common alignment, projection, limits,
panel order, and styling. Each panel, including the target, reports MDS Stress-1
and fixed-path GMDS relative RMSE. The scores use the raw 3D coordinates and all
retained vertex pairs, with separately fitted input-distance scales. For chord
lengths `c` and graph distances `g`, the input scale is
`a = sum(c * g) / sum(g^2)` and Stress-1 is
`sqrt(sum((c - a * g)^2) / sum(c^2))`. This is not the cached
`metric.chord.stress`, which uses `sum((a * g)^2)` as its denominator.
The path label uses the profiled `score.gmds()` result, matching the table.
Small nonzero errors use scientific notation; panel values below `1e-12`
are shown as `< 10^-12` to avoid emphasizing numerical roundoff in the target.

To re-export both figures and their unrounded panel scores from the saved
coordinates, without rerunning the layouts or benchmarks, run:

```sh
Rscript scripts/plot-weighted-saddle.R \
  precomputed/vs_alternatives/benchmark_results.rds generated/saddle-resolutions
Rscript scripts/check-saddle-panel-scores.R \
  precomputed/vs_alternatives/benchmark_results.rds
```

### Three-panel criterion illustration (proposed manuscript figure)

The standalone preview for subsection *Which criterion answers which question*
selects the generating saddle, metric MDS, and weighted GRIP followed by edge-KK
from the same saved comparison. It highlights one retained input-graph shortest
path and its endpoint chord in every panel. The path joins opposite ends of the
grid row closest to, and on the nonnegative side of, `y = 0`; this selection
does not depend on the candidate scores. Scores still use all 4,950 pairs and
the original three-dimensional coordinates, not the highlighted path alone or
the two-dimensional display.

From this directory, generate the preview without recomputing layouts:

```sh
Rscript scripts/plot-saddle-criteria.R \
  precomputed/vs_alternatives/benchmark_results.rds ../build/saddle-criteria
```

The command writes a vector PDF, PNG, proposed caption, retained-path edge list,
and unrounded scores with the highlighted path and chord lengths. The latter
lengths are in each layout's original coordinate units; the reported graph
distance uses median-normalized input edge lengths. They should not be compared
without accounting for scale. The two all-pairs criteria instead fit their
scales separately, as in the full comparison. The script checks the retained
path's endpoints, connectivity, grid row, and input length before plotting.
This preview is not yet inserted in the manuscript; the complete six-panel
comparison is retained in Supplement S3.

The second command checks the panel metrics against an independent least-squares
fit and explicit sums along the stored paths, checks agreement with the cached
scores, and tests invariance to translation, reflection, and uniform scaling.

For a single resolution, source `scripts/weighted-saddle-comparison.R` in R and
call `weighted_saddle_comparison(grid_size = 10L)` or use `15L`. The default is
now 10; `5L` reproduces the former illustrative graph. To update only one
existing case while retaining the others, run
`Rscript scripts/weighted-saddle-comparison.R INPUT.rds OUTPUT.rds`, with an
optional final integer grid size to generate a fresh case. To recheck an existing case independently,
call `check_weighted_saddle(results$weighted_saddle)`, where `results` is the
loaded RDS object. `BENCHMARK_PROVENANCE.md` separates this focused update from
the original long-running benchmark environment.

## Fast article build

The manuscript evaluates the smaller examples during rendering. Its
weighted-saddle comparison, UMB-HMP layout figure, and repeated UMB-HMP runtime
benchmark instead use the supplied results, allowing the article to be
rendered without rerunning the longer computations. Those computations can be
reviewed separately using the inputs, checksums, and scripts above.

The manuscript's hidden `setup` chunk loads the artifact once:

```r
benchmark.results <- read.extdata.rds("vs_alternatives", "benchmark_results.rds")
```

`read.extdata.rds()` and `find.extdata.file()` are defined in that same chunk.
In the submission archive, the lookup first uses
`reproducibility/precomputed/vs_alternatives/benchmark_results.rds`, relative
to the manuscript. Installed-package and repository copies are fallback
locations for development builds. Loading the file does not recompute the
results or validate its checksum; run the `shasum` command above to verify the
supplied artifact before rendering.

## Manuscript-only helpers

The `setup` chunk also defines helpers for rendering and locating inputs:

- `accessible_kable()` and `paper_kable()` format the manuscript tables;
- `find.extdata.file()`, `read.extdata.csv()`, and `read.extdata.rds()` locate
  and read supplied artifacts;
- `plot_layout_panel()` and `plot_common_projection()` draw figure panels;
- `edge.matrix.from.adj()` converts an adjacency list to an edge matrix for
  the manuscript examples.

These helpers are part of the manuscript source, not exported `grip`
functions. Their definitions are included in both `grip-software-paper.Rmd`
and its extracted `grip-software-paper.R` companion. The weighted-saddle figure
instead sources `scripts/plot-weighted-saddle.R`, also included in this archive.
Its similarity alignment uses translation, rotation or reflection, and uniform
scaling solely for display; it is not part of layout generation or scoring.
