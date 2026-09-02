# Reproducing the `grip` software-paper results

This directory makes the precomputed components of the paper inspectable and
reproduces all reported results, including reconstruction of the UMB-HMP-only graph
from the supplied minimal upstream tables. The manuscript evaluates the smaller
examples directly.

## Requirements

- R and `grip` 0.2.0 or later;
- `dgraphs` 0.2.0 or later for the manuscript examples and to rebuild the
  UMB-HMP graph from the upstream count and metadata tables;
- `igraph` and `graphlayouts` for the cross-package benchmark;
- the packages listed in the manuscript-level `_Rpackages.txt`.

The article describes `grip` 0.2.0 and uses `dgraphs` 0.2.0 in its examples.
The corresponding source-release tags resolve to the following commits:

- [`grip` v0.2.0](https://github.com/pgajer/grip/tree/v0.2.0):
  [`b3b532a51b3cf8be400203270c4a3a217382878e`](https://github.com/pgajer/grip/tree/b3b532a51b3cf8be400203270c4a3a217382878e);
- [`dgraphs` v0.2.0](https://github.com/pgajer/dgraphs/tree/v0.2.0):
  [`8733d2a74dc09d57b453b88eff3119610c6440f3`](https://github.com/pgajer/dgraphs/tree/8733d2a74dc09d57b453b88eff3119610c6440f3).

The full commit identifiers pin the source trees independently of the tag
names. The environment used to generate the precomputed benchmark results is
recorded separately in `BENCHMARK_PROVENANCE.md` and the RDS metadata.

## Included artifacts

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
| Fixed weighted graph comparison | `precomputed/vs_alternatives/benchmark_results.rds`, components `weighted_saddle` and `weighted_saddle_resolutions` | `scripts/compare-saddle-resolutions.R` and its numerical/plotting helpers, also used by the full benchmark |
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
selects 10 by 10 for the main manuscript and stores both cases and their shared
plotting limits in `weighted_saddle_resolutions`. To export new comparisons
without updating a benchmark RDS, omit the final two arguments.

The helper computes three-dimensional metric MDS and both edge-KK workflows,
reuses each weighted-GRIP matrix for both its refinements, and independently
checks all four score columns by summing the stored fixed paths.
Both edge-KK calls use identical package-default settings, including five
density-continuation stages with at most 50 iterations per stage. The supplied
artifact retains these settings, refinement traces, fixed paths, and check
results within each case. `plot-weighted-saddle.R` is shared by the standalone
exports and the manuscript figure, with common alignment, projection, limits,
panel order, and styling. Small errors use scientific notation
in the table so they are not presented as exact zero.

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
