# grip 0.2.0.9001 (development version)

- Added opt-in `sparse.metric.mds()` for 2D/3D sparse SGD with pivot-based
  distance preparation, asymmetric endpoint updates, explicit memory checks,
  reproducible sampling, and separately labeled sparse diagnostics.
- Added a sparse/full SGD and weighted-GRIP vignette with 3D examples, repeated
  runs, independent distance-error scores, and measured time and peak memory.
- Avoided repeated copying when collecting edge lengths for edge-only preparation.

- Added `metric.mds(pair_weights = "inverse_squared")` for both SGD and SMACOF,
  using the Zheng et al. (2018) graph-drawing stress weights. Uniform weighting
  remains the default. Fitting, start selection, scaling, and reported stress
  consistently use the selected weights; inverse-squared weighting requires
  positive off-diagonal distances.
- Expanded the SGD–SMACOF vignette with a separate matched comparison of both
  weightings on three 3D graph examples, reporting both error measures and timing.

- Static plot controls now override defaults consistently, including labels,
  limits, aspect ratio, and orthographic vertex symbols/sizes.
- Added `layout.coords()` for unchanged coordinate extraction and
  `estimate.preparation()` for dense-storage lower bounds. Large dense
  preparations warn before searches; the advisory threshold is configurable.
- Native layout and refinement loops support cancellation at safe main-thread
  boundaries. Interrupted construction releases its owned temporary storage.
- Undirected adjacency now requires reciprocal entries, equal multiplicities,
  and matching lengths. Self-loops now fail explicitly instead of being ignored
  or treated differently across representations; remove them before calling.
  Existing parallel-edge handling is preserved and documented in the guide.
- Replaced bending load-order rebinding with explicit base implementations,
  consolidated graph validation, and added a developer map. Numerical objectives
  and defaults are unchanged.
- Full-dependency checks require all tests to execute on Linux, Windows, and
  macOS; a local runner isolates split-library dependencies. The comparison
  article now identifies frozen versions and the distinction between repeated
  timings and independent layout seeds.

## Documentation

* Interactive vignette figures now show their controls continuously, without
  a collapsible "View controls" panel.

* Vignette legends now expand to fit their entries. Clarified the gallery
  Procrustes alignment and verified it with known transformations and the
  saved saddle layouts.
* Simplified the graph gallery to one interactive view per graph, with a
  layout selector for metric-MDS, metric-MDS + edge-KK, or both. Removed
  duplicate static panels and the legend-table disclosure from vignette views.
* Expanded the installed "Graph examples" vignette with interactive 3D
  metric-MDS and edge-KK views of 13 graphs. Added six attributed SuiteSparse
  datasets as `zheng.graphs`, including original matrices and character labels.

* Distinguished CRAN and development installations, including offline-vignette
  building, and added a development notice and direct quick-start navigation.
* Repaired Markdown function links in help and made package help lead to a
  small runnable workflow and the existing installed guides.
* Connected the visual introduction to showcase recipes and a rendered first
  example. Website animations now start with still previews and explicit
  playback controls; guide figures have alternative text.

* Added two installed vignettes: a task-oriented function guide with checked
  export and method coverage, and synthetic graph-family layout examples.
  Linked them from the package overview, README, and website navigation;
  clarified metric, reference-scoring, and projection choices in related guides.

## Breaking changes

* `metric.mds()` now selects native SGD by default. Request
  `backend = "smacof"` to preserve a previous SMACOF analysis, including calls
  that supply `eps`. Explicit `init = "metric_mds"` in `edge.kk()` and
  `kernel.gram.gkk()` follows the new default; their own default initializer
  remains classical MDS. The iteration budget and stopping rules are unchanged.

* Graph-input APIs now reject fractional, nonfinite, or out-of-range vertex
  indices and counts before conversion. Supply either `edges`/`edge_weights`
  or `adj_list`/`weight_list`, not both. A `prepared` graph cannot be combined
  with raw graph inputs; an explicit `n` must match its stored vertex count.
  Previously, some of these inputs were silently truncated or ignored.

* `metric.mds()` now minimizes unweighted raw distance stress using native SGD
  or the optional smacof package. To retain the previous classical-scaling behavior,
  use `classical.mds()`, which also accepts the `add` and `eig` arguments.
  See `help("grip-mds-migration")` for migration instructions.
* In `edge.kk()` and `kernel.gram.gkk()`, `init = "metric_mds"` now requests
  stress minimization. Use `init = "classical_mds"` to retain the previous
  behavior. Both functions still use classical scaling by default.

## New features

* Added the native `metric.mds(backend = "sgd")` backend for uniform
  all-pairs stress and made it the default. SGD records reproducible seeds,
  checkpoint stress, pair-update counts, and failed starts, and reports schedule
  completion without claiming convergence. Its tuning defaults are provisional.
  No Python runtime is required. Use `coords = fit$coords` to pass an SGD result
  to `edge.kk()`.

* Layout printing shows the fitted objective, available termination and start
  information, and separately labeled graph diagnostics without recomputing
  them. The score help explains normalization, units, and unavailable values.

* `metric.mds()` supports multiple starts and reports stress and convergence
  diagnostics, with coordinates returned in the input distance units.
* `score.coordinates()` measures coordinate error against a reference, with
  rigid alignment, similarity alignment, or no alignment.
* `score.surface()` measures area-weighted symmetric RMS distance between
  surfaces.

## Bug fixes

* Fixed `metric.mds()` rejecting valid undirected graphs with nearly equal
  shortest paths. Recreate previously saved graph preparations to obtain
  corrected symmetric distance matrices.

# grip 0.2.0

## Breaking changes

* Removed the 40 long-form compatibility aliases deprecated in the 0.1
  series. See `help("grip-0.2-migration")` for the complete replacement map.
* Restricted the public graph-family API to primitive edge generators,
  configurable mask helpers, and complete graph bundles. Specialized
  edge-list builders, standalone surface embeddings, and parameter-coordinate
  builders remain available internally to those bundles.
* Made `gmds.result()` internal; public layout methods continue to return the
  common `"grip_gmds_layout"` result format.

## Documentation

* Added executable examples for every exported function.

# grip 0.1.3

* Unified topology-first and edge-length-metric layouts under
  `grip(metric = "hop")` and `grip(metric = "edge_length")`.
* Added the same metric selection to `trace.grip()`. Replaced `weighted.grip()`
  and `trace.weighted.grip()` with `grip(metric = "edge_length")` and
  `trace.grip(metric = "edge_length")`, respectively.
* Documented how edge lengths affect layout construction and refinement in
  each metric mode.

# grip 0.1.2

* Made saved graph-family bundle ordering deterministic when files have tied
  modification times, as can occur on Windows.
* Fixed two undersized adjacency allocations in the compiled triangular
  mesh generator.

# grip 0.1.1

* Added geodesic-MDS, geodesic-KK, landmark, MISF, weighted, and interactive
  layout workflows.
* Limited compiled parallel work to two threads to comply with CRAN resource
  limits, including when automatic thread selection is requested.
* Replaced native assertion aborts with recoverable R errors.
* Shortened graph-family save filenames so they remain portable under long
  Windows temporary-directory paths.
* Expanded the pkgdown reference index.
* Declared the C++17 compiler requirement explicitly.

# grip 0.1.0

* First public release of `grip`.
* Introduced two main layout workflows:
  `grip.layout()` for unweighted graphs and
  `grip.layout.weighted()` for weighted graphs.
* Added layout scoring and real-data candidate comparison via
  `grip.score.layout()` and `grip.compare.layouts()`, together with
  multiscale trace diagnostics via `grip.layout.trace()` and
  `grip.layout.trace.weighted()`.
* Included synthetic graph families, weighted real-data examples,
  disconnected-component packing, and plotting helpers.
* Added experimental full geodesic-KK (GKK) and landmark geodesic-KK (LGKK)
  preparation, scoring, and refinement helpers for weighted layouts.
* Included optional Shiny explorers for browsing layout catalogs and synthetic
  graph families.
