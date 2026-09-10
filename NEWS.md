# grip 0.2.0.9001 (development version)

## Breaking changes

* `metric.mds()` now minimizes unweighted raw distance stress using the
  optional smacof package. To retain the previous classical-scaling behavior,
  use `classical.mds()`, which also accepts the `add` and `eig` arguments.
  See `help("grip-mds-migration")` for migration instructions.
* In `edge.kk()` and `kernel.gram.gkk()`, `init = "metric_mds"` now requests
  stress minimization. Use `init = "classical_mds"` to retain the previous
  behavior. Both functions still use classical scaling by default.

## New features

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
