# grip 0.2.0

* Removes the 40 long-form compatibility aliases deprecated in the 0.1
  series. See `help("grip-0.2-migration")` for the complete replacement map.
* Narrows the public graph-family API to primitive edge generators,
  configurable mask helpers, and complete graph bundles. Specialized
  edge-list builders, standalone surface embeddings, and parameter-coordinate
  builders remain available internally to those bundles.
* Makes `gmds.result()` internal; public layout methods continue to return the
  common `"grip_gmds_layout"` result format.

# grip 0.1.3

* Unifies topology-first and edge-length-metric layouts under
  `grip(metric = "hop")` and `grip(metric = "edge_length")`.
* Gives `trace.grip()` the same metric selection and removes the short-lived
  `weighted.grip()` and `trace.weighted.grip()` entry points in favor of the
  unified API.
* Documents precisely how positive edge lengths affect force targets,
  shortest paths, multiscale hierarchy construction, neighborhoods,
  normalization, and optional LGKK stages in each metric mode.

# grip 0.1.2

* Makes saved graph-family bundle ordering deterministic when files have tied
  modification times, as can occur on Windows.
* Corrects two undersized adjacency allocations in the compiled triangular
  mesh generator.

# grip 0.1.1

* Adds geodesic-MDS, geodesic-KK, landmark, MISF, weighted, and interactive
  layout workflows developed since the initial submission.
* Caps compiled parallel work at two threads to comply with CRAN resource
  limits, including when automatic thread selection is requested.
* Replaces native assertion aborts with recoverable errors and strengthens
  queue state handling.
* Removes non-package utilities from the installed package and removes
  machine-specific installed paths.
* Shortens graph-family save filenames so they remain portable under long
  Windows temporary-directory paths.
* Expands the pkgdown reference index and reduces README animation assets for
  a smaller source package.
* Declares the package C++17 requirement explicitly.

# grip 0.1.0

* First public release of `grip`.
* Centers the package around two main layout workflows:
  `grip.layout()` for unweighted graphs and
  `grip.layout.weighted()` for weighted graphs.
* Provides layout scoring and real-data candidate comparison via
  `grip.score.layout()` and `grip.compare.layouts()`, together with
  multiscale trace diagnostics via `grip.layout.trace()` and
  `grip.layout.trace.weighted()`.
* Ships a large synthetic graph-family library, bundled weighted real-data
  examples, disconnected-component packing, and plotting helpers.
* Includes advanced public experimental full geodesic-KK (GKK) and landmark
  geodesic-KK (LGKK) preparation, scoring, and refinement helpers for
  weighted-layout evaluation and polish.
* Includes optional Shiny explorers for browsing layout catalogs and synthetic
  graph families.
