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
