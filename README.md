
<!-- README.md is generated from README.Rmd. Please edit that file -->

# grip

**grip** (Graph dRawing with Intelligent Placement) is an R package for
multiscale graph layout. Its primary unweighted and weighted workflows
target 2D and 3D, while opt-in weighted-GRIP, metric-MDS, and edge-KK
workflows also support higher-dimensional embeddings. The main workflow
is:

- `grip(metric = "hop")` for topology-first layouts,
- `grip(metric = "edge_length")` when edge lengths define the graph
  metric,
- `compare.layouts()` and `score.layout()` for real-data layout
  selection,
- `trace.grip()` with the same metric choice for diagnostics.

The package also includes advanced public experimental geodesic-KK
utilities for weighted-layout scoring and polish. It builds on the GRIP
method described in [Gajer & Kobourov
(2002)](https://doi.org/10.7155/jgaa.00052) and [Gajer, Goodrich &
Kobourov (2004)](https://doi.org/10.1016/j.comgeo.2004.03.014).

## Installation

``` r
# Install from GitHub
install.packages("remotes")
remotes::install_github("pgajer/grip")
```

## Quick start

``` r
library(grip)

# Lay out a small mesh in 2D using the "mesh" preset
edges <- edges.mesh(8, 8)
coords <- grip(edges, n = 64, dim = 2, preset = "mesh", seed = 1)
plot.layout(coords, edges, pch = 16, cex = 0.6, main = "8x8 mesh")
```

## Features

- Multiscale force-directed layout in 2D and 3D via C++ (Rcpp).
- A unified `grip()` interface for hop-metric and edge-length-metric
  layouts.
- Opt-in multiscale weighted layout in dimensions greater than 3 via
  `weighted.grip.nd()`, with higher-dimensional metric-MDS and edge-KK
  workflows available through `classical.mds()`, `metric.mds()`, and
  `edge.kk()`. `classical.mds()` uses classical scaling; `metric.mds()`
  minimizes raw distance stress through the optional **smacof** package.
  In versions through 0.2.0, `metric.mds()` performed classical scaling:
  use `classical.mds()` to preserve that behavior. Edge-KK defaults to
  `init = "classical_mds"`; explicit `init = "metric_mds"` now requests
  stress minimization.
- Layout comparison and quality scoring across seeds and parameter
  settings (`compare.layouts()`, `score.layout()`).
- Multiscale trace diagnostics for both metrics via `trace.grip()`.
- Advanced public experimental geodesic-KK utilities for weighted-layout
  scoring and polish (`prepare.geodesic.kk()`, `score.geodesic.kk()`,
  `prepare.landmark.geodesic.kk()`, `score.landmark.geodesic.kk()`).
- Synthetic graph-family helpers for benchmark and geometry-rich
  examples.
- Handles disconnected graphs automatically (component packing).
- Tuned presets for common graph families (see table below).
- Static 3D projection for vignettes and reports
  (`plot.layout(projection = "ortho")`, `project.3d()`).

## Presets

| Family | Preset | Tuned on |
|:---|:---|:---|
| Rectangular grid or lattice | `preset = "mesh"` | `8x8` and `12x12` meshes |
| Sierpinski carpet | `preset = "carpet"` | Level 3 and 4 carpets |
| Tree-like graph | `preset = "tree"` | Binary trees, depths 5 and 6 |
| 3D torus or cylinder | `preset = "torus"` | Torus sizes `8x8` through `20x20` |

Presets set sensible defaults for the GRIP parameters. Any explicit
argument you pass overrides the preset value.

## Choosing a workflow

- Start with `grip(metric = "hop")`, the default, when topology should
  define the multiscale hierarchy and graph neighborhoods.
- Use `grip(metric = "edge_length")` when positive edge lengths should
  also define shortest-path distances, hierarchy construction,
  neighborhoods, and insertion anchors.
- Use `compare.layouts()` and `score.layout()` when the graph is
  important enough to justify a candidate shortlist rather than a single
  run.
- Use `trace.grip()` with the corresponding `metric` when you need to
  diagnose how a solve evolved.
- Add GKK/LGKK only after you already have weighted candidate layouts
  and need geodesic-aware scoring or polish; these are advanced public
  experimental tools rather than the default starting point.

The historical argument names `edge_weights` and `weight_list` represent
positive edge **lengths**, not connection strengths. With
`metric = "hop"`, supplied lengths set adjacent-edge force targets while
standard GRIP hierarchy and neighborhood searches still count hops. With
`metric = "edge_length"`, the lengths also define weighted shortest
paths throughout the multiscale engine and are median-normalized by
default. See `?grip` for the complete semantics and normalization
options.

## Gallery

These animations show the multiscale refinement recorded by
`trace.grip()`. Starting from a coarse placement, GRIP introduces
vertices and refines their positions through the final layout.

### Sierpinski Carpet (Level 4)

<img src="https://pgajer.github.io/grip/reference/figures/readme-sierpinski-carpet-level-4-trace.gif" alt="Animated multiscale layout of a level-4 Sierpinski carpet" width="600" />

### Sierpinski Triangle (Level 6)

<img src="https://pgajer.github.io/grip/reference/figures/readme-sierpinski-triangle-level-6-trace.gif" alt="Animated multiscale layout of a level-6 Sierpinski triangle" width="600" />

## More examples

### Graph fidelity and reference-surface agreement

The sampled-saddle example compares the original observations,
metric-MDS, and metric-MDS + edge-KK using the same parameter-plane
triangulation. Graph-path fidelity, corresponding-coordinate agreement,
and surface proximity measure different aspects of the result.

[![Three aligned saddle configurations with a common
triangulation](https://pgajer.github.io/grip/reference/figures/readme-saddle-reference.png)](https://pgajer.github.io/grip/supplements/S4-interactive-saddle.html)

[Rotate the configurations and compare surface
overlays](https://pgajer.github.io/grip/supplements/S4-interactive-saddle.html).
The self-contained [offline
HTML](papers/grip-software-paper/supplement/S4-interactive-saddle.html)
is also included with the manuscript. Fits use grip 0.2.0; additional
reference diagnostics use development commit `b72f61d`, and
visualization uses ivue commit `872f9d4`.

**Edge-list input (2D, circle placement)**

``` r
edges <- edges.cycle(18)
coords <- grip(edges, n = 18, dim = 2, placement = "circle", seed = 2)
plot.layout(coords, edges, pch = 16, cex = 0.7)
```

**Edge-length-metric adjacency list (geometry-aware)**

``` r
adj_list <- list(c(2), c(1, 3), c(2, 4), c(3))
weight_list <- list(c(1.0), c(1.0, 2.0), c(2.0, 1.5), c(1.5))
coords <- grip(
  adj_list = adj_list, weight_list = weight_list,
  metric = "edge_length", n = 4, dim = 2, seed = 12
)
plot.layout(coords)
```

**3D layout with static projection**

``` r
edges <- edges.torus(8, 12)
coords <- grip(edges, n = max(edges), dim = 3, preset = "torus", seed = 3)
plot.layout(coords, edges, projection = "ortho", main = "Torus (8x12)")
```

## Layout comparison

For real-world graphs without a known target layout, `compare.layouts()`
compares candidates across seeds and reports quality metrics. Use
`params.from.summary()` to extract the winning parameters for reuse.

``` r
edges <- edges.mesh(10, 10)
cmp <- compare.layouts(edges, n = 100, dim = 2,
                            candidates = c("default", "mesh"),
                            seeds = 1:3)
cmp$summary[, c("candidate", "score.composite", "sampled.stress.mean")]
```

## Documentation

The package ships with four core vignettes:

- **Getting Started with grip** — the shortest path through the default
  unweighted workflow, with guidance on when to switch to weighted,
  trace, or comparison workflows.
- **Weighted Graph Layouts with grip** — geometry-aware layouts,
  geodesic scoring, and 2D-versus-3D decisions for weighted graphs.
- **Choosing Layouts for Real Data** — a step-by-step workflow using the
  Zachary karate club and Krackhardt kite examples, plus a larger
  weighted HMP/U01 case study.
- **Tracing and Diagnosing Layouts** — frame-by-frame tracing for
  understanding how a solve evolves.

The pkgdown site also includes companion articles such as the
interactive explorer guide, the HMP/U01 object-structure note, the
comparison article, and the synthetic-family gallery.

The geodesic-KK helpers are public and documented in the reference
index, but they are intentionally positioned as advanced experimental
tools layered on top of the main weighted workflow.

## Citation

If you use **grip** in published work, please cite the underlying
algorithm:

> Gajer, P. and Kobourov, S.G. (2002). GRIP: Graph dRawing with
> Intelligent Placement. *Journal of Graph Algorithms and Applications*,
> 6(3), 203–224. doi:
> [10.7155/jgaa.00052](https://doi.org/10.7155/jgaa.00052)

> Gajer, P., Goodrich, M.T. and Kobourov, S.G. (2004). A
> multi-dimensional approach to force-directed layouts of large graphs.
> *Computational Geometry*, 29(1), 3–18. doi:
> [10.1016/j.comgeo.2004.03.014](https://doi.org/10.1016/j.comgeo.2004.03.014)

## License

GPL (\>= 3)
