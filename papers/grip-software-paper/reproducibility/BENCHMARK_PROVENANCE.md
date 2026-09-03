# Cross-package benchmark provenance

This record describes the `benchmark_results.rds` artifact used by the
`grip` R Journal software paper. The artifact was originally regenerated on
2026-08-23 by `inst/scripts/precompute-vs-alternatives.R` from the `grip` 0.2.0
source tree. Its weighted-saddle component was extended on 2026-09-02, and a
paired mesh-resolution component was added, as documented below. All other
components, including the UMB-HMP timings and
their original software and hardware metadata, are unchanged.

## Reported scaled benchmark

- graph: bundled UMB-HMP-only Illumina 16S symmetric-kNN giant component
- vertices: 4,391
- edges: 9,067
- methods: `igraph::layout_with_fr()`, `igraph::layout_with_drl()`,
  `graphlayouts::layout_with_stress()`, and default hop-metric `grip::grip()`
- repeats per method: 5
- random-number policy: reset to seed 1 before every layout call
- timing statistic: median elapsed wall time, with interquartile range retained
- timing boundary: the layout function call only; graph construction, garbage
  collection, and all `score.layout()` calculations are excluded
- graph metric: unweighted topology for every method; PCA-space edge lengths
  are not supplied as layout weights
- scoring: hop-distance stress on coordinates from the first repeat, with fixed
  scoring seeds 42, 2,000 sampled pairs for stress, and up to 5,000 sampled
  nonedges
- exact edge crossings: disabled for this graph

Median elapsed times were 0.115 seconds for Fruchterman--Reingold, 2.840
seconds for DrL, 24.530 seconds for stress majorization, and 0.531 seconds for
default `grip`. These values characterize this hardware and software stack;
they are not universal performance constants.

## Original benchmark hardware and software

- hardware model: Apple Mac16,5
- processor: Apple M4 Max
- memory: 64 GiB
- operating system: macOS/Darwin 25.6.0, arm64
- R: R 4.6.1 (2026-06-24), "Happy Hop" release
- R platform: aarch64-apple-darwin23
- BLAS: Apple R framework `<R_HOME>/lib/libRblas.0.dylib`
- package versions: `grip` 0.2.0, `igraph` 2.3.3,
  `graphlayouts` 1.2.5
- runtime source: official CRAN macOS arm64 R 4.6.1 installer, verified against
  published SHA-1 `fc9f4ada15589e8e037b9bf05563d21e97181635`
- package installation: `grip` installed from the 0.2.0 source tree into an
  isolated library; `igraph` and `graphlayouts` installed as CRAN binaries for
  R 4.6

## Other artifact contents

The same RDS file retains layouts and common scores for the Zachary karate
club graph, a 12 by 12 mesh, a level-4 Sierpinski carpet, and a weighted
saddle-surface mesh. These retained internal results are not presented as
additional real-data evidence. The paper uses the weighted saddle results for
its external weighted comparison. The single-run carpet timing retained for
historical comparison is not reported as repeated benchmark evidence.

## Weighted-saddle comparison update, 2026-09-02

`scripts/compare-saddle-resolutions.R` compares seven layouts at each of two
resolutions of the same saddle surface. The 10 by 10 mesh has 100 vertices,
180 edges, and 4,950 unordered pairs; the 15 by 15 mesh has 225 vertices,
420 edges, and 25,200 pairs. All 14 layouts were generated and scored, with no
failed or excluded cases. The main paper uses the 10 by 10 mesh for readability
in a six-panel figure. Both complete cases are retained in
`weighted_saddle_resolutions$cases`; `weighted_saddle` is the selected 10 by 10
case. Figure 8 now includes unrefined metric MDS; LGKK and combinatorial GRIP
remain in the table only. All coordinate matrices have three columns.

- Surface: `mesh.surface.graph(size, size, surface = "saddle", amplitude = 0.8)`
  with the same parameter domain `[-1, 1]^2`, orthogonal connectivity, and
  median-normalized edge lengths. Resolution changes the sampled graph, not the
  continuous surface parameters. The even 10-point grid does not sample the
  central parameter lines exactly; the odd 15-point grid does.
- Initial layouts: recomputed from scratch for each graph. GRIP uses seed 1 and
  the mesh preset, with hop and edge-length metrics respectively; igraph KK
  uses edge weights and seed 1. No baseline coordinates from the former 5 by 5
  illustration are reused.
- Metric MDS: `metric.mds(prepared = prepared, dim = 3)` on the fixed weighted
  graph's distance matrix, not on distances between the generating coordinates.
- Edge-KK: the same settings for the metric-MDS and weighted-GRIP initializers:
  density stiffness with identity transform, continuation mix values
  `c(0, 0.25, 0.5, 0.75, 1)`, at most 50 iterations per stage, profiled scale,
  gradient tolerance `1e-8`, and the C++ engine. The remaining settings and
  iteration traces are recorded in `weighted_saddle$comparison_metadata` and
  `weighted_saddle$refinement_traces`. These are the package defaults; settings
  were not tuned separately for the two initializers.
- LGKK: six rounds with 20 local neighbors and eight landmarks, starting from
  exactly the same weighted-GRIP matrix as its edge-KK refinement.
- Scoring: a common `prepare.geodesic.kk()` object within each resolution retains
  a single shortest path for every unordered vertex pair. All seven layouts
  use that same object. The paths necessarily differ between resolutions
  because the graphs differ. Unweighted edge and fixed-path relative RMSE use
  separate least-squares profiled scales; the two gKK diagnostics use their
  common weighted least-squares scale.
- Independent checks: direct edge-length and fixed-path summation, followed by
  explicit scale and residual calculations, agree with all four reported
  diagnostics to within `1e-12` at both resolutions (observed maximum differences
  `2.78e-17` and `5.55e-17`, respectively). The checks also
  verify graph-path lengths and the generating geometry's edge fidelity.
  Results are saved in `weighted_saddle$validation`.
- Display: one similarity alignment per layout to its generating coordinates,
  followed by the same orthographic view (azimuth 35, elevation 22). Shared
  plotting limits are calculated across both resolutions and all six panels,
  and retained as `weighted_saddle_resolutions$display_limits`. The middle and
  right columns show each initializer before and after edge-KK.
- Update environment: R Under development (unstable), 2026-06-24 r90190,
  platform `aarch64-apple-darwin23`, `grip` 0.2.0, `igraph` 2.3.3. The exact
  generation time and versions are retained in the component metadata.

The edge-KK results under the same settings are:

| Mesh | Initializer | Edge relative RMSE | Fixed-path relative RMSE |
|---|---|---:|---:|
| 10 by 10 | Weighted GRIP | 4.21e-9 | 1.60e-9 |
| 10 by 10 | Metric MDS | 5.66e-7 | 4.92e-7 |
| 15 by 15 | Weighted GRIP | 9.57e-6 | 5.43e-6 |
| 15 by 15 | Metric MDS | 1.67e-5 | 1.44e-5 |

Both resolutions show a reduction from each initializer, but the denser graph
retains larger residuals under the same maximum iteration budget. These are
numerical optimization results, not exact-zero claims or recovery of the
saddle's extrinsic curvature. The former 5 by 5 example's below-1e-9 statement
has been removed from the manuscript.

The full generator, `scripts/precompute-vs-alternatives.R`, uses the same
numerical and plotting helpers. The focused command in `README.md` updates only
the two saddle-related components and checks by exact R-object comparison that
all other benchmark components remain unchanged. This comparison does not
measure runtime or match iteration budgets between edge-KK and the six-round
LGKK illustration.

## Metadata portability update (2026-09-03)

Local DESCRIPTION-file path attributes were removed from saved benchmark
session metadata. Every other R-object field, including all numerical results,
coordinates, timings, and software versions, was retained unchanged. The shared
`scripts/portable-session-info.R` helper prevents reintroducing those attributes.
