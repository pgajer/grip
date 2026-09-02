# Cross-package benchmark provenance

This record describes the `benchmark_results.rds` artifact used by the
`grip` R Journal software paper. The artifact was originally regenerated on
2026-08-23 by `inst/scripts/precompute-vs-alternatives.R` from the `grip` 0.2.0
source tree. Only its weighted-saddle component was extended on 2026-09-02,
as documented below; all other components, including the UMB-HMP timings and
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

`scripts/weighted-saddle-comparison.R` extends the comparison to seven layouts
on the same 25-vertex, 40-edge saddle graph. Figure 8 shows six panels: the
generating geometry plus five layouts; combinatorial GRIP and unrefined metric
MDS are retained in the score table only. All coordinate matrices have three
columns.

- Retained baselines: the original combinatorial GRIP, weighted GRIP, and
  weighted igraph KK coordinates. A fresh call to the shared generator
  reproduced all seven coordinate matrices exactly on the update environment.
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
  the exact retained weighted-GRIP matrix. This reproduces the previous LGKK
  matrix exactly, with maximum coordinate difference zero.
- Scoring: one common `prepare.geodesic.kk()` object retains a single shortest
  path for each of 300 unordered vertex pairs. The object is saved as
  `weighted_saddle$prepared`. Unweighted edge and fixed-path relative RMSE use
  separate least-squares profiled scales; the two gKK diagnostics use their
  common weighted least-squares scale.
- Independent checks: direct edge-length and fixed-path summation, followed by
  explicit scale and residual calculations, agree with all four reported
  diagnostics to a maximum absolute difference of `2.78e-17`. The checks also
  verify graph-path lengths and the generating geometry's edge fidelity.
  Results are saved in `weighted_saddle$validation`.
- Result: weighted-GRIP plus edge-KK has edge/path relative RMSE
  `4.01e-11` / `2.48e-11`; metric-MDS plus edge-KK has
  `5.09e-10` / `1.34e-10`. These are numerical residuals, not assertions of
  mathematically exact zero or recovery of the saddle's extrinsic curvature.
- Update environment: R Under development (unstable), 2026-06-24 r90190,
  platform `aarch64-apple-darwin23`, `grip` 0.2.0, `igraph` 2.3.3. The exact
  generation time and versions are retained in the component metadata.

The full generator, `scripts/precompute-vs-alternatives.R`, sources the same
helper. The focused command in `README.md` updates only this component and
checks by exact R-object comparison that the other benchmark components remain
unchanged. This comparison does not measure runtime or match iteration budgets
between edge-KK and the six-round LGKK illustration.
