# Cross-package benchmark provenance

This record describes the `benchmark_results.rds` artifact used by the
`grip` R Journal software paper. The artifact was regenerated on 2026-08-23 by
`inst/scripts/precompute-vs-alternatives.R` from the `grip` 0.2.0 source tree.

## Reported scaled benchmark

- graph: bundled HMP-only Illumina 16S symmetric-kNN giant component
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

## Hardware and software

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
