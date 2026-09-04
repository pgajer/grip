# Reference-based scoring of the single saddle example

Run from the repository root, using the development source of grip:

```sh
Rscript papers/grip-software-paper/reproducibility/experiments/two-fidelity-pilot/run-saddle-reference.R output/single-saddle-reference
Rscript papers/grip-software-paper/reproducibility/experiments/two-fidelity-pilot/check-saddle-reference.R output/single-saddle-reference
```

An optional second argument changes the sample size (default 500). The run uses
one area-uniform saddle sample with seed 1, C = 0.8, and parameter domain
[-1, 1]^2. It fits metric MDS and metric MDS followed by edge-KK for symmetric
kNN graphs with k = 3,...,20 and component-MST repair. This is a new single-cloud
comparison, not a recomputation of the five-cloud Figure 7 experiment.

The output directory contains the fitted coordinates, a CSV of all scores,
the reference mesh, a sampling/discretization check at k = 10, an environment
record, and a PNG comparison. Existing layout caches are reused only when their
recorded settings match. Use a new output directory to refit with changed code.
The runner does not install grip or modify the manuscript.

## Scoring conventions

* Graph scores are evaluated before coordinate alignment. The path and edge
  relative RMSE values use independently profiled scales, as in Figure 7.
  All vertex pairs are used, with one retained deterministic input shortest
  path per pair and no edge-length stabilizer. MDS Stress-1 uses the
  configuration-distance normalization, not the graph-target normalization.
* Coordinate RMSE is the square root of the mean squared Euclidean vertex
  displacement after alignment. Relative coordinate RMSE divides this by the
  reference RMS radius. Rigid alignment fits translation and rotation/reflection
  but preserves scale; similarity alignment additionally fits a uniform scale.
  Correspondence is supplied by the observation row order.
* Surface RMS is estimated from independent area-uniform samples on each mesh,
  with exact closest-point distances to the other mesh's triangles. It is the
  square root of the equally weighted mean of the two directional mean squared
  distances, in coordinate units. No second surface registration is performed.
  Zero-area or overlapping faces require the interpretation stated in
  `help("score.surface")`; the score is not a test for folds or mesh validity.
* The embedded mesh uses the original parameter-space Delaunay connectivity.
  The reference subdivides those same parameter triangles and lifts the new
  vertices onto the saddle. Thus both meshes represent the same parameter
  footprint; neither fills the unsampled square corners. There are two
  subdivision rounds in the main comparison and three in the k = 10 check.
* The main run uses 2,000 and 8,000 samples per direction. The k = 10 check uses
  8,000 and 32,000 on the denser reference. Monte Carlo standard errors describe
  surface sampling only, not uncertainty across independently sampled clouds.
  The original saddle sample is scored as a discretization baseline; its
  piecewise-planar surface need not agree exactly with the denser reference.

For interactive work with existing `X`, `graphs`, `mds`, and `mds.edge.kk`, source
`score-saddle-reference.R` and call `score.saddle.reference()`. This avoids
refitting the layouts. The added section of `single-saddle-ivue.R` demonstrates
that call. `plot-saddle-reference.R` plots the returned scores without rerunning
the calculations.
