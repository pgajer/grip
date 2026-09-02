# Random-saddle flattening experiment

This standalone experiment is not included in the manuscript. Its purpose is to
measure achieved two-dimensional graph-path and edge-length distortion after
replacing the orthogonal saddle grid with symmetric nearest-neighbor graphs of
independent, surface-area-uniform samples. No existing paper benchmark is changed.

## Protocol fixed before repeated runs

- Surface: `(x, y, 0.8*(x^2-y^2))`, with `(x,y)` in `[-1,1]^2`, without noise.
- Sampling: uniform **surface area**, using rejection sampling in the parameter
  square with acceptance probability
  `sqrt(1 + 4*0.8^2*(x^2+y^2))/sqrt(1 + 8*0.8^2)`.
- Independent sample seed: `1000000 + 1000*n + replicate`. R uses the
  Mersenne-Twister generator, Inversion normals, and Rejection sampling.
- Pilot: 20 samples at `n=250`; one `n=500` timing case. The larger study uses
  100 samples at each of `n=250` and `n=500`, reusing pilot cases without changes
  to sampling, optimizer budgets, or selection rules.
- Graph: `dgraphs::create.sknn.graph()`, exact neighbors, symmetric union rule,
  Euclidean lengths in the original 3D data, no pruning and no MST repair.
  The smallest connected `k` is found by checking every integer from 1 upward.
  Each sample is evaluated at `k_conn`, `k_conn+2`, and `2*k_conn`. If two rules
  coincide, the identical graph is evaluated once and reported under both rules.
- Edge lengths are divided by their median. One deterministic input shortest
  path is retained for every unordered pair by `grip::prepare.geodesic.kk()`.
  That same family is used for every candidate for a given graph.
- All fitted candidates have **two coordinate columns**. Initializations are
  metric MDS (`grip::metric.mds()`, classical scaling of graph distances), PCA of
  the original 3D observations, and two seeded weighted-GRIP configurations.
  PCA is an informed initialization that uses the supplied observations, not a
  graph-only method. Its inclusion helps search for a low-distortion realization;
  this is not a timing or graph-only algorithm competition.
- Each initialization receives the package's five-stage density-continuation
  edge-KK refinement, with at most 50 iterations per stage. These four results,
  plus unrefined MDS and PCA, start six direct fixed-path optimizations.
- Direct refinement uses L-BFGS-B, at most 200 iterations, `factr=1e5`, and
  `pgtol=1e-8`. No repulsion, crossing, or angular penalty is added. Drawings with
  crossings or folds are admissible; this is not a test of graph planarity.
- Candidate selection: retain the configuration with the lowest exact profiled
  path relative RMSE among all 14 initial, edge-refined, and path-refined
  candidates. Also report the edge-KK candidate with the smallest edge error.
  A solver reaching its iteration budget is a usable finite candidate, not a
  converged solution; solver statuses are reported separately.
- Representative figure selection: for each sample size, select the sample
  closest to the median best-found path relative RMSE at `k_conn`, breaking ties
  by replicate number (distances within `1e-14` are treated as numerical ties).
  This choice is made mechanically, not by appearance.

## Scores and the directly optimized objective

For observed lengths `y` and input lengths `g`, the unweighted scale is
`a = sum(y*g)/sum(g^2)`, and relative RMSE is
`sqrt(sum((y-a*g)^2)/sum((a*g)^2))`. Edge and path statistics fit separate scales.
The path optimizer minimizes the square of exactly this statistic. For
`q=sum(y^2)`, `t=sum(y*g)`, and `r=sum(g^2)`, its derivative with respect to `y`
is `2*r/t^2 * (y - (q/t)*g)`. The saved path-to-edge incidence lists distribute
this derivative to edge lengths, then to their endpoint coordinates.
`objective.cpp` is a standalone evaluator, not a change to the package API.
All reported scores use unstabilized Euclidean lengths (`epsilon=0`).

MDS Stress-1 instead uses the configuration-distance denominator:
`sqrt(sum((chords-a*g)^2)/sum(chords^2))`. Scores use all unordered pairs, not a
sample of pairs or projected display coordinates. Numerical values are
dimensionless; plots show errors as percentages where indicated.

The best achieved path error is an **upper bound** on the optimum over 2D
configurations. Positive optimization residuals alone do not prove an exact
realization impossible. Small aggregate path error can coexist with much larger
local edge error away from zero.

## Independent geometric obstruction and validation

Every four-vertex clique supplies six specified edge lengths. Double centering
their squared distance matrix gives a Gram matrix. A positive third eigenvalue
certifies that these six lengths require three dimensions. The numerical check
requires the third/first eigenvalue ratio to exceed `1e-8`; the strongest witness
and its actual distances are saved. This is a conservative floating-point
obstruction check, not an interval-arithmetic proof. Its absence is inconclusive.

Because the supplied weights are Euclidean distances between the original
observations, each edge is a shortest route between its endpoints. A 2D zero-loss
solution for all pairs would therefore have to preserve every clique edge under
the common positive scale. The obstruction also excludes exact zero all-pairs
fixed-path loss. It does not quantify the minimum average error or establish
continuous-surface geodesic accuracy.

Validation includes analytic-gradient finite differences, scaling invariance,
comparison against `grip::score.gmds(..., edge_length_epsilon=0)`, a known planar
zero-loss control, a known tetrahedral obstruction, and a sampling-distribution
check. Every sampled graph also checks all retained input path lengths and the
known zero-loss original 3D configuration. The selected winner's path score is
checked against the package scorer. Completed cases are saved atomically and
resumed only when their protocol matches.

## Reproduction

From the manuscript directory, with `grip` and `dgraphs` 0.2.0, `igraph`, and
`Rcpp` installed:

```sh
Rscript reproducibility/experiments/random-saddle/experiment.R checks build/random-saddle
Rscript reproducibility/experiments/random-saddle/experiment.R run build/random-saddle 250 20 4
Rscript reproducibility/experiments/random-saddle/experiment.R run build/random-saddle 250,500 100 4
Rscript reproducibility/experiments/random-saddle/plot-results.R build/random-saddle
Rscript reproducibility/experiments/random-saddle/check-results.R build/random-saddle
Rscript reproducibility/experiments/random-saddle/plot-results.R build/random-saddle
Rscript reproducibility/experiments/random-saddle/build-report.R build/random-saddle
```

`GRIP_RJOURNAL_PACKAGE_LIBRARY` optionally prepends a package library. Runs record
their actual R/package versions and start/end times. Parallel workers distribute
independent samples; per-graph elapsed times are wall times within workers, not
serial-equivalent CPU times. The runner preserves every candidate, exact input
graph, retained path, score, witness, and solver status in per-sample RDS files.
Figures, CSV summaries, and a Markdown results record are generated separately
from those saved objects. Changing the plotting code never reruns the layouts.

The post-run checker independently recomputes graph shortest-path distances
using `igraph`, checks that all six edges of each obstruction witness are retained
as one-edge paths, and recomputes winner scores by summing saved paths in R. It also
extends the winning candidate for the median-selected sample within each of the
six sample-size/neighborhood-rule groups by up to 2,000 additional L-BFGS-B
iterations. This is a separately reported budget-sensitivity diagnostic: the
original repeated-study results and figure selection are not replaced by these
six more extensively optimized examples.

The package's weighted shortest-path code treats distances within
`sqrt(.Machine$double.eps)` relative tolerance as near-ties. The independent
checker reports discrepancies against `igraph` within this tolerance, their
frequency, and the resulting score change when the stricter distances replace
the package's input distances while the retained routes remain fixed. It does
not equate tolerance-level agreement with bitwise equality.
It separately counts edge-endpoint pairs for which near-tie handling selected
an indirect route. These do not invalidate the obstruction when all six edges
of the independently verified witness retain their direct one-edge routes.

The optional final command requires `pdflatex` and builds `report.pdf`, combining
the methods, summary tables, validation and optimization-budget accounting, and
the four full-size figure pages. It requires the complete 100-sample-per-size
study and the extended-budget diagnostic; it will not silently label partial
results as complete. Figure-only PDF and PNG exports do not require LaTeX.

## Interpretation boundary

This experiment measures graph-to-coordinate fidelity conditional on the
constructed graphs. Connectivity is not a guarantee of accurate surface
geodesics. Increasing `k` both adds geometric constraints and changes the
shortest-path family. Results across `k` are sensitivity analyses, not estimates
against one unchanging reference metric. No continuous-surface curvature or
geodesic-recovery claim follows from a low graph-path score alone.
