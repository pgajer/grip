# Non-Oracle Graph-Selection Benchmark Plan

Status: proposed plan  
Owner: Experiment Engineer  
Created: 2026-05-14  
Scope: planning only; no experiments run here

## Purpose

The benchmark should test whether non-oracle graph-selection scores choose
graph parameters that are close to oracle-best parameters on controlled
quadratic-surface examples where intrinsic/surface geodesic truth is known.

The manuscript question is narrow:

```text
Can a practical selector choose a good graph G_theta without seeing
the surface geodesic distance matrix?
```

The benchmark should not become a new regression-over-graphs paper. Its role is
to support the SIMODS bridge:

```text
data -> graph-selection criterion -> weighted graph -> graph geodesic metric
     -> GMDS
```

## Recommended High-Level Design

Use the existing quadratic-surface oracle benchmark lane as the primary target,
but add a selector layer on top of each graph-parameter grid. The main response
for each selector is:

```text
theta_hat(selector, family, dataset, seed)
```

and the main evaluation compares that choice to:

```text
theta_oracle = argmin_theta Err_oracle(G_theta)
```

The first full benchmark should stay 2D until the selector pipeline is stable.
3D quadratic hypersurfaces should enter as a staged follow-up, not the first
full non-oracle run, because the 3D reference layer is more expensive and still
depends on Delaunay-oracle settings.

## 1. Synthetic Examples

### Primary 2D Examples

Use two graph surfaces in a 2D parameter domain:

```text
paraboloid: y = c1 x1^2 + c2 x2^2
saddle:     y = c1 x1^2 - c2 x2^2
```

Use the requested curvature grid:

```text
(c1, c2) in {(1,1), (1,2), (1,4), (2,4)}
```

Use the requested sample-size grid:

```text
n in {100, 200, 400}
```

Use uniform sampling on the parameter disk for the first selector benchmark.
Square domains and radial center/boundary sampling should remain a supplement
or second wave, because the first non-oracle question is selector recovery, not
sampling-domain stress.

Recommended seeds:

```text
smoke: seeds {1, 2}
full:  seeds {1, 2, 3, 4, 5}
```

The full 2D grid is therefore:

```text
2 surfaces x 4 curvatures x 3 sample sizes x 5 seeds = 120 datasets
```

That is large enough to test selector stability without drifting into a huge
domain/sampling benchmark.

### Staged 3D Quadratic Hypersurfaces

Use the existing 3D quadratic form:

```text
y = sum_{i=1}^k c_i x_i^2 - sum_{i=k+1}^3 c_i x_i^2
```

Parameter-domain shapes:

```text
primary staged domain: ball
secondary stress domain: cube
```

Representative index values should be:

```text
k = 0, 1
```

Indices 2 and 3 generate hypersurfaces identical in shape to those of indices
1 and 0 up to coordinate/sign reflection. Existing 3D scripts use `k = 3` for
positive-curvature cases; that is shape-equivalent to the all-negative `k = 0`
case after reflecting the vertical coordinate. For continuity with existing
reports, the implementation can keep existing `k = 3` case labels, but the
manuscript-facing design should describe the canonical shape classes as
index 0 and index 1.

Candidate 3D coefficient vectors:

```text
isotropic:  (1,1,1)
anisotropic: (1,2,4), (1,4,4)
optional hard anisotropic: (1,1,4)
```

The first three are already used in the 3D quadform/Delaunay benchmark scripts
as `c111`, `c124`, and `c144`.

Candidate 3D sample sizes:

```text
n in {80, 120, 200}
```

Recommended staged 3D smoke after the 2D selector works:

```text
domain: ball
index classes: k = 0-equivalent and k = 1
coefficients: (1,1,1), (1,2,4)
n: 80, 120
seeds: 2
families: adaptive-radius, cKNN, sKNN, iKNN
```

Recommended staged 3D full follow-up:

```text
domain: ball first; add cube only if ball results are interpretable
index classes: k = 0-equivalent, k = 1
coefficients: (1,1,1), (1,2,4), (1,4,4)
n: 80, 120, 200
seeds: 3
```

3D should not enter the first full run. It should remain staged until:

- the 2D selector pipeline produces finite scores for all families;
- selector outputs can be joined cleanly to oracle metrics;
- ambiguity rules for near-oracle optima are implemented;
- 3D reference generation settings are fixed, likely using Delaunay
  `edge.length.factor = 4`.

## 2. Graph Families And Parameter Grids

### Primary Families

Use adaptive-radius and continuous-kNN as the main manuscript families.

Adaptive-radius:

```text
k_scale in {3, 4, 5, 6, 8, 10, 12, 15}
radius_rule in {"max", "geomean"}
radius_factor in {1.0, 1.25, 1.5, 1.75}
```

For main-text simplicity, `radius_rule = "max"` can represent adaptive-radius,
while `radius_rule = "geomean"` is better reported as cKNN/continuous-kNN. If
both are generated through `create.adaptive.radius.graph()`, the results table
must still label them as separate graph families.

cKNN:

```text
k_scale in {3, 4, 5, 6, 8, 10, 12, 15}
delta in {1.0, 1.25, 1.5, 1.75}
```

For each family, order parameter settings by increasing graph density/edge
count in addition to the native parameter order. The stability criteria are
consecutive-graph criteria and need a meaningful sequence.

### Structural Baselines

sKNN:

```text
k in {3, 4, 5, 6, 8, 10, 12, 15}
```

iKNN:

```text
k in {3, 4, 5, 6, 8, 10, 12, 15}
prune.method = "none" for direct comparison
```

iKNN is useful because existing `gflow` selection machinery was first developed
around iKNN. However, the benchmark should not let iKNN-specific conveniences
define the common interface.

### Sparse Sentinels

Fixed-radius and mKNN should not be full-grid main competitors. Current project
state treats them as weaker broad performers.

Recommended handling:

- Smoke run: exclude fixed-radius and mKNN.
- Full 2D run: include a small sentinel grid only if runtime is acceptable.
- Main text: mention only if they clarify why adaptive scale is needed.
- Supplement: place their negative/sentinel result there.

Fixed-radius sentinel:

```text
radius_rank in {4, 8, 12}
```

mKNN sentinel:

```text
k in {4, 8, 12}
```

The sentinel purpose is to show that the selector is not merely selecting dense
graphs and to retain continuity with previous graph-family reports.

## 3. Oracle Target

For each dataset, graph family, and parameter theta, build a weighted graph
`G_theta` and compute its graph geodesic distance matrix:

```text
D_G(theta) = graph shortest-path distances on G_theta
```

Use the surface-reference distance matrix as truth:

```text
D_S = numerical surface geodesic distances
```

Define the oracle error:

```text
Err_oracle(G_theta) = summarize.isometry.deviation(D_G(theta), D_S, scale = TRUE)$rel_rms_error
```

Then:

```text
theta_oracle = argmin_theta Err_oracle(G_theta)
```

Use the MST-repaired unpruned lifecycle stage as the primary graph stage:

```text
stage = "raw.repaired" for 2D constructors with lifecycle branches
stage = "final" only where final is exactly the repaired unpruned graph
```

Primary oracle metrics:

- calibrated relative RMS distance error;
- relative geodesic stress;
- relative absolute residual q90 and q95;
- signed bias;
- shortcut fraction;
- Spearman distance correlation;
- number of graph edges;
- number of MST repair edges added.

The headline oracle optimum should be chosen by calibrated relative RMS error.
The q95 residual, bias, and shortcut fraction are required diagnostics for
failure-mode interpretation, not separate oracle objectives unless explicitly
specified in a sensitivity table.

## 4. Non-Oracle Selection Criteria

Every selector should produce a scalar score for each graph parameter setting
and a selected parameter `theta_hat`. Lower scores are better unless stated
otherwise.

### A. Graph Smoothing / GCV

Primary non-oracle selector:

```text
GCV_score(theta) = average_r GCV(y_r; G_theta)
theta_hat_GCV = argmin_theta GCV_score(theta)
```

Preferred mathematical formulation:

```text
f_hat = argmin_f ||y - f||_2^2 + lambda f^T L_G f
```

Equivalently, using graph spectral filtering:

```text
f_hat = V h_lambda(Lambda) V^T y
GCV(lambda) = RSS(lambda) / (n - trace(S_lambda))^2
```

where `S_lambda = V h_lambda(Lambda) V^T`.

Existing `gflow` has graph-regression/GCV machinery based on spectral
smoothing and `refit.rdgraph.regression()`, including per-column GCV helpers.
The missing benchmark piece is a graph-family-general wrapper that accepts a
prebuilt graph object or adjacency/weight lists, computes a Laplacian smoother,
and returns comparable GCV scores for adaptive-radius, cKNN, sKNN, and iKNN.

For the first implementation, use one common smoother across families rather
than each family's historical smoother. The simplest acceptable benchmark
wrapper is:

```text
1. Extract adj_list and weight_list from G_theta.
2. Build a symmetric weighted graph Laplacian.
3. Compute enough eigenpairs, or the dense eigendecomposition for n <= 400.
4. Evaluate a shared lambda/eta grid.
5. Return the minimum GCV for each feature.
6. Average over selected features.
```

Do not compare adaptive-radius GCV from one smoother to iKNN GCV from a
different regression pipeline in the headline analysis.

### B. Degree-Distribution Stabilization

Compute Jensen-Shannon divergence between degree distributions of consecutive
graphs in an ordered parameter sequence:

```text
JS_degree(theta_j) =
  JS(P_degree(G_theta_j), P_degree(G_theta_{j+1}))
```

Select either:

- the first parameter after a stabilization elbow;
- or the smallest parameter whose JS score is within `eps` of the minimum over
  the connected/repaired candidate region.

Recommended first rule:

```text
theta_hat_JS = smallest theta in eligible set with
JS_degree(theta) <= min(JS_degree) * (1 + eps)
eps = 0.05
```

This mirrors the existing iKNN stability-selection style and favors the
sparsest stable graph.

### C. Graph Edit / Edge Symmetric-Difference Stability

Compute edge-set symmetric difference or Jaccard distance between consecutive
graphs:

```text
Edit(theta_j) =
  |E_j symmetric_difference E_{j+1}| / |E_j union E_{j+1}|
```

Use normalized Jaccard distance for cross-family comparability. Raw symmetric
difference counts can be reported as a diagnostic but should not drive
cross-family conclusions.

Selection rule:

```text
theta_hat_edit = smallest eligible theta within 5% of minimum Edit(theta)
```

### D. Connectivity / MST Repair Burden

Track for every graph:

```text
n_components_raw
n_mst_edges_added
mst_edge_fraction = n_mst_edges_added / max(n_edges_raw_repaired, 1)
lcc_fraction_raw
```

This should be a veto/penalty rather than a standalone selector:

- ineligible if the raw graph has too many components after the parameter is
  expected to be in a connected regime;
- penalize large MST repair burden;
- report the first parameter where `n_mst_edges_added = 0` or where the burden
  remains below a small threshold.

Recommended penalty:

```text
repair_penalty(theta) = n_mst_edges_added / max(n_edges_raw_repaired, 1)
eligible if repair_penalty <= 0.02
```

Keep the threshold visible and run a sensitivity check at 0.01 and 0.05.

### E. Optional Layout Stability

Layout stability should be optional and supplement-only unless it is cheap and
fully deterministic.

A possible definition:

```text
1. Run a fixed-seed weighted GRIP layout for theta_j and theta_{j+1}.
2. Procrustes-align layouts.
3. Compute RMS displacement after alignment.
```

Do not use layout stability in the first smoke run. It is downstream of graph
selection and can confuse graph-geodesic quality with layout behavior.

## 5. Synthetic Feature-Signal Design For GCV

GCV requires graph signals. These signals should be generated from the known
parameter/surface coordinates but must not use the oracle distance matrix.

For each sampled point with parameter coordinates `x1, x2` and surface height
`y`, generate a feature matrix with four blocks.

### Block 1: Ambient Coordinates

```text
x1, x2, y
```

These should be used first in smoke because they are simple, deterministic, and
debuggable. They may be too easy, so they should not be the only full-run
feature set.

### Block 2: Smooth Polynomial Signals

Examples:

```text
x1
x2
y
x1 * x2
x1^2 - x2^2
x1^2 + x2^2
sin(pi*x1)
cos(pi*x2)
sin(pi*x1) * cos(pi*x2)
```

These are the primary full-run smooth features. They test whether the graph
supports smooth functions over the sampled geometry.

### Block 3: Noisy Smooth Signals

For each smooth signal `s_r`, add:

```text
s_r_noisy = standardize(s_r) + sigma * epsilon
epsilon ~ N(0, 1)
sigma in {0.05, 0.15}
```

The first full run should use `sigma = 0.10` or `0.15`. Use `0.05` only as a
sensitivity check because it may make GCV too easy.

### Block 4: Distractor / Noise Features

Add independent features:

```text
noise_r ~ N(0, 1)
```

Recommended feature counts:

```text
10 smooth deterministic
20 noisy smooth
30 pure noise distractors
```

Then run top-variable-feature selection:

```text
top 20 and top 30 most variable features
```

First benchmark order:

1. Smoke: ambient coordinates plus 10 smooth polynomial features, no pure noise.
2. First full: 10 smooth deterministic + 20 noisy smooth + 30 pure noise, select
   top 20 variable features.
3. Sensitivity: top 30 variable features and a no-distractor version.

The main GCV result should use top 20 variable features, because that mirrors
the intended practical real-data workflow while keeping runtime controlled.

## 6. Selector Evaluation Metrics

For every `(dataset, seed, family, selector)` report:

- selected parameter `theta_hat`;
- oracle-best parameter `theta_oracle`;
- oracle rank of `theta_hat`, where rank 1 is best oracle error;
- oracle error ratio:

```text
Err_oracle(G_theta_hat) / min_theta Err_oracle(G_theta)
```

- absolute and normalized parameter distance from the oracle optimum;
- edge-count ratio relative to oracle graph:

```text
n_edges(theta_hat) / n_edges(theta_oracle)
```

- whether `theta_hat` is in the near-oracle set:

```text
Err_oracle(G_theta_hat) <= 1.05 * min_theta Err_oracle(G_theta)
Err_oracle(G_theta_hat) <= 1.10 * min_theta Err_oracle(G_theta)
```

- selector stability across seeds:
  - modal selected parameter;
  - median selected parameter;
  - interquartile range;
  - probability of near-oracle selection;
  - median oracle error ratio;
- failure flags:
  - disconnected or high MST burden;
  - non-finite GCV;
  - selector chooses boundary of parameter grid;
  - multiple near-oracle optima;
  - selector score flat or multi-modal;
  - oracle optimum itself is boundary-of-grid.

For 2D single-parameter curves, normalized parameter distance can be:

```text
abs(rank(theta_hat) - rank(theta_oracle)) / (n_candidates - 1)
```

For adaptive-radius and cKNN multi-parameter grids, use:

```text
native grid rank distance in density-ordered sequence
and
componentwise normalized distance over (k_scale, radius_factor/delta)
```

Do not overinterpret componentwise distance when two different parameter
combinations produce similar edge counts and near-identical oracle errors.

## 7. Expected Manuscript Figures And Tables

### Main Candidate Figure

For representative 2D cases, show two aligned curves over a graph-density or
native parameter axis:

1. oracle error curve:

```text
Err_oracle(G_theta)
```

2. non-oracle score curve:

```text
GCV_score(theta), JS_degree(theta), or Edit(theta)
```

Mark:

- `theta_oracle`;
- `theta_hat`;
- near-oracle band, e.g. <= 1.05 x oracle minimum;
- connectivity/MST repair ineligible region if applicable.

Recommended panels:

- paraboloid, `(c1,c2) = (1,1)`, `n = 200`;
- saddle, `(c1,c2) = (1,4)`, `n = 200`;
- one adaptive-radius panel;
- one cKNN panel.

If only one selector succeeds clearly, use it in the main figure and move the
others to supplement.

### Main Summary Table

Rows:

```text
family x selector
```

Columns:

- median oracle rank;
- median oracle error ratio;
- percent within 5% of oracle;
- percent within 10% of oracle;
- median edge-count ratio;
- failure rate;
- boundary-selection rate.

Aggregate over all primary 2D full-run datasets.

### Supplement Tables

Supplement tables should include:

- per-surface and per-curvature selector recovery;
- per-sample-size recovery;
- seed-level selected parameters;
- sentinel fixed-radius/mKNN results if run;
- 3D staged follow-up results if approved and run;
- sensitivity to top 20 versus top 30 variable features;
- sensitivity to GCV noise level and repair-penalty threshold.

### Supplement Figures

Recommended supplement figures:

- oracle error curves for all representative curvature cases;
- all selector score curves by graph family;
- scatterplot of selector score rank versus oracle error rank;
- seed stability violin/strip plots for selected parameter rank;
- edge-count versus oracle error with selected and oracle points highlighted;
- failure-mode examples where selectors disagree.

## 8. Smoke And Full Run Design

### Minimal Smoke Run

Purpose: test the full data flow, not establish conclusions.

Datasets:

```text
surfaces: paraboloid, saddle
curvature: (1,1)
n: 100
seeds: 1, 2
domain: parameter disk
```

Graph families:

```text
adaptive-radius
cKNN
sKNN
iKNN
```

Parameter grid:

```text
adaptive-radius:
  k_scale in {4, 8}
  radius_rule in {"max"}
  radius_factor in {1.0, 1.5}

cKNN:
  k_scale in {4, 8}
  delta in {1.0, 1.5}

sKNN/iKNN:
  k in {4, 8, 12}
```

Selectors:

```text
GCV on ambient + smooth polynomial features
JS degree stabilization
edge Jaccard/edit stabilization
MST repair burden
```

Outputs expected:

```text
results/datasets.csv
results/settings.csv
results/oracle_metrics.csv
results/non_oracle_scores.csv
results/selector_choices.csv
results/selector_evaluation.csv
figures/representative_selector_curves/*.png
run_config.json
progress.log
```

No dashboard should be built for smoke. A compact HTML or Markdown report can
be proposed only after smoke outputs are validated.

### Full 2D Run

Run only after smoke verifies joins, selectors, and failure flags.

Datasets:

```text
surfaces: paraboloid, saddle
curvatures: (1,1), (1,2), (1,4), (2,4)
n: 100, 200, 400
seeds: 1, 2, 3, 4, 5
domain: parameter disk
```

Families:

```text
primary: adaptive-radius, cKNN
baselines: sKNN, iKNN
sentinels: fixed-radius, mKNN only if approved and cheap
```

Outputs should live under a new manuscript-local or geodesic-data-geometry run
directory, not in `gflow` source code. Recommended location:

```text
/Users/pgajer/current_projects/geodesic_data_geometry/experiments/non-oracle-graph-selection/
```

or, if the run is SIMODS-specific:

```text
/Users/pgajer/current_projects/grip/papers/data_to_geodesic_embedding_paper/experiments/non-oracle-graph-selection/
```

The final decision on run home should be made before implementation. Reusable
functions should be proposed for `gflow`, but this handoff does not authorize
package edits.

### 3D Staged Follow-Up

Do not include 3D in the first full run. After the 2D selector benchmark is
stable, run:

```text
domain: ball
index classes: k = 0-equivalent, k = 1
coefficients: (1,1,1), (1,2,4)
n: 80, 120
seeds: 2
```

Use `quadform.delaunay.geodesic.distances()` with the established reference
settings:

```text
edge.length.factor = 4
delaunay.backend = "cpp" if parity remains acceptable
```

Expand to `(1,4,4)`, `n = 200`, and cube domains only if the staged smoke is
computationally and scientifically clean.

## 9. Existing gflow Assets

### Reusable Graph Constructors

The following are directly reusable:

- `create.adaptive.radius.graph()` in `gflow/R/radius_graphs.R`;
- `create.cknn.graph()` in `gflow/R/radius_graphs.R`;
- `create.sknn.graph()` in `gflow/R/sknn_graphs.R`;
- `create.single.iknn.graph()` and `create.iknn.graphs()` in
  `gflow/R/iknn_graphs.R`;
- `create.mknn.graph()` in `gflow/R/mknn_graphs.R`;
- `create.radius.graph()` in `gflow/R/radius_graphs.R`.

These constructors already carry component/MST repair diagnostics such as
`n_components_before`, `n_components_after`, `mst_edge_matrix`, and
`n_mst_edges_added`.

### Reusable Geodesic And Oracle Utilities

2D utilities in `gflow/R/quadform_geodesics.R`:

- `quadform.embed()`;
- `quadform.gradient()`;
- `quadform.metric()`;
- `quadform.edge.length()`;
- `quadform.edge.lengths()`;
- `quadform.reference.geodesics()`;
- `quadform.grid.geodesic.distances()`;
- `quadform.grid.geodesic.calibration()`;
- `quadform.sample.dataset()`.

3D utilities:

- internal epsilon-net/reference helpers in `quadform_geodesics.R`;
- `quadform.delaunay.geodesic.distances()`.

Graph-distance utility:

- `graph.geodesic.distances()` in `gflow/R/graph_geodesic_distances.R`.

Oracle metric utilities:

- `summarize.isometry.deviation()`;
- `isometry.scale()`;
- `isometry.rel.rms.error()`;
- `isometry.geodesic.diagnostics()`.

### Reusable Non-Oracle Diagnostics

Graph-summary stability:

- `compute.graph.summary.pmf()`;
- `graph.summary.divergence()`;
- `compute.graph.summary.stability()`;
- `jensen.shannon.divergence()`.

iKNN stability:

- `compute.stability.metrics()`;
- `compute.edit.distances()`;
- `internal.compute.edit.distances()`;
- `find.optimal.k()`;
- `build.iknn.graphs.and.selectk()`.

Graph edit:

- `graph.edit.distance()`;
- `calculate.edit.distances()`.

Smoothing/GCV:

- `data.smoother()`;
- `fit.rdgraph.regression()`;
- `refit.rdgraph.regression()`;
- internal `select.eta.gcv.single()` and `select.eta.gcv.single.fast()`;
- `k.diagnostics.plots()` and `get.rcx.optimal.k()` for historical k-selection
  diagnostics.

### Reusable Benchmark Scripts / Reports

2D:

- `gflow/dev/data-geodesic-reconstruction/quadform-curvature-tier1-benchmark/run_quadform_curvature_tier1_benchmark.R`;
- `gflow/dev/data-geodesic-reconstruction/quadform-tier2-domain-sampling-benchmark/run_quadform_tier2_domain_sampling_benchmark.R`;
- `gflow/dev/data-geodesic-reconstruction/quadform-surface-unpruned-benchmark/`.

3D:

- `gflow/dev/data-geodesic-reconstruction/quadform-3d-smoke-benchmark/run_quadform_3d_smoke_benchmark.R`;
- `gflow/dev/data-geodesic-reconstruction/quadform-3d-delaunay-reference/run_quadform_3d_delaunay_reference_report.R`;
- `gflow/dev/data-geodesic-reconstruction/quadform-3d-delaunay-oracle-stress-test/run_quadform_3d_delaunay_oracle_stress_test.R`.

Tests:

- `gflow/tests/testthat/test-quadform-geodesics.R`;
- graph-constructor tests for sKNN, mKNN, radius, component-MST repair, graph
  summary divergence, and graph geodesic distances.

### Missing Implementation Pieces

The main missing pieces are benchmark glue, not new scientific machinery:

1. A graph-family-general parameter grid builder for adaptive-radius, cKNN,
   sKNN, iKNN, and optional sentinels.
2. A common graph extractor that returns `adj_list`, `weight_list`, edge table,
   raw/repaired stage, and MST burden for every family.
3. A common dense/spectral Laplacian GCV scorer that accepts any graph object
   and a feature matrix.
4. Feature-signal generator for 2D and 3D quadform samples.
5. Ordered parameter-sequence definitions for multi-parameter families.
6. Normalized edge symmetric-difference/Jaccard stability for arbitrary graph
   sequences.
7. Selector choice rules with explicit tie and near-optimum handling.
8. A results schema joining datasets, settings, oracle metrics, non-oracle
   scores, selector choices, and selector-vs-oracle evaluation.
9. Plot builders for oracle-vs-selector curves.
10. A smoke report that validates counts and joins without creating a large
    dashboard.

## 10. Risks And Failure Modes

### Feature Signals Too Easy

Ambient coordinates may make GCV look good simply because every reasonable
geometric graph smooths them. Mitigation: use them for smoke only, then add
noisy smooth functions and distractors.

### GCV Selects Smoothness, Not Geometry

GCV may prefer graphs that oversmooth signals rather than graphs that recover
surface geodesics. Mitigation: compare GCV-selected oracle error ratios, edge
counts, shortcut fraction, and signed bias; include noise/distractor feature
sensitivity.

### Multiple Near-Oracle Optima

The oracle error curve may have a broad flat basin. A selector should not be
penalized harshly for choosing a different parameter inside the near-oracle
set. Mitigation: report 5% and 10% near-oracle success in addition to exact
oracle rank.

### Instability Across Seeds

Selectors may work only for favorable samples. Mitigation: report stability
across at least five seeds in the full 2D run and flag high selected-parameter
IQR.

### Selector Disagreement

GCV, JS stabilization, edit stability, and repair burden may disagree for valid
reasons. Mitigation: treat GCV as primary candidate and the stability criteria
as diagnostics unless empirical recovery shows a combined score is better.

### Graph-Family Parameter-Scale Mismatches

`k`, `k_scale`, `delta`, and `radius_factor` are not directly comparable.
Mitigation: compare selector recovery within each family first; use graph edge
count/density as a secondary common axis.

### Boundary Optima

If oracle or selector optima occur at the edge of the parameter grid, the grid
is too narrow. Mitigation: flag boundary optima and do not treat such cases as
final evidence until the grid is expanded.

### Connectivity Repair Dominates The Choice

If MST repair adds many edges, the selected graph may mostly reflect repair
behavior. Mitigation: report raw components and repair burden; use repair
burden as an eligibility filter or penalty.

### 3D Reference Cost And Boundary Behavior

3D Delaunay reference generation may dominate runtime and boundary cases may
need tolerant interpretation. Mitigation: keep 3D staged and use established
reference settings only after the 2D selector pipeline is stable.

## Go / No-Go Criteria After Smoke

Proceed to the full 2D run only if:

- every planned output table is produced and joinable by dataset and setting;
- all selectors return finite scores for at least adaptive-radius and cKNN;
- oracle metrics reproduce expected fields from `summarize.isometry.deviation()`;
- failure flags correctly identify disconnected/high-repair cases;
- at least one representative curve plot shows both `theta_hat` and
  `theta_oracle`;
- the smoke run reveals no systematic boundary optimum requiring grid redesign.

If smoke fails, revise the selector interface and parameter grids before any
full run.

## Proposed Immediate Next Task

After review and approval, implement only the smoke pipeline first. Do not run
the full 2D grid or any 3D selector benchmark until smoke output schemas and
plots are reviewed.
