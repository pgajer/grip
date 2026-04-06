# MISF Geodesic KK Design

Date: 2026-04-05

Status: **Implemented baseline public family. Paper-grade benchmarking still pending.**

## Purpose

This note turns the proposed `grip.*misf*.geodesic.kk()` family into a concrete
design for `grip`.

The motivating package-level conclusion is:

- `grip` already has a strong software story around fast multiscale 2D and 3D
  layout,
- the package already includes exact geodesic-KK scoring and optimization,
- the package already includes standalone landmark geodesic-KK tooling,
- the package already includes MISF-based multiscale GMDS,
- but it does **not** yet expose a standalone MISF-based geodesic-KK family.

The next method-family addition was therefore a multiscale
**MISF-geodesic-KK** pipeline that reuses GRIP's maximal independent set
filtration as a scaffold for scalable weighted geodesic-KK optimization.

## Scope

The work covered by this design is:

1. add a new public MISF-GKK family with prepare, optimize, and score entry
   points,
2. reuse exact full GKK on small or coarse active sets,
3. reuse LGKK as the sparse approximation path once active sets become too
   large for practical all-pairs refinement,
4. keep the API close to the current `grip.prepare/optimize/score.geodesic.kk()`
   and `grip.optimize.misf.geodesic.mds()` families.

This design does **not** aim to:

- replace `grip.optimize.geodesic.kk()` for small graphs,
- merge MISF-GKK into `grip.layout*()` before the new family is validated,
- add a new graph-first `grip.prepare.graph.geodesic.kk()` layer, since
  `grip.prepare.geodesic.kk()` already serves that role for weighted graphs,
- introduce a compiled exact GKK optimizer immediately.

## Relevant Current Code

- exact full geodesic-KK preparation, scoring, and optimization in
  [grip_quality.R](/Users/pgajer/current_projects/grip/R/grip_quality.R#L1790)
- standalone landmark geodesic-KK preparation, scoring, and optimization in
  [grip_quality.R](/Users/pgajer/current_projects/grip/R/grip_quality.R#L1713)
- MISF-based GMDS preparation, optimization, and scoring in
  [grip_geodesic_misf.R](/Users/pgajer/current_projects/grip/R/grip_geodesic_misf.R#L1658)
- MISF extraction API in
  [grip_misf.R](/Users/pgajer/current_projects/grip/R/grip_misf.R)
- multiscale LGKK integration in the main GRIP solver in
  [grip_layout.R](/Users/pgajer/current_projects/grip/R/grip_layout.R#L849)

## Core Design Decision

The MISF-GKK family should **reuse MISF as a multiscale scaffold**, but it
should **not** be implemented as more conditionals inside the existing
`grip.layout*()` force-directed solver.

Instead, the recommended architecture is:

1. keep `grip.prepare.geodesic.kk()` as the full-graph weighted geodesic cache,
2. layer MISF metadata and per-level active-set metadata on top of that cache,
3. solve the coarsest level with exact GKK when affordable,
4. insert finer levels with a KK-aware warm start,
5. refine each active level under either exact GKK or sparse LGKK according to
   an explicit pair policy,
6. finish with a final full-graph geodesic-KK polish under exact or sparse
   scoring depending on graph size.

This keeps the new family mathematically tied to geodesic-KK rather than
turning it into another force-directed GRIP variant.

## Proposed Public Surface

The recommended exported family is:

- `grip.prepare.misf.geodesic.kk()`
- `grip.optimize.misf.geodesic.kk()`
- `grip.score.misf.geodesic.kk()`

Recommended S3 classes:

- `grip_misf_gkk_prepared`
- `grip_misf_gkk_fit`

Recommended class layering:

- prepared objects should inherit from
  `c("grip_misf_gkk_prepared", "grip_gkk_prepared", "grip_geodesic_kk_prepared", "list")`
- fit objects should inherit from
  `c("grip_misf_gkk_fit", "list")`

This allows the final full-graph object to remain compatible with existing
exact GKK validators and scorers where possible.

## Algorithm To Implement

Let

`V = V_0 ⊃ V_1 ⊃ ... ⊃ V_k`

be the MIS filtration of the weighted graph.

The proposed MISF-GKK algorithm is:

1. prepare the full graph with `grip.prepare.geodesic.kk()`;
2. extract MISF on the same weighted graph;
3. build the induced coarse weighted graph on the top MISF level;
4. solve the top level under exact GKK if its active-set size is below a
   configurable exact limit, otherwise solve it under LGKK;
5. for each level from coarse to fine:
   - insert newly activated vertices with a KK-aware warm start,
   - build the active-level weighted graph,
   - refine the active level under exact GKK or LGKK according to the level
     pair policy;
6. after reaching the full graph, run a short final geodesic-KK polish under
   exact GKK or LGKK according to the final pair policy;
7. score the result with exact full-graph GKK by default.

The main approximation policy should be explicit and user-visible:

- `pair_mode = "full"` means exact all-pairs GKK on the active set,
- `pair_mode = "landmark"` means LGKK on the active set,
- `pair_mode = "auto"` means use exact GKK when the active set is at most a
  configured size and switch to LGKK above that size.

## Proposed Exported Functions

## 1. `grip.prepare.misf.geodesic.kk()`

### Proposed signature

```r
grip.prepare.misf.geodesic.kk(
  edges = NULL,
  n = NULL,
  adj_list = NULL,
  weight_list = NULL,
  edge_weights = NULL,
  tie_mode = c("single", "average"),
  num_init = 24L,
  num_nbrs = 20L,
  dim = 2L,
  top_level_mode = c("solve", "skip"),
  top_level_pair_mode = c("auto", "full", "landmark"),
  top_level_full_limit = 512L,
  top_level_local_nbrs = 20L,
  top_level_landmark_count = 8L,
  top_level_restarts = 8L,
  top_level_max_iter = 16L,
  top_level_init = c("cmdscale", "random"),
  seed = 6L
)
```

### Argument rationale

Graph inputs:

- `edges`, `n`, `adj_list`, `weight_list`, `edge_weights`
- exactly mirror the existing GKK and MISF-GMDS graph-first entry points

Base geodesic behavior:

- `tie_mode`
- should match `grip.prepare.geodesic.kk()`

MISF scaffold:

- `num_init`
- `num_nbrs`
- these should match the current MISF extraction controls

Top-level multiscale solve:

- `dim`
- `top_level_mode`
- `top_level_pair_mode`
- `top_level_full_limit`
- `top_level_local_nbrs`
- `top_level_landmark_count`
- `top_level_restarts`
- `top_level_max_iter`
- `top_level_init`
- `seed`

### Recommended defaults

- `tie_mode = "average"` is likely the stronger default for weighted geometry
  tasks, but the public function should still expose both values
- `top_level_pair_mode = "auto"`
- `top_level_full_limit = 512L`
- `top_level_local_nbrs = 20L`
- `top_level_landmark_count = 8L`
- `top_level_init = "cmdscale"`
- `top_level_restarts = 8L`
- `top_level_max_iter = 16L`

### Proposed return object

Return class:

- `grip_misf_gkk_prepared`

Return fields:

Full-graph exact GKK cache inherited from `grip_gkk_prepared`:

- `n`
- `edges`
- `edge_targets`
- `adj_list`
- `weight_list`
- `pair_matrix`
- `pair_graph_distance`
- `path_vertices`
- `path_edges`
- `path_edge_weights`
- `pair_path_count_log`
- `flat_pair_edge_offsets`
- `flat_edge_u`
- `flat_edge_v`
- `flat_edge_coeff`
- `graph_diameter`
- `distance_matrix`
- `pair_mode`
- `tie_mode`

MISF scaffold metadata:

- `misf`
- `level_vertices`
- `active_levels`
- `insertion_order`
- `top_level_index`
- `top_level_level`
- `top_level_vertices`

Top-level graph and prepared objects:

- `top_level_graph`
- `top_level_pair_mode`
- `top_level_effective_pair_mode`
- `top_level_full_limit`
- `top_level_local_nbrs`
- `top_level_landmark_count`
- `top_level_prepared_full`
- `top_level_prepared_sparse`
- `top_level_prepared`
- `top_level_dim`
- `top_level_mode`
- `top_level_restarts`
- `top_level_max_iter`
- `top_level_init`
- `top_level_fit`

Stored defaults for later optimization:

- `multiscale_mode`
- `insertion_anchor_policy`
- `insertion_anchor_count`
- `insertion_anchor_weight_mode`
- `insertion_max_iter`
- `insertion_mode`
- `insertion_layout_k`
- `insertion_weighted_preset`
- `insertion_grip_args`
- `insertion_weighted_args`
- `insertion_fr_niter`
- `refinement_pair_mode`
- `refinement_full_limit`
- `refinement_local_nbrs`
- `refinement_landmark_count`
- `refinement_max_iter`
- `final_pair_mode`
- `final_full_limit`
- `final_local_nbrs`
- `final_landmark_count`
- `final_max_iter`
- `stiffness`
- `distance_floor`
- `edge_length_epsilon`
- `initial_step`
- `step_shrink`
- `armijo_factor`
- `grad_tol`
- `min_step`
- `recenter`
- `misf_seed`

### Notes

- `top_level_prepared` should point to whichever cache is actually active at the
  top level.
- `top_level_prepared_full` and `top_level_prepared_sparse` should be optional
  fields so later scoring and diagnostics can compare exact and approximate
  stage behavior without rebuilding from scratch.

## 2. `grip.optimize.misf.geodesic.kk()`

### Proposed signature

```r
grip.optimize.misf.geodesic.kk(
  prepared = NULL,
  edges = NULL,
  n = NULL,
  adj_list = NULL,
  weight_list = NULL,
  edge_weights = NULL,
  tie_mode = NULL,
  num_init = 24L,
  num_nbrs = 20L,
  dim = NULL,
  top_level_pair_mode = NULL,
  top_level_full_limit = NULL,
  top_level_local_nbrs = NULL,
  top_level_landmark_count = NULL,
  top_level_restarts = NULL,
  top_level_max_iter = NULL,
  top_level_init = NULL,
  insertion_anchor_policy = NULL,
  insertion_anchor_count = NULL,
  insertion_anchor_weight_mode = NULL,
  insertion_max_iter = NULL,
  insertion_mode = NULL,
  insertion_layout_k = NULL,
  insertion_weighted_preset = NULL,
  insertion_grip_args = NULL,
  insertion_weighted_args = NULL,
  insertion_fr_niter = NULL,
  refinement_pair_mode = NULL,
  refinement_full_limit = NULL,
  refinement_local_nbrs = NULL,
  refinement_landmark_count = NULL,
  refinement_max_iter = NULL,
  final_pair_mode = NULL,
  final_full_limit = NULL,
  final_local_nbrs = NULL,
  final_landmark_count = NULL,
  final_max_iter = NULL,
  stiffness = 1.0,
  distance_floor = 1e-8,
  edge_length_epsilon = 1e-8,
  initial_step = 1.0,
  step_shrink = 0.5,
  armijo_factor = 1e-4,
  grad_tol = 1e-8,
  min_step = 1e-8,
  recenter = TRUE,
  return_trace = FALSE,
  return_frames = FALSE,
  seed = 6L
)
```

### Argument rationale

Preparation reuse:

- `prepared`
- raw graph inputs when `prepared` is omitted
- `tie_mode`
- `num_init`
- `num_nbrs`
- `dim`

These should mirror the current MISF-GMDS family so users can build or reuse a
prepared object in the same way.

Top-level stage:

- `top_level_pair_mode`
- `top_level_full_limit`
- `top_level_local_nbrs`
- `top_level_landmark_count`
- `top_level_restarts`
- `top_level_max_iter`
- `top_level_init`

Insertion stage:

- `insertion_anchor_policy`
- `insertion_anchor_count`
- `insertion_anchor_weight_mode`
- `insertion_max_iter`
- `insertion_mode`
- `insertion_layout_k`
- `insertion_weighted_preset`
- `insertion_grip_args`
- `insertion_weighted_args`
- `insertion_fr_niter`

Recommended insertion modes:

- `"geodesic"` for anchor-based weighted trilateration
- `"weighted_kk"` for local weighted-KK warm starts on the active level
- `"weighted_grip"` for weighted GRIP warm starts
- `"grip"` for combinatorial GRIP warm starts
- `"fr"` for a simple force-directed fallback

Default recommendation:

- `insertion_mode = "weighted_kk"`

Active-level refinement:

- `refinement_pair_mode`
- `refinement_full_limit`
- `refinement_local_nbrs`
- `refinement_landmark_count`
- `refinement_max_iter`

Final polish:

- `final_pair_mode`
- `final_full_limit`
- `final_local_nbrs`
- `final_landmark_count`
- `final_max_iter`

Shared GKK/LGKK optimizer controls:

- `stiffness`
- `distance_floor`
- `edge_length_epsilon`
- `initial_step`
- `step_shrink`
- `armijo_factor`
- `grad_tol`
- `min_step`
- `recenter`

Diagnostics:

- `return_trace`
- `return_frames`
- `seed`

### Recommended defaults

Top level:

- `top_level_pair_mode = "auto"`
- `top_level_full_limit = 512L`
- `top_level_local_nbrs = 20L`
- `top_level_landmark_count = 8L`
- `top_level_restarts = 8L`
- `top_level_max_iter = 16L`
- `top_level_init = "cmdscale"`

Insertion:

- `insertion_anchor_policy = "prev_level_spread"`
- `insertion_anchor_count = 4L` in 2D, `6L` in 3D
- `insertion_anchor_weight_mode = "inverse_graph_distance_sq"`
- `insertion_mode = "weighted_kk"`
- `insertion_layout_k = 6L`
- `insertion_max_iter = 32L`

Refinement:

- `refinement_pair_mode = "auto"`
- `refinement_full_limit = 256L`
- `refinement_local_nbrs = 8L`
- `refinement_landmark_count = 4L`
- `refinement_max_iter = 8L`

Final polish:

- `final_pair_mode = "auto"`
- `final_full_limit = 1024L`
- `final_local_nbrs = 20L`
- `final_landmark_count = 8L`
- `final_max_iter = 8L`

Shared optimizer controls:

- keep the same defaults as `grip.optimize.geodesic.kk()`

### Proposed return object

Return class:

- `grip_misf_gkk_fit`

Return fields:

High-level outputs:

- `coords`
- `prepared`
- `score`

Top-level stage:

- `top_level_fit`
- `top_level_coords_full`
- `top_level_pair_mode`
- `top_level_effective_pair_mode`

Insertion stage:

- `insertion`

Recommended insertion subfields:

- `coords`
- `level_results`
- `level_trace`
- `vertex_trace`

Refinement stage:

- `refinement`

Recommended refinement subfields:

- `coords`
- `level_results`
- `level_trace`

Final polish stage:

- `final_polish`
- `final_pair_mode`
- `final_effective_pair_mode`

Shared diagnostics:

- `stage_trace`
- `trace`
- `frames`
- `timing`
- `resolved_controls`

### Proposed `stage_trace` fields

The stage-level summary table should mirror the GMDS MISF family where
possible, but report GKK-specific metrics:

- `stage`
- `level`
- `active_n`
- `inserted_n`
- `pair_n`
- `pair.mode`
- `energy`
- `weighted.rmse`
- `weighted.rel.rmse`
- `mean.rel.path.error`
- `max_grad_norm`
- `all_converged`
- `elapsed_sec`
- `trace_rows`
- `frame_count`

### Proposed `trace` structure

When `return_trace = TRUE`, return a named list:

- `top_level_trace`
- `top_restart_summary`
- `insertion_level_trace`
- `insertion_vertex_trace`
- `refinement_level_trace`
- `final_polish_trace`

### Proposed `frames` structure

When `return_frames = TRUE`, return a named list:

- `top_level`
- `after_top_level`
- `insertion_levels`
- `after_insertion`
- `refinement_levels`
- `after_refinement`
- `final_polish`
- `final`

## 3. `grip.score.misf.geodesic.kk()`

### Proposed signature

```r
grip.score.misf.geodesic.kk(
  fit = NULL,
  coords = NULL,
  prepared = NULL,
  stiffness = 1.0,
  distance_floor = 1e-8,
  edge_length_epsilon = 1e-8,
  score_pair_mode = c("full", "landmark", "auto"),
  score_full_limit = 2048L,
  score_local_nbrs = 20L,
  score_landmark_count = 8L,
  return_trace = FALSE
)
```

### Argument rationale

Two entry modes:

- `fit` for scoring a full MISF-GKK fit object
- `coords + prepared` for rescoring external coordinates against the same MISF
  prepared object

Shared geodesic-KK controls:

- `stiffness`
- `distance_floor`
- `edge_length_epsilon`

Scoring pair policy:

- `score_pair_mode`
- `score_full_limit`
- `score_local_nbrs`
- `score_landmark_count`

Default recommendation:

- when `fit` is supplied, score the final coordinates with exact full GKK if
  the graph is at most `score_full_limit`, otherwise score with LGKK unless the
  user explicitly forces exact full GKK

### Proposed return object

Return type:

- one-row `data.frame`

Proposed scalar fields:

- `multiscale.mode`
- `n`
- `dim`
- `top.level`
- `top.level.n`
- `top.level.pair.mode`
- `top.level.energy`
- `top.level.weighted.rmse`
- `top.level.weighted.rel.rmse`
- `insertion.level.count`
- `inserted.vertex.count`
- `refinement.level.count`
- `refinement.pair.mode`
- `final.polish.pair.mode`
- `final.polish.trace.rows`
- `elapsed.top.level`
- `elapsed.insertion`
- `elapsed.refinement`
- `elapsed.final.polish`
- `elapsed.total`
- `final.n.pairs`
- `final.pair.mode`
- `final.scale.mode`
- `final.scale.L0`
- `final.scale.L`
- `final.gkk.energy`
- `final.gkk.weighted.rmse`
- `final.gkk.weighted.rel.rmse`
- `final.gkk.mean.abs.path.error`
- `final.gkk.mean.rel.path.error`

Optional list columns when `return_trace = TRUE`:

- `stage.trace`
- `top.restart.summary`
- `insertion.level.trace`
- `insertion.vertex.trace`
- `refinement.level.trace`
- `final.polish.trace`

### Notes

- `grip.score.misf.geodesic.kk()` should summarize the multiscale pipeline in
  the same way `grip.score.misf.geodesic.mds()` summarizes MISF-GMDS.
- The final score should use the exact full-graph GKK scorer whenever practical,
  even if some intermediate levels used LGKK.

## Suggested Internal Helpers

The family will likely need these non-exported helpers:

- `grip.validate.misf.geodesic.kk.prepared()`
- `grip.validate.misf.geodesic.kk.fit()`
- `grip.resolve.misf.geodesic.kk.prepared()`
- `grip.geodesic.misf.resolve.pair.policy()`
- `grip.geodesic.misf.induced_level_graph()`
- `grip.geodesic.misf.prepare.level.gkk()`
- `grip.geodesic.misf.solve.top.level.gkk()`
- `grip.geodesic.misf.insert.level.gkk()`
- `grip.geodesic.misf.refine.level.gkk()`
- `grip.geodesic.misf.final.polish.gkk()`
- `grip.geodesic.misf.build.stage.trace.gkk()`

## Phased Implementation Plan

## Phase 1. Prepared-object and policy layer

### Goal

Create `grip.prepare.misf.geodesic.kk()` and the supporting object model.

### Tasks

- reuse `grip.prepare.geodesic.kk()` as the base full-graph cache
- reuse `grip.build.misf()` for filtration extraction
- add exact-vs-landmark pair-policy resolution helpers
- build the top-level induced weighted graph
- optionally build both top-level exact and sparse prepared caches

### Verification

- prepared objects inherit from `grip_gkk_prepared`
- top-level graph and vertex ids are deterministic
- exact and sparse top-level cache selection is deterministic

## Phase 2. Top-level MISF-GKK solve

### Goal

Solve the coarsest level under geodesic-KK with restart support.

### Tasks

- add top-level initializers:
  - `cmdscale`
  - random
- reuse `grip.optimize.geodesic.kk()` for exact top-level solves
- reuse `grip.optimize.landmark.geodesic.kk()` for sparse top-level solves
- record restart summaries and the chosen best top-level fit

### Verification

- finite coordinates on small weighted graphs
- top-level energy improves over raw initializer energies
- exact and sparse top-level paths agree on tiny graphs where both are feasible

## Phase 3. KK-aware insertion stage

### Goal

Insert finer-level vertices with a warm start that respects weighted graph
geometry better than raw Euclidean GRIP insertion.

### Tasks

- reuse the current MISF anchor-selection patterns
- add `weighted_kk` insertion warm starts on active-level subgraphs
- keep `geodesic` anchor-based trilateration as a simpler fallback
- expose the same optional `grip`, `weighted_grip`, and `fr` warm-start modes
  already used by MISF-GMDS where helpful

### Verification

- inserted vertices are finite immediately after insertion
- KK-based insertion lowers post-insertion active-level GKK energy relative to
  a naive barycentric fallback on weighted mesh tests

## Phase 4. Active-level refinement

### Goal

Run geodesic-KK refinement after each insertion level.

### Tasks

- implement levelwise exact GKK refinement for small active sets
- implement levelwise LGKK refinement for larger active sets
- keep the exact/LGKK choice explicit in level traces

### Verification

- active-level energy decreases on each refinement stage
- `pair_mode = "full"` matches direct exact GKK optimization on tiny levels
- `pair_mode = "landmark"` produces deterministic sparse refinements

## Phase 5. Final full-graph polish and scoring

### Goal

Finish the multiscale pipeline and make its output inspectable.

### Tasks

- add final exact or sparse polish
- add `grip.score.misf.geodesic.kk()`
- align stage summaries with the MISF-GMDS scorer
- attach optional traces and frames

### Verification

- fit objects produce stable score summaries
- exact final scoring agrees with standalone `grip.score.geodesic.kk()` on the
  final coordinates
- stage trace tables are consistent with the stored timings and frames

## Testing Plan

### Structural tests

- prepared object class and field presence
- deterministic MISF level extraction on weighted graphs
- deterministic pair-policy resolution

### Exactness tests

- weighted path graphs
- small weighted rectangles
- tiny trees with analytically simple geodesic targets

### Multiscale behavior tests

- top-level energy lower than initializer energy
- levelwise refinement does not increase stage energy on tiny deterministic
  cases
- final score lower than score from the top-level-only embedding on small
  weighted meshes

### Approximation boundary tests

- `pair_mode = "auto"` matches `full` below the configured limit
- `pair_mode = "auto"` matches `landmark` above the configured limit

## Open Questions

1. Should `top_level_init = "cmdscale"` remain the default if weighted meshes
   show strong basin sensitivity even after multiscale insertion?
2. Should `score_pair_mode = "full"` be the scoring default regardless of graph
   size, or should the scorer default to `auto` for consistency with the
   optimizer?
3. Should later versions expose stage-specific optimizer controls such as
   `scale_mode`, or should the first public version keep those fixed for
   simplicity?
4. Should the family eventually add a compiled exact GKK kernel, or is the
   intended large-graph path always going to be MISF plus LGKK approximation?

## Recommendation

The implementation should start with the three exported functions listed above
and no additional public entry points. The first public milestone should be a
usable, well-scored, traceable MISF-GKK pipeline whose exact-vs-sparse
approximation policy is explicit in the API and in the returned diagnostics.
