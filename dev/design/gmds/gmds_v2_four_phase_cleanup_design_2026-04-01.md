# GMDS v2 Four-Phase Cleanup Design

Date: 2026-04-01

## Purpose

This note translates the revised manuscript formulation of Geodesic MDS (GMDS)
into an implementation roadmap for `grip`.

The key manuscript-level clarifications are now:

- **GMDS proper** is the fixed-path all-pairs graph-geodesic objective on a
  connected weighted graph together with a chosen geodesic family
  `(G, Gamma)`.
- **Regularized GMDS** is a separate object of study:
  `E_geo(Z) + lambda R(Z)`.
- **Edge relaxation** is a surrogate objective, not GMDS itself.

That distinction should now be made explicit in the implementation, public API,
tests, and benchmark reports.

This note proposes a four-phase cleanup to get there.

## Current State

The current implementation already contains most of the raw ingredients:

- graph-native all-pairs geodesic preparation:
  [`grip.prepare.geodesic.kk()`](https://github.com/pgajer/grip/blob/main/R/grip_quality.R#L1820)
- data-native `k`-NN preparation:
  [`grip.prepare.geodesic.mds()`](https://github.com/pgajer/grip/blob/main/R/grip_quality.R#L3363)
- pure fixed-path GMDS scoring:
  [`grip.score.geodesic.mds()`](https://github.com/pgajer/grip/blob/main/R/grip_quality.R#L3472)
- pure fixed-path GMDS optimization:
  [`grip.optimize.geodesic.mds()`](https://github.com/pgajer/grip/blob/main/R/grip_quality.R#L3649)
- flattened compiled path-stress kernel:
  [`geodesic_mds_rcpp.cpp`](https://github.com/pgajer/grip/blob/main/src/geodesic_mds_rcpp.cpp#L834)

The generic graph terms proposed earlier also already exist in prototype form:

- graph-edge springs:
  [`grip.geodesic.mds.edge.spring.stats()`](https://github.com/pgajer/grip/blob/main/R/grip_quality.R#L2788)
- graph-aware repulsion:
  [`grip.geodesic.mds.repulsion.stats()`](https://github.com/pgajer/grip/blob/main/R/grip_quality.R#L2862)
- graph-term cache builder:
  [`grip.geodesic.mds.ensure.graph.term.cache()`](https://github.com/pgajer/grip/blob/main/R/grip_quality.R#L2720)

However, these additions are currently mixed into the same public GMDS API, and
the spring/repulsion terms still force a fallback to the R engine rather than
running in the compiled path.

## Cleanup Goal

The goal is not to throw away the current implementation. It is to make it
match the manuscript:

1. **pure GMDS** becomes the stable, graph-first core implementation,
2. **regularized GMDS** becomes an explicitly named extension of that core,
3. **edge relaxation** becomes an explicitly named surrogate method,
4. the benchmark and test suite compare these three objects honestly.

## Phase 1: Make The Public API Match The Manuscript

### Objective

Expose the graph-first structure of the method directly in the public API and
stop conflating graph-native GMDS with the `k`-NN convenience wrapper.

### Why This Comes First

The revised manuscript now defines GMDS first for an arbitrary connected
weighted graph `(G, Gamma)` and only later specializes to the `k`-NN graph
built from data. The code should reflect that same hierarchy.

### Proposed Changes

#### 1. Add a graph-first prepare entry point

Add a public function with a name such as:

- `grip.prepare.graph.geodesic.mds()`

This function should accept:

- `edges`
- `n`
- `adj_list`
- `weight_list`
- `edge_weights`
- `tie_mode = c("single", "average")`

Internally, it can largely wrap the same lower-level graph preparation used by
[`grip.prepare.geodesic.kk()`](https://github.com/pgajer/grip/blob/main/R/grip_quality.R#L1820),
but it should return an object whose semantics are clearly GMDS rather than KK.

#### 2. Recast the current data-native function as a convenience wrapper

Keep:

- [`grip.prepare.geodesic.mds()`](https://github.com/pgajer/grip/blob/main/R/grip_quality.R#L3363)

but explicitly define it as:

1. build a symmetric weighted `k`-NN graph,
2. choose the geodesic family `Gamma_k`,
3. call the graph-first GMDS prepare function.

That keeps the existing user-facing convenience while aligning the package with
the manuscript.

#### 3. Make `Gamma` explicit in the docs

The implementation does not need a separate explicit `Gamma` object in the API
if `tie_mode` remains the operative choice, but the documentation should state
clearly that:

- GMDS is defined on `(G, Gamma)`,
- `tie_mode = "single"` means a deterministic chosen geodesic family,
- `tie_mode = "average"` means the exact average over all tied shortest-path
  families.

#### 4. Tighten naming around the objective

The current names are fairly good, but the documentation of
[`grip.score.geodesic.mds()`](https://github.com/pgajer/grip/blob/main/R/grip_quality.R#L3472)
and
[`grip.optimize.geodesic.mds()`](https://github.com/pgajer/grip/blob/main/R/grip_quality.R#L3649)
should emphasize:

- these are the **pure fixed-path GMDS** functions by default,
- extra Euclidean regularizers are opt-in extensions,
- edge-only surrogates are not the same thing.

### Deliverables

- graph-first GMDS prepare function
- updated docs for graph-first and data-first preparation
- clear explanation of `(G, Gamma)` in the help text
- compatibility-preserving wrappers so existing tests and reports still run

## Phase 2: Separate Pure GMDS, Regularized GMDS, and Edge Relaxation

### Objective

Make the three objective families explicit in code rather than mixing them under
one name.

### Why This Matters

The manuscript now makes a clean conceptual distinction:

- `E_geo` is GMDS,
- `E_geo + lambda R` is regularized GMDS,
- `E_K` is a relaxation/surrogate.

The implementation should mirror that distinction.

### Proposed Changes

#### 1. Keep pure GMDS as the default meaning of `grip.optimize.geodesic.mds()`

The default path should remain:

- pure path stress,
- optional `tie_mode`,
- optional initialization strategy,
- optional anchor only if explicitly requested.

In particular, edge springs and graph-aware repulsion should not quietly turn
`grip.optimize.geodesic.mds()` into a different method without the user noticing.

#### 2. Create a dedicated regularized-GMDS interface

Add an explicitly named interface such as:

- `grip.optimize.regularized.geodesic.mds()`
- `grip.score.regularized.geodesic.mds()`

with scalarized objective

`E_reg(Z) = E_geo(Z) + lambda_anchor E_anchor(Z) + lambda_rep E_rep(Z) + lambda_smooth E_smooth(Z)`

where the public naming emphasizes that these are regularizers, not core GMDS.

The existing anchor, smoothness, and later generic graph terms can be wired
through this interface rather than overloading the base GMDS name.

#### 3. Create a dedicated edge-relaxation interface

Add an explicitly named surrogate interface such as:

- `grip.optimize.edge.relaxed.geodesic.mds()`
- `grip.score.edge.relaxed.geodesic.mds()`

for the edge-only objective

`E_K(Z) = sum_{(u,v) in E} k_uv (||z_u-z_v|| - w_uv)^2`

This function should be documented as:

- a computationally cheaper graph-drawing-style relaxation,
- useful as a surrogate or warm start,
- not identical to all-pairs GMDS.

#### 4. Keep mesh-specific regularizers clearly out of the core path

The bending and mesh-only code in:

- [`zz_geodesic_mds_bending_core.R`](https://github.com/pgajer/grip/blob/main/R/zz_geodesic_mds_bending_core.R)
- [`zz_geodesic_mds_bending_optimize.R`](https://github.com/pgajer/grip/blob/main/R/zz_geodesic_mds_bending_optimize.R)
- [`geodesic_mds_bending_rcpp.cpp`](https://github.com/pgajer/grip/blob/main/src/geodesic_mds_bending_rcpp.cpp)

should remain available for specialized experiments, but they should not define
what GMDS means in the package.

### Deliverables

- cleanly separated function families
- updated help files and roxygen docs
- clearer naming in reports and tests

## Phase 3: Finish The Generic Graph-Regularized Engine In C++

### Objective

Take the graph-generic spring/repulsion regularizers that already exist in R and
port them into the compiled flat kernel so regularized GMDS is practical on
nontrivial graphs.

### Why This Matters

Right now the graph-generic terms are conceptually the right direction, but they
are operationally limited because requesting them forces an R-engine fallback.
That makes them too slow for broad comparisons and too easy to dismiss as merely
prototype code.

### Current Pieces Already In Place

R already contains:

- graph-edge spring gradient and energy:
  [`grip.geodesic.mds.edge.spring.stats()`](https://github.com/pgajer/grip/blob/main/R/grip_quality.R#L2788)
- graph-aware repulsion gradient and energy:
  [`grip.geodesic.mds.repulsion.stats()`](https://github.com/pgajer/grip/blob/main/R/grip_quality.R#L2862)
- automatic repulsion-pair cache construction:
  [`grip.build.graph.repulsion.cache()`](https://github.com/pgajer/grip/blob/main/R/grip_quality.R#L2660)

The flat compiled GMDS path already supports:

- flattened pair-edge cache
- anchor schedule
- smoothness schedule
- multithreading

in:

- [`grip_optimize_geodesic_mds_flat_cpp()`](https://github.com/pgajer/grip/blob/main/src/geodesic_mds_rcpp.cpp#L1365)

### Proposed C++ Additions

#### 1. Extend the flat cache view

Add flattened arrays for:

- graph edges:
  - `flat_graph_edge_u`
  - `flat_graph_edge_v`
  - `flat_graph_edge_target`
- repulsion pairs:
  - `flat_repulsion_u`
  - `flat_repulsion_v`
  - `flat_repulsion_target`

These should be prepared once on the R side and then passed to the compiled
optimizer just as the flattened path cache already is.

#### 2. Extend `evaluate_flat_state`

Add compiled accumulation functions parallel to the existing path and smoothness
terms:

- `accumulate_flat_edge_springs(...)`
- `accumulate_flat_repulsion(...)`

Then define the total compiled state as:

`E_total = E_path + E_anchor + E_smooth + E_edge + E_rep`

with per-term trace columns.

#### 3. Add schedules for the new weights

Reuse the same scheduling approach already used for:

- `anchor_weights`
- `smooth_weights`

so the flat kernel can support:

- `edge_spring_weights`
- `repulsion_weights`

with constant / linear / geometric continuation.

#### 4. Preserve R as the reference implementation

The R implementation should remain the audit/reference path.
Every compiled addition should be checked against the corresponding R version
with finite-difference and numeric-agreement tests.

### Deliverables

- compiled graph-edge spring term
- compiled graph-aware repulsion term
- no more forced R fallback for generic regularized GMDS
- numeric-agreement tests between R and C++ paths

## Phase 4: Rebuild The Evaluation Suite Around The New Hierarchy

### Objective

Reorganize the experiments so they answer the right questions for the revised
formulation.

This phase should not start by mixing all graph families again. It should focus
first on **parabolic meshes**, because that is where the pathology story is
clearest and where the difference between pure GMDS, regularized GMDS, and the
edge relaxation is easiest to interpret.

### Phase 4A: Validate Pure GMDS On Foundation Cases

Before aesthetic or performance comparisons, validate the mathematical core:

#### Exact realizability cases

Use small graphs with exact zero-stress embeddings, such as:

- path graphs,
- shared-edge triangle examples,
- small trees,
- small cycles with exact geometric realizations.

Check:

- pure GMDS can reach near-zero raw stress,
- different exact realizations have the same objective value,
- the implementation matches the manuscript’s nonuniqueness claims.

#### Tie-handling cases

Use graphs with tied shortest paths, such as:

- squares,
- diamonds,
- orthogonal lattice patches.

Check:

- `tie_mode = "single"` and `tie_mode = "average"` are genuinely different,
- `tie_mode = "average"` restores symmetry on symmetric graphs,
- behavior under relabeling is understood and documented.

This foundation work keeps the implementation anchored to the actual GMDS
definition before any regularization is added.

### Phase 4B: Test Regularized GMDS As A Separate Object Of Study

Once pure GMDS is validated, study regularized GMDS explicitly rather than
treating it as a hidden tweak.

#### Focus first on parabolic meshes

Use orthogonal paraboloid meshes at increasing densities:

- `12 x 12`
- `15 x 15`
- `20 x 20`

and compare:

- classical MDS baseline,
- pure GMDS,
- regularized GMDS with anchor only,
- regularized GMDS with anchor plus graph-aware repulsion,
- optional anchor plus graph-edge springs if helpful.

#### Keep `lambda` explicit

The report and tables should not hide the regularization strength.
For every run, record:

- `lambda_anchor`
- `lambda_repulsion`
- `lambda_edge`
- schedule type
- runtime
- raw GMDS stress
- any regularizer penalties
- Procrustes error to the reference paraboloid
- qualitative notes about singularity / fold formation

#### Main question

The question here is:

`Can explicit regularization repair the paraboloid pathology while still preserving the geodesic objective well enough to be scientifically meaningful?`

### Phase 4C: Test The Edge Relaxation As A Surrogate, Not As “GMDS”

Only after the pure and regularized GMDS experiments should the edge-only
surrogate be benchmarked.

#### Again focus first on parabolic meshes

Use the same paraboloid benchmark family and compare:

- pure GMDS,
- regularized GMDS,
- edge-relaxed surrogate,
- optionally edge-relaxed plus repulsion.

#### Keep the interpretation honest

The edge relaxation is not supposed to “win the GMDS benchmark” by definition.
Instead, study it as:

- a possible warm start,
- a faster approximate surrogate,
- a possible way to avoid some pathologies by controlling only local edge
  distortion.

#### Main question

The question here is:

`Does the cheaper edge-only objective produce better-behaved paraboloid embeddings, and if so, what geodesic fidelity is lost compared with pure and regularized GMDS?`

### Phase 4D: Expand Beyond Paraboloids Only After The Hierarchy Is Stable

Once the paraboloid story is clear, extend the same framework to:

- flat orthogonal meshes,
- ripple surfaces in 3D,
- other simple graph families.

But this should come after the parabolic-mesh comparison framework is stable and
well-instrumented.

### Deliverables

- a reorganized benchmark/report pipeline with separate sections for:
  - pure GMDS,
  - regularized GMDS,
  - edge relaxation
- parabolic-mesh-first reports
- per-method runtime recording
- interactive 3D inspection assets where useful

## Recommended Immediate Order

If we start implementation after review of this plan, I would do the phases in
this order:

1. Phase 1 API cleanup
2. Phase 2 objective-family separation
3. Phase 3 compiled generic regularization
4. Phase 4 evaluation rebuild, starting with paraboloid families

That order keeps the semantics clear before we invest more into performance or
benchmark generation.

## Summary

The main idea is simple:

- keep the existing pure GMDS implementation as the mathematical core,
- stop overloading that name with every later modification,
- make regularized GMDS and edge relaxation explicit alternative objects,
- and evaluate them first on parabolic meshes, where the pathology question is
  most informative.

That gives us an implementation that is finally aligned with the revised
manuscript rather than merely inspired by earlier experimental variants.
