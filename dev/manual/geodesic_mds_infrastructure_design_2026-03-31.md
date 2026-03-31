# Geodesic MDS in `grip`: Detailed Infrastructure Design

Date: 2026-03-31

## Goal

Implement a first-class geodesic MDS workflow in `grip` that reuses the package's existing geodesic path infrastructure wherever possible, while keeping the implementation additive and low-risk.

The target behavior is:

1. Build a weighted symmetric `k`-NN graph from an input data matrix.
2. Choose one deterministic shortest path for every unordered vertex pair.
3. Optimize an embedding by minimizing the geodesic MDS numerator
   `sum_{i < j} (h_ij(Z) - delta_ij)^2`,
   where `delta_ij` is the input-side graph geodesic and `h_ij(Z)` is the embedded length of the stored path.
4. Report Stress-1 style diagnostics using the fixed denominator `sum_{i < j} delta_ij^2`.

This matches the manuscript formulation in
[geodesic_mds.tex](/Users/pgajer/current_projects/geodesic_MDS/manuscript/geodesic_mds.tex#L614).

## Existing Assets in `grip`

The package already contains most of the reusable mathematical machinery.

### R path-cache and scoring infrastructure

The recent geodesic-KK code in
[grip_quality.R](/Users/pgajer/current_projects/grip/R/grip_quality.R#L804)
already provides:

- deterministic shortest-path-tree construction for weighted and unweighted graphs,
- deterministic path reconstruction,
- full all-pairs path caches,
- embedded path-length evaluation,
- gradient accumulation by back-propagating residuals along stored path edges,
- deterministic gradient-descent polishing with Armijo backtracking.

The critical reusable pieces are:

- `grip.prepare.geodesic.kk.base()`
- `grip.build.geodesic.kk.path.cache()`
- `grip.geodesic.kk.path.lengths()`
- `grip.geodesic.kk.energy.gradient()`
- `grip.optimize.geodesic.kk()`

### Compiled geodesic layout infrastructure

The multiscale LGKK implementation inside `DrawGraph` already proves that the core geodesic kernel is practical in compiled form.

Relevant code:

- deterministic active-set shortest paths:
  [MishSupport.cpp](/Users/pgajer/current_projects/grip/src/MishSupport.cpp#L201)
- sparse path-cache construction:
  [MishSupport.cpp](/Users/pgajer/current_projects/grip/src/MishSupport.cpp#L399)
- compiled geodesic gradient-descent refinement:
  [MishSupport.cpp](/Users/pgajer/current_projects/grip/src/MishSupport.cpp#L504)

The package also already has the low-level graph and point abstractions that a standalone compiled geodesic-MDS optimizer can reuse:

- `Graph`:
  [Graph.h](/Users/pgajer/current_projects/grip/src/Graph.h)
- `Point`:
  [Point.h](/Users/pgajer/current_projects/grip/src/Point.h)

## Architectural Decision

Geodesic MDS should be implemented as a new sibling solver, not as another branch inside `DrawGraph::mish_engine()`.

### Why not put GMDS inside `DrawGraph`?

`DrawGraph` is a multiscale graph-drawing engine with GRIP-specific concerns:

- MISF scheduling,
- anchor-based insertion,
- FR / KK force phases,
- 2D/3D drawing heuristics,
- final-stage layout polish hooks.

Those are valuable for graph drawing, but geodesic MDS is a different optimization problem. Binding GMDS too tightly to `DrawGraph` would make the implementation harder to reason about and harder to test independently.

### Chosen structure

The first implementation should be split as follows:

#### R layer

The R layer will handle:

- input validation,
- pairwise Euclidean distances on the original data,
- symmetric `k`-NN graph construction,
- optional graph-connection correction via MST augmentation,
- classical-MDS initialization via `cmdscale`,
- score formatting and public API ergonomics.

#### Compiled layer

The compiled layer will handle:

- deterministic all-pairs shortest-path trees on the weighted graph,
- all-pairs chosen-path cache construction,
- repeated embedded path-length evaluation,
- geodesic-MDS energy and gradient accumulation,
- Armijo backtracking gradient descent,
- optional trace / frame capture.

This keeps the data-geometry logic close to the R API and the expensive repeated optimization in C++.

## Public API

The first implementation should add these user-facing functions:

- `grip.prepare.geodesic.mds()`
- `grip.score.geodesic.mds()`
- `grip.optimize.geodesic.mds()`

### `grip.prepare.geodesic.mds()`

Purpose:

- accept an input data matrix `X`,
- build a deterministic symmetric `k`-NN graph with Euclidean edge weights,
- optionally augment it to connectedness using the Euclidean MST,
- reuse the full all-pairs geodesic-KK preparation machinery to obtain a reusable fixed-path cache.

Output:

- a prepared object that contains:
  - the input graph,
  - the weighted adjacency data,
  - all-pairs graph geodesics,
  - the chosen shortest-path family,
  - the graph distance matrix,
  - metadata such as `k` and the connectivity policy.

### `grip.score.geodesic.mds()`

Purpose:

- evaluate an embedding against the fixed-path geodesic-MDS objective.

Primary summary fields:

- raw stress numerator,
- half-energy used by the optimizer,
- normalized geodesic stress,
- RMSE,
- mean absolute path error,
- mean relative path error.

### `grip.optimize.geodesic.mds()`

Purpose:

- optimize a 2D or 3D embedding using the fixed-path geodesic-MDS objective.

Capabilities:

- initialize from classical MDS on the graph geodesic distance matrix,
- accept user-supplied starting coordinates,
- run either a compiled optimizer or an R fallback,
- return optional per-iteration trace and accepted frames.

## Objective and Kernel Reuse

The geodesic-KK and geodesic-MDS objectives differ only in weights and targets.

### Existing geodesic-KK form

Current GKK/LGKK uses

- target: `L0 * g_ij`
- weight: `K / g_ij^2`

### Geodesic-MDS form

Geodesic MDS uses

- target: `g_ij`
- weight: `1`

Therefore the path kernel is the same:

1. compute each stored path's embedded length,
2. form the residual against a target path length,
3. distribute that residual back across the path edges.

The first implementation can therefore reuse the structure of the current geodesic-KK path kernel directly.

## Detailed Implementation Plan

### Phase 1: additive implementation

This phase should be implemented now.

1. Add R helpers for deterministic symmetric `k`-NN graph construction and MST augmentation.
2. Add `grip.prepare.geodesic.mds()` on top of the existing full all-pairs geodesic-KK preparation code.
3. Add `grip.score.geodesic.mds()` using the same path cache fields as `grip.score.geodesic.kk()`, but with unit weights and no fitted scale.
4. Add a dedicated compiled optimizer with a standalone all-pairs shortest-path cache builder.
5. Add `grip.optimize.geodesic.mds()` that uses classical-MDS initialization and the compiled optimizer by default.
6. Add tests for exactness, determinism, and optimizer improvement.

### Phase 2: generic geodesic path kernel refactor

This phase is desirable but not required for the first working implementation.

Refactor the duplicated path machinery into a generic internal kernel parameterized by:

- pair targets,
- pair weights,
- path cache.

That would let:

- geodesic KK,
- landmark geodesic KK,
- geodesic MDS

share one internal evaluator in both R and C++.

## Connectivity Policy

The full geodesic-MDS path cache requires a connected graph.

For the first implementation, the default policy should be:

- `connect = "mst"`

That means:

1. build the symmetric `k`-NN graph,
2. test connectivity,
3. if disconnected, augment it with the Euclidean minimum spanning tree.

This keeps the prepared metric connected and practical for downstream optimization, while remaining close to the manuscript's discussion of graph-based connectivity correction.

## Initialization Policy

The default initialization should be classical MDS on the graph geodesic distance matrix:

- `cmdscale(as.dist(prepared$distance_matrix), k = dim)`

This follows the manuscript and places the optimizer near the Isomap solution.

If the embedding rank is lower than the requested output dimension, the initializer should pad missing columns with zeros.

## Dimensionality Support

The package's compiled geometry infrastructure is currently built around 2D and 3D points.

Therefore the first compiled geodesic-MDS implementation should support:

- `dim = 2`
- `dim = 3`

This is consistent with the rest of `grip`.

The underlying manuscript is dimension-agnostic, but a generic `d`-dimensional compiled kernel should be treated as a later extension.

## Testing Strategy

The first implementation should add targeted tests for:

1. exact zero stress on a weighted path realization,
2. deterministic `k`-NN + MST preparation,
3. optimizer decrease from a perturbed initialization,
4. consistency between graph-derived prepared objects and geodesic-MDS scoring.

## Risks and Guardrails

### Acceptable initial duplication

The compiled shortest-path and path-cache logic may initially duplicate small portions of the LGKK code in file-local helpers.

That is acceptable in phase 1 because:

- it keeps the new feature additive,
- it reduces regression risk in the existing GRIP engine,
- it gets a working GMDS path into the package quickly.

### Explicit limitations

The first implementation should document these constraints:

- only 2D and 3D optimization are compiled,
- `k`-NN construction is dense and intended for moderate `n`,
- the optimizer is gradient descent, not yet path-augmented SMACOF,
- the graph is fixed during optimization.

## Recommended Next Refactor After Phase 1

After the first implementation passes tests, the next cleanup should be:

1. move the deterministic shortest-path-tree logic into a shared compiled helper,
2. move all-pairs and sparse path-cache building into a shared compiled helper,
3. parameterize the compiled path objective by pair weights and targets,
4. let LGKK and GMDS call that common kernel.

That will align the codebase with the long-term design without making the first implementation harder than necessary.
