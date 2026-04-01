# Generic Spring-Repulsion GMDS Design

Date: 2026-04-01

## Goal

Design a **generic graph** replacement for plain GMDS that:

- keeps the geodesic-path objective as the core signal,
- discourages graph-distant vertices from collapsing together,
- prevents unlimited stretching by attaching a spring to each graph edge,
- does **not** rely on mesh-specific notions such as cells, triangle orientation, or surface normals.

This note treats the current `grip` GMDS implementation as the starting point:

- path-stress and R optimizer in [/Users/pgajer/current_projects/grip/R/grip_quality.R](/Users/pgajer/current_projects/grip/R/grip_quality.R)
- flattened C++ kernels in [/Users/pgajer/current_projects/grip/src/geodesic_mds_rcpp.cpp](/Users/pgajer/current_projects/grip/src/geodesic_mds_rcpp.cpp)
- current prepared geodesic object in [/Users/pgajer/current_projects/grip/R/grip_quality.R#L2956](/Users/pgajer/current_projects/grip/R/grip_quality.R#L2956)

## Why Change The Objective

The current GMDS objective only penalizes mismatch between graph geodesic targets and embedded cached-path lengths:

`E_path(Z) = 1/2 * sum_{i<j} (h_ij(Z) - g_ij)^2`

where the implemented embedded path length is

`h_ij(Z) = sum_{e in cached_path(i,j)} coeff_ij(e) * ||z_u - z_v||_eps`

This is exactly what the current R and C++ code does:

- R path term accumulation in [/Users/pgajer/current_projects/grip/R/grip_quality.R#L2760](/Users/pgajer/current_projects/grip/R/grip_quality.R#L2760)
- C++ flat path term accumulation in [/Users/pgajer/current_projects/grip/src/geodesic_mds_rcpp.cpp#L834](/Users/pgajer/current_projects/grip/src/geodesic_mds_rcpp.cpp#L834)

That term controls **path lengths**, but it does not directly control:

- whether graph-distant vertices collapse together in the embedding,
- whether different parts of the graph crowd into the same ambient-space region,
- whether local graph edges stretch far beyond their original target lengths.

So the generic fix should add:

1. an **edge spring term** to keep local graph geometry from stretching too far,
2. a **graph-aware repulsion term** to keep graph-distant vertices from collapsing too close,
3. an optional anchor continuation term for stable optimization.

## Proposed Generic Objective

Let `G = (V, E)` be a weighted graph with edge targets `l_uv > 0`.

Let:

- `g_ij` be the graph geodesic distance,
- `hbar_ij(Z)` be the current cached-path or tie-averaged embedded path length,
- `d_ij(Z) = ||z_i - z_j||_eps`,
- `l_uv` be the target embedded length for graph edge `(u,v)`, usually the graph edge weight.

The proposed generic objective is:

`E_total(Z) = alpha * E_path(Z) + beta * E_edge(Z) + gamma * E_rep(Z) + lambda * E_anchor(Z)`

with

`E_path(Z) = 1/2 * sum_{(i,j) in P} w_ij^path * (hbar_ij(Z) - g_ij)^2`

`E_edge(Z) = 1/2 * sum_{(u,v) in E} w_uv^edge * (d_uv(Z) - l_uv)^2`

`E_rep(Z) = 1/2 * sum_{(i,j) in R} w_ij^rep * [tau_ij - d_ij(Z)]_+^2`

`E_anchor(Z) = ||Z - A||_F^2`

where:

- `P` is the set of path-stress pairs, typically all pairs or a sampled subset,
- `R` is a set of repulsion pairs chosen from graph-distant, nonadjacent pairs,
- `[x]_+ = max(x, 0)`.

## Interpretation

This is the graph-generic spring picture:

- `E_path` says: preserve the graph's intrinsic path geometry.
- `E_edge` says: every graph edge is a spring and resists excessive stretching or shrinking.
- `E_rep` says: graph-distant vertices should not collapse into the same ambient-space region.
- `E_anchor` says: start from a sensible global shape and relax away from it gradually.

This is no longer "pure" GMDS in the narrow manuscript sense, but it is a natural graph-generic extension of GMDS rather than a mesh-specific patch.

## Recommended Default Parameterization

### Scaling

First normalize graph edge targets so that the median edge target is `1`.

That makes the following defaults dimensionless and easier to tune across graphs.

### Path Term

Default:

- `alpha = 1`
- `w_ij^path = 1`

Rationale:

- preserve the current GMDS meaning,
- do not add KK-style inverse-square weighting in the first version.

### Edge Spring Term

Default:

- `beta = 0.25`
- `w_uv^edge = 1`
- `l_uv = original graph edge weight`

Rationale:

- this is strong enough to resist runaway distortion of local graph edges,
- but weaker than the geodesic path term so the method still behaves like a geodesic layout method rather than a plain spring embedder.

### Repulsion Term

Repulsion should be **graph-aware**, not all-pairs Coulomb repulsion.

Recommended active set:

- `R = {(i,j): graph_hop(i,j) >= 3}` if an unweighted hop metric is available, or
- `R = {(i,j): g_ij >= q_0.6(g)}` using the upper 40% of graph distances.

Recommended lower-bound target:

`tau_ij = rho * min(g_ij, q_0.9(g))`

with default:

- `rho = 0.20`
- `gamma = 0.10`
- `w_ij^rep = 1`

Rationale:

- only graph-distant pairs repel,
- the lower bound scales with graph distance,
- clipping at the 90th percentile avoids impossible long-range lower bounds dominating the energy.

### Anchor Term

Optional but recommended for optimization stability:

- anchor `A = cmdscale(as.dist(D_g))`
- start `lambda = 0.10`
- end `lambda_end = 0.01`
- linear continuation over iterations

Rationale:

- the anchor stabilizes the early trajectory,
- but the final layout is not forced to stay too close to classical MDS.

## Continuation Schedule

Recommended default optimization schedule:

1. Warm start:
   - `lambda = 0.10`
   - `gamma = 0.02`
   - `beta = 0.25`

2. Stretch phase:
   - linearly decrease `lambda` to `0.01`
   - linearly increase `gamma` to `0.10`
   - keep `beta` fixed

3. Polish phase:
   - keep `lambda = 0.01`
   - keep `gamma = 0.10`
   - optionally reduce `initial_step`

A simple implementation is:

- `lambda_t = lambda_0 + (t/T) * (lambda_T - lambda_0)`
- `gamma_t = gamma_0 + (t/T) * (gamma_T - gamma_0)`
- `beta_t = beta`

This fits naturally with the existing schedule machinery already used for anchor and smoothness weights in [/Users/pgajer/current_projects/grip/R/grip_quality.R#L3264](/Users/pgajer/current_projects/grip/R/grip_quality.R#L3264).

## Gradients

### Path Term

Keep the current gradient unchanged:

`grad_path(u) = sum_{(i,j): e=(u,v) in path_ij} (hbar_ij - g_ij) * coeff_ij(e) * (z_u - z_v) / ||z_u - z_v||_eps`

This is already implemented in:

- [/Users/pgajer/current_projects/grip/R/grip_quality.R#L2772](/Users/pgajer/current_projects/grip/R/grip_quality.R#L2772)
- [/Users/pgajer/current_projects/grip/src/geodesic_mds_rcpp.cpp#L849](/Users/pgajer/current_projects/grip/src/geodesic_mds_rcpp.cpp#L849)

### Edge Spring Term

For edge `(u,v)`:

`grad_edge(u) = beta * w_uv^edge * (d_uv - l_uv) * (z_u - z_v) / d_uv`

`grad_edge(v) = -grad_edge(u)`

This is structurally identical to the current path-edge derivative, so it is cheap to add.

### Repulsion Term

For an active repulsion pair `(i,j)` with `d_ij < tau_ij`:

`grad_rep(i) = -gamma * w_ij^rep * (tau_ij - d_ij) * (z_i - z_j) / d_ij`

`grad_rep(j) = -grad_rep(i)`

If `d_ij >= tau_ij`, the pair contributes zero.

This is a one-sided spring that only pushes apart when graph-distant vertices get too close.

## Why This Stays Generic

This formulation depends only on:

- graph edges,
- graph shortest-path distances,
- graph-distant vertex pairs.

It does **not** depend on:

- mesh cells,
- triangulations,
- surface normals,
- cell orientations,
- manifold assumptions.

So it is meaningful for:

- meshes,
- trees,
- social / biological graphs,
- k-NN graphs from data,
- generic weighted graphs.

## How To Integrate This Into `grip`

### 1. Extend The Prepared Object

Add the following fields to the prepared GMDS object produced in [/Users/pgajer/current_projects/grip/R/grip_quality.R#L2989](/Users/pgajer/current_projects/grip/R/grip_quality.R#L2989):

- `graph_edge_matrix`
- `graph_edge_target`
- `flat_graph_edge_u`
- `flat_graph_edge_v`
- `flat_graph_edge_target`
- `repulsion_pair_matrix`
- `repulsion_pair_target`
- `flat_repulsion_u`
- `flat_repulsion_v`
- `flat_repulsion_target`

The graph edge target can be read directly from `base$edges` and `base$weight_list`.

The repulsion pairs can be built from `prepared$distance_matrix` and `prepared$adj_list`.

### 2. Add New R-Level API Arguments

Extend [/Users/pgajer/current_projects/grip/R/grip_quality.R#L3176](/Users/pgajer/current_projects/grip/R/grip_quality.R#L3176) with:

- `edge_spring_weight`
- `edge_spring_weight_end`
- `edge_spring_continuation`
- `repulsion_weight`
- `repulsion_weight_end`
- `repulsion_continuation`
- `repulsion_quantile`
- `repulsion_scale`
- `repulsion_cap_quantile`
- `repulsion_hop_min`

Suggested first defaults:

- `edge_spring_weight = 0.25`
- `edge_spring_weight_end = 0.25`
- `repulsion_weight = 0.02`
- `repulsion_weight_end = 0.10`
- `repulsion_quantile = 0.60`
- `repulsion_scale = 0.20`
- `repulsion_cap_quantile = 0.90`

### 3. R Prototype First

Add the new terms first in the R evaluation path:

- extend [/Users/pgajer/current_projects/grip/R/grip_quality.R#L2723](/Users/pgajer/current_projects/grip/R/grip_quality.R#L2723)
- extend [/Users/pgajer/current_projects/grip/R/grip_quality.R#L2807](/Users/pgajer/current_projects/grip/R/grip_quality.R#L2807)
- extend [/Users/pgajer/current_projects/grip/R/grip_quality.R#L2902](/Users/pgajer/current_projects/grip/R/grip_quality.R#L2902)

This gives:

- simple debugging,
- finite-difference gradient checks,
- a reference implementation for the C++ kernel.

### 4. Then Add The C++ Flat Kernel

Extend the flat state evaluator in [/Users/pgajer/current_projects/grip/src/geodesic_mds_rcpp.cpp#L911](/Users/pgajer/current_projects/grip/src/geodesic_mds_rcpp.cpp#L911):

- keep `accumulate_flat_pair_range()` for path stress,
- add `accumulate_flat_edge_springs()`,
- add `accumulate_flat_repulsion()`,
- add schedules for edge and repulsion weights parallel to `resolve_anchor_schedule()` in [/Users/pgajer/current_projects/grip/src/geodesic_mds_rcpp.cpp#L707](/Users/pgajer/current_projects/grip/src/geodesic_mds_rcpp.cpp#L707).

This should plug into the current flattened optimizer path used by
`grip_optimize_geodesic_mds_flat_cpp()`.

### 5. Keep The Old Objective As A Named Mode

Do not replace current GMDS silently.

Instead expose:

- `objective = "plain_gmds"`
- `objective = "spring_repulsion_gmds"`

This keeps comparability with the manuscript objective.

## Recommended Validation Plan

### Phase 1: Algebra And Unit Tests

Add tests for:

- edge-spring gradient finite differences,
- repulsion gradient finite differences,
- R vs C++ equality on a tiny graph,
- zero-contribution behavior when no repulsion pair is active.

### Phase 2: Generic Graph Families

Test on:

- path graphs,
- cycles,
- grids,
- trees,
- clustered random geometric graphs,
- symmetric k-NN graphs from point clouds.

Primary metrics:

- path stress,
- edge spring distortion,
- graph-distant crowding rate,
- Procrustes error where ground truth exists,
- symmetry / regularity where expected.

### Phase 3: Pathology Regression

Re-run the problematic cases:

- flat orthogonal meshes,
- paraboloids,
- hemispheres,
- ripple surfaces in 3D.

The question is not whether stress becomes minimal at all costs. The question is whether the new objective preserves low path stress **without** the dramatic crowding and fold-over seen under plain GMDS.

## Recommended First Implementation Scope

Keep the first implementation small:

1. all-pairs path term exactly as now,
2. edge spring term on all graph edges,
3. repulsion on graph-distance-selected vertex pairs,
4. anchor continuation retained,
5. no mesh-specific penalties.

That is already enough to test the core hypothesis:

> The missing generic ingredient in GMDS is not mesh geometry, but graph-aware self-separation plus local edge stiffness.

## Bottom Line

The most practical generic next objective for `grip` is:

`E_total(Z) = E_path(Z) + 0.25 E_edge(Z) + gamma_t E_rep(Z) + lambda_t E_anchor(Z)`

with:

- `gamma_t` increasing from `0.02` to `0.10`,
- `lambda_t` decreasing from `0.10` to `0.01`,
- graph-distant repulsion pairs chosen from the upper-distance quantiles,
- edge targets equal to the original graph edge weights.

This keeps the spirit of GMDS, stays graph-generic, and directly addresses the collapse mode that plain path stress does not see.
