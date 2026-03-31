# Hybrid GMDS with MDS Tether and Tie-Averaged Shortest Paths

## Goal

This note records the Phase 1 follow-up changes motivated by the orthogonal-mesh diagnostics:

1. add an MDS-tethered GMDS mode,
2. add a continuation mode that gradually relaxes the tether,
3. add tie-averaged multi-path GMDS so rectangular meshes are not forced through one deterministic shortest-path selection.

The key point is that these are complementary:

- the tie-averaged path objective removes symmetry-breaking artifacts caused by selecting one arbitrary shortest path out of many tied shortest paths,
- the MDS tether discourages aesthetically undesirable fold-overs even when the GMDS term alone would permit them.

## Base GMDS Objective

For an embedding `Z = (z_1, ..., z_n)` and graph-geodesic distances `g_ij`, the fixed-path GMDS objective is

`E_GMDS(Z) = (1/2) sum_{i<j} (h_ij(Z) - g_ij)^2`

where `h_ij(Z)` is the embedded length of the cached graph-shortest path between vertices `i` and `j`.

For a cached path `p = (v_0, ..., v_m)`, its embedded length is

`h_p(Z) = sum_{r=1}^m || z_{v_r} - z_{v_{r-1}} ||`.

## Tie-Averaged Multi-Path GMDS

When there are several tied graph-shortest paths between `i` and `j`, let `P_ij` denote the full tied family.

Instead of choosing one deterministic path, define the uniform tie-averaged embedded path length

`bar(h)_ij(Z) = (1 / |P_ij|) sum_{p in P_ij} h_p(Z)`.

This can be written edgewise as

`bar(h)_ij(Z) = sum_{e in E} pi_ij(e) * ell_e(Z)`

where

- `ell_e(Z)` is the embedded Euclidean length of edge `e`,
- `pi_ij(e)` is the probability that edge `e` lies on a uniformly random shortest path in `P_ij`.

The tie-averaged GMDS objective is then

`E_avg(Z) = (1/2) sum_{i<j} (bar(h)_ij(Z) - g_ij)^2`.

The implementation computes `pi_ij(e)` exactly from the shortest-path DAG using forward and backward path counts, so no path enumeration is required.

## MDS-Tethered Hybrid Objective

Let `Z_cmd` be the classical-MDS embedding of the graph-geodesic distance matrix. The tethered objective is

`E_hyb(Z; lambda) = E_GMDS_or_avg(Z) + lambda || Z - Z_cmd ||_F^2`

where `|| . ||_F` is the Frobenius norm.

The tether gradient is

`nabla_Z [ lambda || Z - Z_cmd ||_F^2 ] = 2 lambda (Z - Z_cmd)`.

This term does not replace the GMDS objective. It acts as a shape prior that penalizes drift away from the globally symmetric `cmdscale` solution.

## Continuation Schedule

Instead of keeping `lambda` fixed, continuation uses an iteration-dependent tether

`E_t(Z) = E_GMDS_or_avg(Z) + lambda_t || Z - Z_cmd ||_F^2`,

with `t = 0, 1, ..., T`.

Two schedules are implemented.

### Linear continuation

`lambda_t = (1 - s_t) lambda_0 + s_t lambda_T`

with

`s_t = t / T`.

This is the recommended schedule when the final tether should be exactly zero.

### Geometric continuation

`lambda_t = lambda_0 (lambda_T / lambda_0)^(t / T)`

for `lambda_0 > 0` and `lambda_T > 0`.

This keeps the tether strong early, then relaxes it multiplicatively. It is useful when the tether should remain positive throughout the run.

## Practical Interpretation

The intended workflow is:

1. initialize with `cmdscale`,
2. start with a moderate or strong tether to preserve the good global shape,
3. gradually relax the tether so GMDS can improve intrinsic distances,
4. use `tie_mode = "average"` on orthogonal rectangular meshes so the objective respects all tied Manhattan shortest paths rather than one arbitrary route.

## Implementation Summary

The package changes expose:

- `tie_mode = c("single", "average")` in `grip.prepare.geodesic.kk()`, `grip.prepare.geodesic.mds()`, `grip.score.geodesic.mds()`, and `grip.optimize.geodesic.mds()`,
- `anchor_mode = c("none", "cmdscale", "initial", "user")`,
- `anchor_coords`,
- `anchor_weight`,
- `anchor_weight_end`,
- `continuation = c("constant", "linear", "geometric")`.

The compiled optimizer now consumes the prepared path cache directly, so it can optimize:

- ordinary single-path GMDS,
- exact tie-averaged GMDS,
- tethered GMDS,
- continuation GMDS,
- and combinations of the above.

## Rectangular-Mesh Diagnostic

On the orthogonal flat `10 x 10` mesh:

- single-path GMDS at 50 iterations gave `rho ≈ 0.213`,
- tie-averaged GMDS at 50 iterations gave `rho ≈ 1.9e-4`.

So for rectangular meshes, the main artifact was not an unavoidable limitation of orthogonal connectivity. It was largely the consequence of collapsing a large tied shortest-path family down to one deterministic route per vertex pair.
