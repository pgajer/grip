# GMDS vs Classical MDS Comparison Test Suite Proposal

Date: 2026-03-31

This note proposes a phased comparison suite for geodesic MDS (GMDS) versus
classical MDS on the synthetic graph families documented in
[graph_families_generated_in_thread_2026-03-31.md](https://github.com/pgajer/grip/blob/main/dev/design/graph_families_generated_in_thread_2026-03-31.md).

## Goal

Build a comparison suite that answers three questions:

1. On which graph families does GMDS materially outperform classical MDS?
2. When are the two methods effectively tied?
3. How do topology, holes, curvature, irregular sampling, and ambient
   dimension change the comparison?

The suite should be incremental. Phase 1 should stay deliberately simple and
focus on mesh-like families before expanding to higher-genus, recursive, and
3D families.

## Proposed Baselines

### GMDS

Primary GMDS method:

- `grip.optimize.geodesic.mds()`
- default optimizer engine: `engine = "cpp"`
- default initialization: `init = "cmdscale"`

This is the method we want to evaluate as the `grip` geodesic-MDS
implementation.

### Classical MDS

Primary classical-MDS baseline:

- `stats::cmdscale(as.dist(Dg), k = dim, eig = TRUE)`

where `Dg` is the all-pairs graph-geodesic distance matrix on the weighted test
graph and `dim` is the target embedding dimension.

This is the classical MDS implementation I propose to use in the comparison
tests.

Why this choice:

- it is the standard classical-MDS implementation in base R,
- it is deterministic,
- it has no extra package dependency,
- it is already used as the default initializer for
  `grip.optimize.geodesic.mds()`,
- it gives a clean closed-form baseline rather than another iterative solver.

Scope note:

- I do not propose `smacof` or `isoMDS` as the primary baseline here because
  they are not classical MDS.
- If we later want a stronger non-classical baseline, we can add a separate
  metric-MDS phase without changing the primary classical-MDS comparison.
- For non-Euclidean dissimilarity matrices we should record the negative
  eigenvalue mass from `cmdscale`. As a later sensitivity analysis we can add
  `cmdscale(..., add = TRUE)`, but I would not make that the phase-1 baseline.

## Common Comparison Protocol

For each test case:

1. Generate a weighted graph bundle from a family helper in `R/graph_helpers.R`.
2. Compute the weighted graph-geodesic matrix `Dg`.
3. Fit a classical-MDS embedding with `stats::cmdscale()`.
4. Fit a GMDS embedding with `grip.optimize.geodesic.mds()`.
5. Score both embeddings under both objectives.
6. When ground-truth coordinates are available, compare each embedding to the
   known synthetic geometry.

## Primary Metrics

Every case should report:

- `gmds_raw_stress`: `sum_{i<j} (h_ij(Z) - g_ij)^2`
- `gmds_stress`: Stress-1 style normalized GMDS stress
- `euclidean_raw_stress`: `sum_{i<j} (||z_i-z_j|| - g_ij)^2`
- `euclidean_stress`: normalized classical stress against `Dg`
- `mean_abs_path_error`
- `mean_rel_path_error`
- Pearson correlation between `g_ij` and embedded Euclidean distances
- Spearman correlation between `g_ij` and embedded Euclidean distances
- runtime
- GMDS iteration count / accepted steps
- classical-MDS positive-eigen fraction and negative-eigen fraction

When a family supplies a true embedding, also report:

- Procrustes RMSE to `coords_surface`, `coords_solid`, or canonical coordinates
- edge-length RMSE against the true synthetic geometry
- local `k`-NN preservation in the true geometry

Important comparison rule:

- Both embeddings should be scored under the GMDS path objective and under the
  direct Euclidean-distance objective.

That gives a cleaner picture than judging classical MDS only by GMDS’s loss or
judging GMDS only by classical MDS’s loss.

## Normalization and Dimensions

Primary weighting convention:

- For `*.surface.graph()` and `*.solid.graph()` wrappers, use
  `normalize = "median"` in the main benchmark suite.

Rationale:

- this makes step sizes and stress magnitudes more comparable across families,
- it reduces arbitrary scale effects,
- it keeps the comparison focused on geometry and topology rather than raw
  unit differences.

Sensitivity convention:

- rerun a smaller subset with `normalize = "none"` later to make sure the
  conclusions are not only a normalization artifact.

Target dimensions:

- 2D for planar-layout comparisons and visualization-oriented tests,
- 3D for families with known 3D geometry,
- for selected families run both 2D and 3D, because the GMDS/MDS gap may be
  small in 3D but large in 2D.

## Test Architecture

The suite should have two layers.

### 1. `testthat` regression layer

Purpose:

- deterministic correctness and non-regression checks,
- small graphs only,
- quick enough for routine test runs.

What belongs here:

- smoke cases for each phase,
- a few direct GMDS-vs-MDS assertions where the expected ordering is stable,
- metric-computation checks,
- reproducibility checks.

### 2. benchmark/report layer

Purpose:

- larger parameter grids,
- timing,
- summary tables and plots,
- phase-by-phase evidence rather than binary pass/fail.

Suggested outputs:

- case-level CSV
- RDS with full metrics
- per-phase markdown summary
- representative plot grids

## Phased Rollout

## Phase 1: Meshes and Perforated Meshes

This should be the first implemented phase.

Families:

- `edges.mesh()`
- `mesh.surface.graph()`
- `edges.occupied.mesh()`
- `occupied.mesh.surface.graph()`
- `keep.periodic.holes()`
- `keep.staggered.windows()`
- `keep.slit.channels()`
- `keep.asymmetric.notches()`

Why start here:

- meshes are easy to reason about,
- they give us both control cases and failure cases,
- they let us separate curvature effects from hole/bottleneck effects,
- they are inexpensive enough for early benchmarking.

Phase-1 subgroups:

1. Flat regular meshes
2. Curved regular meshes
3. Flat perforated meshes
4. Curved perforated meshes

### Phase-1 case grid

Flat regular meshes:

- `edges.mesh(10, 10)`, `dim = 2`
- `edges.mesh(20, 20)`, `dim = 2`
- `edges.mesh(35, 35)`, `dim = 2`

Curved regular meshes:

- `mesh.surface.graph(15, 15, surface = "saddle", amplitude = 0.35)`, `dim = 2`
- `mesh.surface.graph(15, 15, surface = "paraboloid", amplitude = 0.35)`, `dim = 2`
- `mesh.surface.graph(15, 15, surface = "ripple", amplitude = 0.50, freq_u = 2, freq_v = 2)`, `dim = 2`
- repeat the same three with `dim = 3`

Flat perforated meshes:

- `occupied.mesh.surface.graph(keep.periodic.holes(25, 25), surface = "saddle", amplitude = 0)`, `dim = 2`
- `occupied.mesh.surface.graph(keep.staggered.windows(25, 25), surface = "saddle", amplitude = 0)`, `dim = 2`
- `occupied.mesh.surface.graph(keep.slit.channels(25, 25, orientation = "vertical"), surface = "saddle", amplitude = 0)`, `dim = 2`
- `occupied.mesh.surface.graph(keep.asymmetric.notches(25, 25), surface = "saddle", amplitude = 0)`, `dim = 2`

Curved perforated meshes:

- one `ripple` and one `paraboloid` variant for each of the four occupancy
  patterns above, `dim = 2`
- one smaller 3D subset, especially for `periodic.holes` and `slit.channels`

### Phase-1 hypotheses

- On flat regular meshes in 2D, classical MDS should be very competitive and
  may tie GMDS.
- On curved regular meshes in 3D, classical MDS should remain strong.
- On perforated meshes, especially slit/channel and window families, GMDS
  should improve graph-geodesic fidelity more clearly than classical MDS.
- The first robust GMDS win should appear in 2D perforated-mesh cases with
  bottlenecks and holes.

### Phase-1 regression assertions

Small stable `testthat` checks:

- GMDS and classical MDS both achieve near-zero stress on a tiny rectangular
  path-like mesh control.
- On a small perforated mesh, GMDS stress is strictly lower than classical-MDS
  stress under the GMDS path objective.
- On a flat small mesh, the GMDS advantage is small, confirming that the suite
  is not biased toward GMDS on every case.

## Phase 2: Wrapped Lattice Surfaces

Families:

- `cylinder.surface.graph()`
- `torus.surface.graph()`
- `sphere.surface.graph()`

Main questions:

- How much does periodic topology hurt classical MDS?
- Does GMDS better preserve long wrapped geodesics?
- Is the gap small in 3D but large in 2D?

Recommended cases:

- standard / mildly deformed / strongly wavy variants
- small and medium sample sizes
- both `dim = 2` and `dim = 3`

Expected behavior:

- cylinder: moderate GMDS advantage in 2D, smaller in 3D
- torus: stronger GMDS advantage, especially in 2D
- sphere: mixed, with less dramatic separation than torus

## Phase 3: Recursive Planar and Near-Planar Fractals

Families:

- `edges.recursive.mask.grid()`
- `sierpinski.carpet.surface.graph()`
- `vicsek.surface.graph()`
- `recursive.triangle.mask.surface.graph()`
- `sierpinski.triangle.surface.graph()`

Main questions:

- How well do the methods handle repeated bottlenecks?
- Does self-similar structure amplify classical-MDS distortions?
- Are bridge and hole motifs especially favorable to GMDS?

Expected behavior:

- GMDS should outperform classical MDS on carpet-like and bridge-like families
  where short Euclidean shortcuts are topologically wrong.

## Phase 4: Triangulated Manifolds with and without Boundary

Families:

- `triangulated.polyhedron.surface.graph()`
- `triangulated.annulus.surface.graph()`
- `triangulated.pair.of.pants.surface.graph()`

Main questions:

- How does GMDS behave on better-shaped triangulations than grid-derived
  families?
- Are boundary components and multiple holes enough to separate GMDS from
  classical MDS?

Expected behavior:

- small gap on closed low-distortion polyhedra in 3D,
- larger gap on annulus and pair-of-pants families in 2D.

## Phase 5: Irregular Point-Sampled Manifolds

Families:

- `irregular.annulus.surface.graph()`
- `irregular.sphere.surface.graph()`
- `irregular.pair.of.pants.surface.graph()`
- `irregular.torus.surface.graph()`
- `irregular.double.torus.surface.graph()`

Main questions:

- Are the GMDS gains robust to nonuniform sampling?
- Does irregular sampling destabilize classical MDS more than GMDS?
- How much harder is genus-2 than genus-1?

Expected behavior:

- irregular annulus and irregular torus should be strong mid-suite
  differentiators,
- irregular double torus should be one of the clearest late-phase GMDS wins.

## Phase 6: Recursive 3D Fractals and Porous Cubes

Families:

- `recursive.tetrahedron.mask.surface.graph()`
- `sierpinski.tetrahedron.surface.graph()`
- `recursive.cube.mask.surface.graph()`
- `menger.sponge.surface.graph()`
- `cube.periodic.tunnels.surface.graph()`
- `cube.asymmetric.cavities.surface.graph()`
- `cube.channel.network.surface.graph()`

Main questions:

- How do the methods compare on volumetric or cavity-rich 3D families?
- Does GMDS better preserve tunnel and cavity structure?
- Which porous-cube family is the strongest GMDS stress test?

Expected behavior:

- Menger sponge and channel-network cubes should strongly stress classical MDS,
- 3D GMDS should show clearer advantages than in simple smooth surfaces because
  the topology is harder and shortcuts are more harmful.

## Phase 7: Irregular 3D Solids and Intrinsic Weighted Trees

Families:

- `irregular.ball.solid.graph()`
- `irregular.shell.solid.graph()`
- `kary.tree.weighted.graph()`

Why keep this phase late:

- the solids are larger and more expensive,
- weighted trees are conceptually simple but use a different intrinsic regime
  from surface families.

Main questions:

- How much does GMDS help in layered nonuniform 3D solids?
- On weighted trees, does GMDS substantially reduce path distortion versus
  classical MDS?

Expected behavior:

- weighted trees should be a strong intrinsic control where GMDS has a clear
  objective advantage,
- irregular shell should be more challenging than irregular ball because of the
  cavity.

## Recommended Reporting Structure

For each phase, produce:

- a case table ranked by GMDS-minus-MDS stress improvement
- representative side-by-side plots
- a short summary of wins, ties, and failures
- a note on where classical MDS remains competitive

This matters because the suite should not become a one-sided “GMDS always wins”
artifact. The most useful result is a map of where each method is appropriate.

## Immediate Next Step

Implement only Phase 1 first.

Recommended immediate deliverables:

- a small `testthat` file with deterministic mesh/perforated-mesh regression
  checks,
- a benchmark script that writes a case-level CSV for the Phase-1 grid,
- a markdown summary with plot panels for the most representative cases.

## Short Answer on Classical MDS

The classical-MDS implementation I propose to use in these tests is:

- `stats::cmdscale(as.dist(Dg), k = dim, eig = TRUE)`

That should be the primary baseline throughout the suite.
