# GMDS Algorithm, Implementation, and Pathology Overview

## Purpose of This Note

This note records the current state of the geodesic-MDS (GMDS) work in `grip`.
It has four goals:

1. describe the GMDS algorithm now implemented in the package,
2. point to the files where each part is implemented,
3. summarize the experiment sequence and point to the corresponding reports,
4. make explicit the main pathology questions that motivated those experiments.

The original mathematical motivation lives in the companion manuscript
[`geodesic_mds.tex`](/Users/pgajer/current_projects/geodesic_MDS/manuscript/geodesic_mds.tex).
The implementation work described here lives in the `grip` package.

## Central Questions

The experiments were driven by a small set of concrete questions.

### Question 1: Why do simple rectangular meshes look wrong under naive GMDS?

On flat orthogonal meshes, single-path GMDS could achieve very low GMDS stress while producing embeddings that looked asymmetric, faceted, or folded relative to the obvious square reference geometry.

There turned out to be two different sources for that behavior:

- **Shortest-path tie breaking.**
  On an orthogonal grid there are many tied shortest paths between the same pair of vertices.
  If GMDS uses one deterministic chosen path per pair, the objective itself becomes asymmetric.

- **Metric mismatch.**
  Even without arbitrary tie breaking, the graph-geodesic metric on a 4-neighbor grid is much closer to an `L1` / Manhattan metric than to the Euclidean `L2` metric of the square parameter domain.
  So a low-stress geodesic embedding does not automatically have to look like the square parameter plane.

### Question 2: Are the curved-surface singularities mainly a convergence problem?

On paraboloid meshes, classical MDS was smooth but globally distorted, whereas GMDS reduced geodesic stress strongly but often created a visually unsettling fold or singular concentration of curvature.

The core question was:

`Is this mainly a gradient-descent tuning problem, or is the GMDS objective itself permitting the bad geometry?`

### Question 3: If the objective is the problem, what kind of extra structure is missing?

Several hypotheses were tested:

- a stronger tether to the classical-MDS shape,
- a first-order local smoothness regularizer,
- a second-order rectangular-grid bending penalty.

The underlying question was:

`What geometric information is absent from pure geodesic path stress?`

### Question 4: Which pathologies are intrinsic, and which are artifacts of the test setup?

Two examples made this especially clear:

- **Flat orthogonal meshes:** some artifacts were caused by deterministic single-path selection and were largely fixed by tie-averaging.
- **Ripple surfaces in 2D:** a faithful 2D embedding is impossible in principle, so poor-looking 2D ripple plots mostly demonstrate a dimensionality limitation, not a GMDS bug.

## GMDS Algorithm

### 1. Input Graph and Graph-Geodesic Distances

GMDS works on a weighted graph `G = (V, E)`.
The graph can be supplied directly, or built from data as a symmetric weighted `k`-nearest-neighbor graph, optionally augmented to ensure connectivity.

The algorithm first computes:

- the all-pairs graph-geodesic distances `g_ij`,
- a representation of shortest paths between every vertex pair,
- and a cache that makes repeated objective and gradient evaluation practical.

In `grip`, this graph-preparation layer is implemented primarily in:

- [`R/grip_quality.R`](/Users/pgajer/current_projects/grip/R/grip_quality.R)

The public entry points are:

- [`grip.prepare.geodesic.kk()`](/Users/pgajer/current_projects/grip/R/grip_quality.R)
- [`grip.prepare.geodesic.mds()`](/Users/pgajer/current_projects/grip/R/grip_quality.R)

### 2. Base Fixed-Path GMDS Objective

For an embedding `Z = (z_1, ..., z_n)`, the basic fixed-path GMDS objective is

`E_GMDS(Z) = (1/2) * sum_{i<j} (h_ij(Z) - g_ij)^2`

where:

- `g_ij` is the graph-geodesic distance between vertices `i` and `j`,
- `h_ij(Z)` is the embedded length of the cached graph-shortest path between `i` and `j`.

If the cached shortest path from `i` to `j` is `p = (v_0, ..., v_m)`, then

`h_p(Z) = sum_{r=1}^m || z_{v_r} - z_{v_{r-1}} ||`.

This is the geometric core of GMDS: preserve graph-geodesic path lengths in the embedding.

### 3. Tie-Averaged Multi-Path GMDS

The first major refinement was to replace deterministic single-path caching with exact tie-averaging over all tied shortest paths.

For a pair `(i, j)` with a tied shortest-path family `P_ij`, define

`hbar_ij(Z) = (1 / |P_ij|) * sum_{p in P_ij} h_p(Z)`.

This can also be written edgewise as

`hbar_ij(Z) = sum_e pi_ij(e) * ell_e(Z)`

where:

- `ell_e(Z)` is the embedded Euclidean length of edge `e`,
- `pi_ij(e)` is the probability that edge `e` belongs to a uniformly random shortest path from `i` to `j`.

The resulting objective is

`E_avg(Z) = (1/2) * sum_{i<j} (hbar_ij(Z) - g_ij)^2`.

This was crucial for orthogonal meshes because it removed symmetry-breaking artifacts caused by choosing one arbitrary Manhattan route per pair.

The relevant implementation files are:

- [`R/grip_quality.R`](/Users/pgajer/current_projects/grip/R/grip_quality.R)
- [`src/geodesic_mds_rcpp.cpp`](/Users/pgajer/current_projects/grip/src/geodesic_mds_rcpp.cpp)

The main public control is:

- `tie_mode = c("single", "average")`

### 4. Classical-MDS Initialization

All later GMDS experiments use classical MDS on the graph-geodesic distance matrix as the default initialization.

The baseline embedding is:

`Z_cmd = cmdscale(D_g)`

where `D_g` is the all-pairs graph-geodesic distance matrix.

In practice, `grip` uses the internal helper `grip.classical.mds.embedding`, which wraps the classical-MDS baseline consistently for the benchmark scripts.

Implementation:

- [`R/grip_quality.R`](/Users/pgajer/current_projects/grip/R/grip_quality.R)

Benchmark scripts consistently start from:

- `grip:::grip.classical.mds.embedding(prepared, dim = ...)`

### 5. MDS-Tethered Hybrid GMDS

To keep the optimizer near the globally smooth classical-MDS shape, a quadratic anchor term was added:

`E_lambda(Z) = E_GMDS_or_avg(Z) + lambda * || Z - Z_cmd ||_F^2`

This does not replace GMDS. It adds a shape prior that penalizes large drift away from the MDS embedding.

The tether can be:

- absent,
- constant,
- or scheduled over the course of optimization.

Public controls:

- `anchor_mode = c("none", "cmdscale", "initial", "user")`
- `anchor_coords`
- `anchor_weight`
- `anchor_weight_end`
- `continuation = c("constant", "linear", "geometric")`

Implementation:

- [`R/grip_quality.R`](/Users/pgajer/current_projects/grip/R/grip_quality.R)
- [`src/geodesic_mds_rcpp.cpp`](/Users/pgajer/current_projects/grip/src/geodesic_mds_rcpp.cpp)

The detailed design note is:

- [`geodesic_mds_hybrid_tether_and_tie_average_2026-03-31.md`](/Users/pgajer/current_projects/grip/dev/design/geodesic_mds_hybrid_tether_and_tie_average_2026-03-31.md)

### 6. First-Order Smoothness Regularization

The first regularization experiment added a local neighborhood-average penalty:

`E_{lambda, mu}(Z) = E_GMDS_or_avg(Z) + lambda * || Z - Z_cmd ||_F^2 + mu * sum_i || z_i - mean_{j in N(i)} z_j ||^2`

Interpretation:

- `lambda` preserves the global MDS shape,
- `mu` discourages local departures from neighborhood smoothness.

Implementation:

- base smoothness support in [`R/grip_quality.R`](/Users/pgajer/current_projects/grip/R/grip_quality.R)
- compiled support in [`src/geodesic_mds_rcpp.cpp`](/Users/pgajer/current_projects/grip/src/geodesic_mds_rcpp.cpp)

Public controls:

- `smoothness_weight`
- `smoothness_weight_end`
- `smoothness_continuation`

### 7. Second-Order Rectangular-Grid Bending Penalty

The second regularization experiment targeted rectangular meshes more directly by penalizing row-and-column second differences.

For a rectangular-grid stencil `(a, b, c)`, the discrete bending residual is

`r_abc(Z) = z_a - 2 z_b + z_c`.

The bending penalty is

`R_bend(Z) = (1 / |S|) * sum_{(a,b,c) in S} || r_abc(Z) ||^2`

and the hybrid objective is

`E_{lambda, beta}(Z) = E_GMDS_or_avg(Z) + lambda * || Z - Z_cmd ||_F^2 + beta * R_bend(Z)`.

This was meant to penalize concentrated second-order curvature more directly than the first-order neighborhood-average term.

Implementation files:

- stencil helpers: [`R/zz_geodesic_mds_bending.R`](/Users/pgajer/current_projects/grip/R/zz_geodesic_mds_bending.R)
- bending statistics: [`R/zz_geodesic_mds_bending_stats.R`](/Users/pgajer/current_projects/grip/R/zz_geodesic_mds_bending_stats.R)
- bending-aware score / optimize overrides:
  - [`R/zz_geodesic_mds_bending_core.R`](/Users/pgajer/current_projects/grip/R/zz_geodesic_mds_bending_core.R)
  - [`R/zz_geodesic_mds_bending_score.R`](/Users/pgajer/current_projects/grip/R/zz_geodesic_mds_bending_score.R)
  - [`R/zz_geodesic_mds_bending_opt_r.R`](/Users/pgajer/current_projects/grip/R/zz_geodesic_mds_bending_opt_r.R)
  - [`R/zz_geodesic_mds_bending_optimize.R`](/Users/pgajer/current_projects/grip/R/zz_geodesic_mds_bending_optimize.R)
- compiled optimizer:
  - [`src/geodesic_mds_bending_rcpp.cpp`](/Users/pgajer/current_projects/grip/src/geodesic_mds_bending_rcpp.cpp)
- generated exports:
  - [`R/RcppExports.R`](/Users/pgajer/current_projects/grip/R/RcppExports.R)
  - [`src/RcppExports.cpp`](/Users/pgajer/current_projects/grip/src/RcppExports.cpp)

Public controls:

- `bending_stencils`
- `bending_weight`
- `bending_weight_end`
- `bending_continuation`

### 8. Optimization Method

The optimizer is not SMACOF. It is deterministic gradient descent with Armijo backtracking.

At each iteration:

1. evaluate the current objective and gradient,
2. try a step with the current trial step size,
3. shrink the step until Armijo decrease is satisfied,
4. accept the new coordinates,
5. optionally recenter the embedding.

This is implemented in two layers:

- an R path for validation and fallback,
- a compiled C++ path for the real benchmark runs.

Key files:

- base R path and scheduling logic:
  - [`R/grip_quality.R`](/Users/pgajer/current_projects/grip/R/grip_quality.R)
- compiled flat-cache optimizer with tie-averaging, anchors, and smoothness:
  - [`src/geodesic_mds_rcpp.cpp`](/Users/pgajer/current_projects/grip/src/geodesic_mds_rcpp.cpp)
- compiled flat-cache optimizer with bending:
  - [`src/geodesic_mds_bending_rcpp.cpp`](/Users/pgajer/current_projects/grip/src/geodesic_mds_bending_rcpp.cpp)

## Implementation Map

### Core Public API

The main public functions are:

- [`grip.prepare.geodesic.kk()`](/Users/pgajer/current_projects/grip/R/grip_quality.R)
- [`grip.prepare.geodesic.mds()`](/Users/pgajer/current_projects/grip/R/grip_quality.R)
- [`grip.score.geodesic.mds()`](/Users/pgajer/current_projects/grip/R/grip_quality.R)
- [`grip.optimize.geodesic.mds()`](/Users/pgajer/current_projects/grip/R/grip_quality.R)

Conceptually:

- `prepare` builds graph distances and shortest-path caches,
- `score` evaluates an embedding,
- `optimize` runs the Armijo solver.

### Base R Implementation

Most of the base GMDS logic lives in:

- [`R/grip_quality.R`](/Users/pgajer/current_projects/grip/R/grip_quality.R)

This file contains:

- graph construction from data,
- shortest-path preparation,
- tie-averaged cache building,
- classical-MDS initialization,
- anchor handling and continuation schedules,
- score calculation,
- the R optimizer,
- and the first-order smoothness support.

### Bending Overlay

The bending implementation was added as a clean overlay rather than a wholesale rewrite of the base file. The `zz_*.R` files load after `grip_quality.R` and extend or override the necessary functions.

That bending layer lives in:

- [`R/zz_geodesic_mds_bending.R`](/Users/pgajer/current_projects/grip/R/zz_geodesic_mds_bending.R)
- [`R/zz_geodesic_mds_bending_core.R`](/Users/pgajer/current_projects/grip/R/zz_geodesic_mds_bending_core.R)
- [`R/zz_geodesic_mds_bending_score.R`](/Users/pgajer/current_projects/grip/R/zz_geodesic_mds_bending_score.R)
- [`R/zz_geodesic_mds_bending_opt_r.R`](/Users/pgajer/current_projects/grip/R/zz_geodesic_mds_bending_opt_r.R)
- [`R/zz_geodesic_mds_bending_optimize.R`](/Users/pgajer/current_projects/grip/R/zz_geodesic_mds_bending_optimize.R)
- [`R/zz_geodesic_mds_bending_stats.R`](/Users/pgajer/current_projects/grip/R/zz_geodesic_mds_bending_stats.R)

### Compiled Kernels

The compiled path is what all serious benchmarks use.

Files:

- base flat-cache optimizer and tie-averaged cache support:
  - [`src/geodesic_mds_rcpp.cpp`](/Users/pgajer/current_projects/grip/src/geodesic_mds_rcpp.cpp)
- bending-aware optimizer:
  - [`src/geodesic_mds_bending_rcpp.cpp`](/Users/pgajer/current_projects/grip/src/geodesic_mds_bending_rcpp.cpp)

Generated interface files:

- [`R/RcppExports.R`](/Users/pgajer/current_projects/grip/R/RcppExports.R)
- [`src/RcppExports.cpp`](/Users/pgajer/current_projects/grip/src/RcppExports.cpp)

### Tests

The focused regression coverage is in:

- base GMDS: [`tests/testthat/test-geodesic-mds.R`](/Users/pgajer/current_projects/grip/tests/testthat/test-geodesic-mds.R)
- tie-averaging and hybrid behavior: [`tests/testthat/test-geodesic-mds-hybrid.R`](/Users/pgajer/current_projects/grip/tests/testthat/test-geodesic-mds-hybrid.R)
- smoothness regularization: [`tests/testthat/test-geodesic-mds-smoothness.R`](/Users/pgajer/current_projects/grip/tests/testthat/test-geodesic-mds-smoothness.R)
- bending regularization: [`tests/testthat/test-geodesic-mds-bending.R`](/Users/pgajer/current_projects/grip/tests/testthat/test-geodesic-mds-bending.R)
- Phase 1 GMDS/MDS comparisons: [`tests/testthat/test-gmds-mds-phase1.R`](/Users/pgajer/current_projects/grip/tests/testthat/test-gmds-mds-phase1.R)

### Benchmark and Report Generators

The main experiment drivers are:

- Phase 1 cross-family comparison:
  - [`tools/benchmark-gmds-mds-phase1.R`](/Users/pgajer/current_projects/grip/tools/benchmark-gmds-mds-phase1.R)
- paraboloid pathology report:
  - [`tools/benchmark-paraboloid-gmds-pathology.R`](/Users/pgajer/current_projects/grip/tools/benchmark-paraboloid-gmds-pathology.R)
- paraboloid first-order regularization report:
  - [`tools/benchmark-paraboloid-gmds-regularization.R`](/Users/pgajer/current_projects/grip/tools/benchmark-paraboloid-gmds-regularization.R)
- paraboloid second-order bending report:
  - [`tools/benchmark-paraboloid-gmds-bending.R`](/Users/pgajer/current_projects/grip/tools/benchmark-paraboloid-gmds-bending.R)

## Experiment Sequence

## Phase 0: Infrastructure and Test-Suite Design

Before the pathology studies, the first task was to decide how GMDS should fit into `grip` and how to compare it to classical MDS in a reproducible way.

Key design notes:

- infrastructure design:
  - [`geodesic_mds_infrastructure_design_2026-03-31.md`](/Users/pgajer/current_projects/grip/dev/design/geodesic_mds_infrastructure_design_2026-03-31.md)
- infrastructure action plan:
  - [`geodesic_mds_infrastructure_action_plan_2026-03-31.md`](/Users/pgajer/current_projects/grip/dev/design/geodesic_mds_infrastructure_action_plan_2026-03-31.md)
- cross-family comparison test-suite design:
  - [`gmds_mds_comparison_test_suite_2026-03-31.md`](/Users/pgajer/current_projects/grip/dev/design/gmds_mds_comparison_test_suite_2026-03-31.md)

These notes established the main experimental baseline:

- use classical MDS as the baseline,
- compare against GMDS on synthetic graph families,
- and focus first on simple mesh families where visual failures are easy to see.

## Phase 1: GMDS vs Classical MDS Across Simple Families

Primary report:

- PDF: [`gmds_mds_comparison_report_2026-03-31.pdf`](/Users/pgajer/current_projects/grip/dev/design/pdf/gmds_mds_comparison_report_2026-03-31.pdf)
- LaTeX: [`gmds_mds_comparison_report_2026-03-31.tex`](/Users/pgajer/current_projects/grip/dev/design/pdf/gmds_mds_comparison_report_2026-03-31.tex)

Driver:

- [`tools/benchmark-gmds-mds-phase1.R`](/Users/pgajer/current_projects/grip/tools/benchmark-gmds-mds-phase1.R)

Main questions:

- Does GMDS reduce geodesic stress much more aggressively than classical MDS?
- On which simple graph families does that stress reduction still produce a geometrically reasonable embedding?
- Which bad pictures are due to the algorithm, and which are due to the setup of the benchmark itself?

Main findings:

- GMDS strongly reduced its own objective relative to classical MDS across the Phase 1 cases.
- On **flat orthogonal meshes**, single-path GMDS produced asymmetric or folded solutions even when stress was low.
- On **flat orthogonal meshes**, exact tie-averaging over all shortest paths largely fixed the symmetry-breaking artifact.
- On **ripple meshes in 2D**, the issue was largely a dimensionality obstruction: the surface cannot be faithfully represented in 2D, so poor-looking 2D embeddings are not a fair diagnostic.
- On **ripple and paraboloid meshes in 3D**, tie-averaged GMDS plus an MDS tether improved the comparison, but paraboloid singularities remained visually concerning.

This phase also clarified an important interpretation point:

- A low-GMDS-stress embedding is not necessarily visually faithful to the obvious Euclidean reference geometry.

## Phase 1 Follow-Up: Tie-Averaging and MDS Tether

Design note:

- [`geodesic_mds_hybrid_tether_and_tie_average_2026-03-31.md`](/Users/pgajer/current_projects/grip/dev/design/geodesic_mds_hybrid_tether_and_tie_average_2026-03-31.md)

This note answered the first pathology question on orthogonal grids:

`Was the flat-mesh artifact a fundamental failure of rectangular meshes, or an artifact of single-path tie breaking?`

Answer:

- On flat orthogonal meshes, the major artifact was not rectangular connectivity by itself.
- It was largely the deterministic selection of one shortest path out of a huge tied family.
- Exact tie-averaging fixed that problem on the flat mesh benchmark.

So after this point the remaining pathology story shifted away from flat grids and toward curved 3D surfaces, especially the paraboloid.

## Paraboloid Pathology Report

Primary report:

- PDF: [`paraboloid_gmds_pathology_report_2026-03-31.pdf`](/Users/pgajer/current_projects/grip/dev/design/pdf/paraboloid_gmds_pathology_report_2026-03-31.pdf)
- LaTeX: [`paraboloid_gmds_pathology_report_2026-03-31.tex`](/Users/pgajer/current_projects/grip/dev/design/pdf/paraboloid_gmds_pathology_report_2026-03-31.tex)

Driver:

- [`tools/benchmark-paraboloid-gmds-pathology.R`](/Users/pgajer/current_projects/grip/tools/benchmark-paraboloid-gmds-pathology.R)

Main questions:

- Is the paraboloid singularity mainly caused by the gradient-descent step schedule?
- Does the bad geometry appear gradually, or almost immediately?
- Can stronger constant tethering to MDS suppress the failure?

Key findings:

- The collapse happens very early, often within the first one or two GMDS steps away from the MDS initialization.
- Accepted steps become tiny quickly, yet the surface keeps sharpening into the bad basin.
- That points away from a simple overshoot story and toward an **objective-driven** pathology.
- A stronger constant MDS tether helps somewhat, but it does not eliminate the singularity.

This was the first strong piece of evidence that the missing ingredient is not just better line-search tuning.

## First-Order Regularization Report

Primary report:

- PDF: [`paraboloid_gmds_regularization_report_2026-03-31.pdf`](/Users/pgajer/current_projects/grip/dev/design/pdf/paraboloid_gmds_regularization_report_2026-03-31.pdf)
- LaTeX: [`paraboloid_gmds_regularization_report_2026-03-31.tex`](/Users/pgajer/current_projects/grip/dev/design/pdf/paraboloid_gmds_regularization_report_2026-03-31.tex)

Driver:

- [`tools/benchmark-paraboloid-gmds-regularization.R`](/Users/pgajer/current_projects/grip/tools/benchmark-paraboloid-gmds-regularization.R)

Main question:

`If the problem is missing regularization, is a simple local smoothness penalty enough?`

Key findings:

- Bending of the surface was only mildly reduced.
- The best smoothness-only and anchor-plus-smoothness runs still converged to the same singular-looking basin.
- Therefore the pathology was **not** just due to the absence of a generic first-order smoothness term.

This ruled out the simplest regularization hypothesis.

## Second-Order Bending Report

Primary report:

- PDF: [`paraboloid_gmds_bending_report_2026-03-31.pdf`](/Users/pgajer/current_projects/grip/dev/design/pdf/paraboloid_gmds_bending_report_2026-03-31.pdf)
- LaTeX: [`paraboloid_gmds_bending_report_2026-03-31.tex`](/Users/pgajer/current_projects/grip/dev/design/pdf/paraboloid_gmds_bending_report_2026-03-31.tex)

Driver:

- [`tools/benchmark-paraboloid-gmds-bending.R`](/Users/pgajer/current_projects/grip/tools/benchmark-paraboloid-gmds-bending.R)

Main question:

`If first-order smoothing is too weak, is an isotropic second-order rectangular-grid bending penalty enough?`

Key findings:

- Bending-only GMDS was almost indistinguishable from untethered GMDS on the `12 x 12` paraboloid.
- Anchor-plus-bending improved diagnostics only modestly.
- The same singular fold remained visible.
- The collapse still happened almost immediately after leaving the MDS initialization.

This is currently the strongest negative result:

- simple isotropic second-order bending regularization is **not enough**.

## Performance Study

The performance concern was addressed in the later section of the Phase 1 report and in the related action-plan note:

- [`gmds_performance_action_plan_2026-03-31.md`](/Users/pgajer/current_projects/grip/dev/design/gmds_performance_action_plan_2026-03-31.md)
- Phase 1 report performance section:
  - [`gmds_mds_comparison_report_2026-03-31.tex`](/Users/pgajer/current_projects/grip/dev/design/pdf/gmds_mds_comparison_report_2026-03-31.tex)
  - [`gmds_mds_comparison_report_2026-03-31.pdf`](/Users/pgajer/current_projects/grip/dev/design/pdf/gmds_mds_comparison_report_2026-03-31.pdf)

The main point was:

- GMDS is a compiled C++ optimizer in the current implementation, not an R-only optimizer.
- Even so, it remains much slower than classical MDS because it repeatedly evaluates and differentiates a nonconvex all-pairs path objective.

The performance pass improved cache construction and optimization speed materially, but it did not erase the fundamental algorithmic gap to classical MDS.

## What the Pathologies Seem to Mean

At this point, the pathologies appear to come from different sources in different examples.

### Flat orthogonal meshes

The dominant issue was:

- deterministic single-path selection on huge tied shortest-path families.

Tie-averaging largely fixed this.

### Ripple in 2D

The dominant issue was:

- dimensionality mismatch.

A curved ripple surface should not be judged by a 2D fidelity standard.

### Paraboloid in 3D

The remaining pathology appears to come from the **incompleteness of the pure geodesic objective**.

GMDS is very good at matching graph-geodesic path lengths, but path-length fidelity alone does not forbid:

- fold-over,
- local cell inversion,
- self-approach,
- concentrated curvature,
- or aesthetically poor but geodesically acceptable singular structures.

The experiments so far strongly suggest:

- the paraboloid pathology is **not primarily a line-search tuning issue**,
- **not primarily a single-path tie-breaking issue**,
- and **not removable by simple first-order or isotropic second-order regularization**.

## Current Best Answer to the Main Question

`Where are these GMDS pathologies coming from?`

Current answer:

1. On simple orthogonal meshes, some of the earliest failures came from an artificial source: deterministic shortest-path tie breaking. Tie-averaging solved that part.
2. On some families, especially ripple in 2D, the problem is not GMDS at all but a mismatch between the target geometry and the embedding dimension.
3. On curved 3D meshes such as the paraboloid, the remaining pathology seems to come from the objective itself: geodesic path stress does not contain enough geometric information to rule out fold-over and singular concentration of curvature.

That is why the next likely direction is not just more step-size tuning, and not just more isotropic smoothing. The next experiments should probably use a regularizer or constraint that reacts directly to fold-over, cell orientation, self-intersection, or topology preservation.

## Most Relevant Files at a Glance

### Algorithm and implementation

- base GMDS implementation:
  - [`R/grip_quality.R`](/Users/pgajer/current_projects/grip/R/grip_quality.R)
- compiled optimizer:
  - [`src/geodesic_mds_rcpp.cpp`](/Users/pgajer/current_projects/grip/src/geodesic_mds_rcpp.cpp)
- bending optimizer:
  - [`src/geodesic_mds_bending_rcpp.cpp`](/Users/pgajer/current_projects/grip/src/geodesic_mds_bending_rcpp.cpp)
- bending R overlay:
  - [`R/zz_geodesic_mds_bending.R`](/Users/pgajer/current_projects/grip/R/zz_geodesic_mds_bending.R)
  - [`R/zz_geodesic_mds_bending_core.R`](/Users/pgajer/current_projects/grip/R/zz_geodesic_mds_bending_core.R)
  - [`R/zz_geodesic_mds_bending_score.R`](/Users/pgajer/current_projects/grip/R/zz_geodesic_mds_bending_score.R)
  - [`R/zz_geodesic_mds_bending_opt_r.R`](/Users/pgajer/current_projects/grip/R/zz_geodesic_mds_bending_opt_r.R)
  - [`R/zz_geodesic_mds_bending_optimize.R`](/Users/pgajer/current_projects/grip/R/zz_geodesic_mds_bending_optimize.R)
  - [`R/zz_geodesic_mds_bending_stats.R`](/Users/pgajer/current_projects/grip/R/zz_geodesic_mds_bending_stats.R)

### Tests

- [`tests/testthat/test-geodesic-mds.R`](/Users/pgajer/current_projects/grip/tests/testthat/test-geodesic-mds.R)
- [`tests/testthat/test-geodesic-mds-hybrid.R`](/Users/pgajer/current_projects/grip/tests/testthat/test-geodesic-mds-hybrid.R)
- [`tests/testthat/test-geodesic-mds-smoothness.R`](/Users/pgajer/current_projects/grip/tests/testthat/test-geodesic-mds-smoothness.R)
- [`tests/testthat/test-geodesic-mds-bending.R`](/Users/pgajer/current_projects/grip/tests/testthat/test-geodesic-mds-bending.R)
- [`tests/testthat/test-gmds-mds-phase1.R`](/Users/pgajer/current_projects/grip/tests/testthat/test-gmds-mds-phase1.R)

### Report generators

- [`tools/benchmark-gmds-mds-phase1.R`](/Users/pgajer/current_projects/grip/tools/benchmark-gmds-mds-phase1.R)
- [`tools/benchmark-paraboloid-gmds-pathology.R`](/Users/pgajer/current_projects/grip/tools/benchmark-paraboloid-gmds-pathology.R)
- [`tools/benchmark-paraboloid-gmds-regularization.R`](/Users/pgajer/current_projects/grip/tools/benchmark-paraboloid-gmds-regularization.R)
- [`tools/benchmark-paraboloid-gmds-bending.R`](/Users/pgajer/current_projects/grip/tools/benchmark-paraboloid-gmds-bending.R)

