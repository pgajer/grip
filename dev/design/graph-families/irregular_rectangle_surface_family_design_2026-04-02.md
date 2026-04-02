# Irregular Rectangle Surface Family Design

Date: 2026-04-02

Status: Planned only. Not started.

## Purpose

This note proposes a new synthetic family for `grip`:

- irregular rectangles,
- with squares as the special case `h = w`,
- preserving ordinary rectangular-mesh topology,
- but breaking the strong symmetry of the regular lattice.

The motivation is the recent MISF-GMDS paraboloid work. We currently have:

- regular rectangular meshes via
  [mesh.surface.graph()](/Users/pgajer/current_projects/grip/R/graph_helpers.R#L3702),
- occupied/perforated rectangular meshes via
  [occupied.mesh.surface.graph()](/Users/pgajer/current_projects/grip/R/graph_helpers.R#L3095),
- and several genuinely irregular non-rectangular families such as
  [irregular.annulus.surface.graph()](/Users/pgajer/current_projects/grip/R/graph_helpers.R#L5796)
  and
  [irregular.sphere.surface.graph()](/Users/pgajer/current_projects/grip/R/graph_helpers.R#L5166).

What is still missing is a family that is:

- topologically just a rectangle,
- simply connected,
- boundary-preserving,
- but no longer a perfectly regular lattice.

That is exactly the right control family for MISF-GMDS. It lets us ask:

- does the bad basin mainly come from lattice symmetry,
- or from the GMDS objective itself?

Unlike occupied/perforated meshes, this family would not change topology.
Unlike annulus/sphere families, it would stay very close to the rectangular
mesh setting that exposed the basin problem.

## Core Design Decision

The new family should **keep the graph adjacency equal to the ordinary
rectangular mesh**, and it should introduce irregularity only through the
canonical geometry used to assign coordinates and edge weights.

That means:

- same `edges.mesh(h, w, connectivity = ...)`,
- same vertex count `n = h * w`,
- same rectangular boundary combinatorics,
- but nonuniform parameter coordinates and optional smooth interior warp.

This is important because it isolates the effect of symmetry breaking without
mixing in holes, slits, or nontrivial topology.

## Proposed Public API

I recommend the following exported helpers.

### 1. Canonical planar parameterization

- `irregular.rectangle.param.coords(h, w = h, x_scale = 1, y_scale = 1, row_irregularity = 0.20, col_irregularity = 0.20, row_phase = 0.35, col_phase = 0.65, interior_warp = 0.08, shear = 0, min_step_ratio = 0.30)`

This would return an `n x 2` numeric matrix with columns `u` and `v`.

### 2. 3D embedding wrapper

- `irregular.rectangle.surface.embedding(h, w = h, surface = c("flat", "saddle", "paraboloid", "ripple"), amplitude = 0.75, freq_u = 1, freq_v = 1, x_scale = 1, y_scale = 1, row_irregularity = 0.20, col_irregularity = 0.20, row_phase = 0.35, col_phase = 0.65, interior_warp = 0.08, shear = 0, min_step_ratio = 0.30)`

This would return an `n x 3` numeric matrix with columns `x`, `y`, and `z`.

### 3. Weighted graph bundle

- `irregular.rectangle.surface.graph(h, w = h, surface = c("flat", "saddle", "paraboloid", "ripple"), amplitude = 0.75, freq_u = 1, freq_v = 1, x_scale = 1, y_scale = 1, row_irregularity = 0.20, col_irregularity = 0.20, row_phase = 0.35, col_phase = 0.65, interior_warp = 0.08, shear = 0, min_step_ratio = 0.30, connectivity = c("orthogonal", "diagonal"), normalize = c("median", "mean", "none"))`

This would return the usual weighted-graph bundle:

- `edges`
- `n`
- `edge_weights`
- `coords_surface`
- `coords_param`
- `weight_scale`
- `family = "irregular.rectangle"`
- `surface`
- `connectivity`
- irregularity metadata
- `label`

I would not add separate square-specific functions. A square is simply
`irregular.rectangle.*(h = s, w = s, ...)`.

## Exact Construction Rules

The family should be deterministic and geometry-first.

### Step 1. Start from the regular rectangular parameter grid

Let

- `u0_j`, `j = 1, ..., w`,
- `v0_i`, `i = 1, ..., h`

be the regular 1D grid coordinates produced by the same centered axis logic as
the current mesh helpers.

Equivalently:

`u0_j in [-x_scale, x_scale]`, monotone increasing  
`v0_i in [-y_scale, y_scale]`, monotone increasing

with rectangular indexing exactly matching
[mesh.surface.graph()](/Users/pgajer/current_projects/grip/R/graph_helpers.R#L3702).

### Step 2. Perturb interval lengths, not vertex order

This is the key design choice.

Instead of perturbing vertex coordinates directly, define positive horizontal
and vertical interval weights:

For columns `j = 1, ..., w - 1`, let

`s_j = (j - 0.5) / (w - 1)`

and define

`dx_j = 1 + col_irregularity * [0.65 cos(2 pi s_j + col_phase) + 0.35 sin(4 pi s_j + 0.5 col_phase)]`

For rows `i = 1, ..., h - 1`, let

`t_i = (i - 0.5) / (h - 1)`

and define

`dy_i = 1 + row_irregularity * [0.65 sin(2 pi t_i + row_phase) + 0.35 cos(4 pi t_i + 0.5 row_phase)]`

Then clip:

`dx_j <- max(dx_j, min_step_ratio)`  
`dy_i <- max(dy_i, min_step_ratio)`

This guarantees:

- strict monotonicity,
- no cell inversion at the parameter level,
- deterministic irregular spacing,
- boundary order preserved.

### Step 3. Renormalize cumulative sums so the boundary stays rectangular

Define cumulative coordinates:

`u_1 = -x_scale`  
`u_{j+1} = -x_scale + 2 x_scale * sum_{m <= j} dx_m / sum_m dx_m`

and similarly

`v_1 = -y_scale`  
`v_{i+1} = -y_scale + 2 y_scale * sum_{m <= i} dy_m / sum_m dy_m`

This fixes the outer rectangle exactly:

- left boundary at `-x_scale`,
- right boundary at `x_scale`,
- bottom boundary at `-y_scale`,
- top boundary at `y_scale`.

So the irregularity changes internal spacing, not the existence or location of
the rectangle boundary as a whole.

### Step 4. Optional smooth interior warp that vanishes on the boundary

To avoid the family being merely “nonuniform row spacing,” add a small
boundary-vanishing interior warp.

Let regular normalized coordinates be

`ubar_ij = u0_j / x_scale`  
`vbar_ij = v0_i / y_scale`

Define the envelope

`B_ij = (1 - ubar_ij^2) (1 - vbar_ij^2)`

and then warp:

`u'_ij = u_j + interior_warp * x_scale * B_ij * sin(pi vbar_ij) * cos(pi ubar_ij + col_phase)`

`v'_ij = v_i + interior_warp * y_scale * B_ij * sin(pi ubar_ij) * cos(pi vbar_ij + row_phase)`

Finally, optionally apply a small shear:

`u''_ij = u'_ij + shear * v'_ij`  
`v''_ij = v'_ij`

The envelope `B_ij` forces the warp to vanish on the boundary, so:

- corners stay corners,
- boundary vertices stay on the four boundary arcs,
- topology stays rectangular.

### Step 5. Lift into 3D exactly like the mesh family

After constructing `coords_param = (u'', v'')`, define

- `flat`: `z = 0`
- `saddle`: `z = amplitude * (u^2 - v^2)`
- `paraboloid`: `z = amplitude * (u^2 + v^2)`
- `ripple`: `z = amplitude * sin(pi * freq_u * u) * cos(pi * freq_v * v)`

So this family stays aligned with the existing
[mesh.surface.embedding()](/Users/pgajer/current_projects/grip/R/graph_helpers.R#L3665)
logic, except that it adds `flat` and replaces regular parameter coordinates by
deterministically irregular ones.

### Step 6. Use the induced Euclidean edge lengths as weights

Use the existing edge-weight helper:

`edge_weights = .edge.weights.from.embedding(edges, coords_surface, normalize = ...)`

and keep the bundle format parallel to
[mesh.surface.graph()](/Users/pgajer/current_projects/grip/R/graph_helpers.R#L3702).

## Why These Parameters

I recommend these defaults:

- `row_irregularity = 0.20`
- `col_irregularity = 0.20`
- `row_phase = 0.35`
- `col_phase = 0.65`
- `interior_warp = 0.08`
- `shear = 0`
- `min_step_ratio = 0.30`

Reasoning:

- `0.20` spacing irregularity is large enough to visibly break symmetry.
- `0.08` interior warp is enough to avoid a purely separable row/column grid.
- `min_step_ratio = 0.30` prevents near-collapsed cells.
- distinct row/column phases avoid accidental bilateral symmetry.

These values should produce a family that is clearly irregular but still
recognizably rectangular.

## Expected Fields In The Graph Bundle

In addition to the standard fields, I recommend storing:

- `coords_regular_param`
- `row_breaks`
- `col_breaks`
- `row_irregularity`
- `col_irregularity`
- `row_phase`
- `col_phase`
- `interior_warp`
- `shear`
- `min_step_ratio`

This will help later diagnostics and report generation.

## Recommended File-By-File Implementation

### 1. Core helper internals

File:
- [graph_helpers.R](/Users/pgajer/current_projects/grip/R/graph_helpers.R)

Add private helpers:

- `.irregular.rectangle.axis.breakpoints()`
- `.irregular.rectangle.param.coords()`

These should:

- validate dimensions and irregularity parameters,
- build monotone perturbed axis breakpoints,
- apply the boundary-vanishing interior warp,
- return deterministic parameter coordinates.

### 2. Exported public helpers

File:
- [graph_helpers.R](/Users/pgajer/current_projects/grip/R/graph_helpers.R)

Add public functions:

- `irregular.rectangle.param.coords()`
- `irregular.rectangle.surface.embedding()`
- `irregular.rectangle.surface.graph()`

These should be documented under a new helper block, parallel to:

- [mesh_surface_helpers](/Users/pgajer/current_projects/grip/R/graph_helpers.R#L3616)
- [irregular_annulus_surface_helpers](/Users/pgajer/current_projects/grip/R/graph_helpers.R#L5708)

### 3. Documentation

Files to regenerate:

- [NAMESPACE](/Users/pgajer/current_projects/grip/NAMESPACE)
- new `.Rd` files under
  [man](/Users/pgajer/current_projects/grip/man)

Update the catalog note:

- [graph_families_generated_in_thread_2026-03-31.md](/Users/pgajer/current_projects/grip/dev/design/graph_families_generated_in_thread_2026-03-31.md)

### 4. Tests

Files:

- [test-graph-helpers.R](/Users/pgajer/current_projects/grip/tests/testthat/test-graph-helpers.R)
- or a new
  [test-irregular-rectangle-helpers.R](/Users/pgajer/current_projects/grip/tests/testthat/test-irregular-rectangle-helpers.R)

Add tests for:

- deterministic output for repeated calls,
- strict monotonicity of row and column coordinates,
- boundary preservation,
- same edge set as `edges.mesh(h, w, connectivity = ...)`,
- positive finite edge weights,
- `flat` family giving zero `z`,
- square special case `h = w`.

### 5. Benchmark integration

Files to extend after the family exists:

- [benchmark-gmds-misf-paraboloid.R](/Users/pgajer/current_projects/grip/tools/benchmark-gmds-misf-paraboloid.R)
- future MISF-GMDS reports

The first benchmark use should compare:

- regular flat rectangle,
- irregular flat rectangle,
- regular paraboloid rectangle,
- irregular paraboloid rectangle

for sizes such as:

- `12 x 12`
- `15 x 15`
- `20 x 20`

## Recommended Experimental Use

This family should be used for two separate questions.

### Question 1. Does irregularity alone improve basin selection?

Compare:

- `mesh.surface.graph(..., surface = "paraboloid")`
- `irregular.rectangle.surface.graph(..., surface = "paraboloid")`

with the same mesh adjacency and nearly the same overall domain.

This isolates symmetry breaking from topology change.

### Question 2. Does MISF help more on regular or irregular rectangles?

Compare:

- `cmdscale -> GMDS`
- `MISF-GMDS`

on both regular and irregular rectangles.

This is the cleanest test of whether MISF is compensating for a bad global
spectral basin or whether it is adding value even after the lattice symmetry is
already broken.

## Why This Family Is Better Than Using Occupied Meshes For This Question

Occupied/perforated meshes are still valuable, but they change too many things
at once:

- topology,
- boundary shape,
- local valence,
- path families,
- and symmetry.

The irregular-rectangle family changes much less:

- same topology,
- same adjacency pattern,
- same rectangular boundary class,
- but lower symmetry.

That makes it a much better diagnostic family for MISF-GMDS.

## Proposed Execution Order

1. Implement the core parameter-coordinate helper.
2. Add the public `*.surface.embedding` and `*.surface.graph` wrappers.
3. Add unit tests and documentation.
4. Add one small gallery/inspection script for `flat` and `paraboloid`.
5. Add the family to the next MISF-GMDS comparison report.

## Recommendation

If we add only one new 2D-ish family next for MISF-GMDS, it should be this one.

It is the cleanest missing bridge between:

- regular rectangular meshes,
- occupied/perforated rectangular meshes,
- and the already implemented irregular annulus/sphere families.
