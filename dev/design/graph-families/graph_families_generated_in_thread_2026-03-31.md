# Historical Graph-Family Design Catalog

This catalog records the March 2026 graph-family designs implemented in
[R/graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R).
The tables include both public constructors and internal helpers. The API-status
column identifies exports in grip 0.2.0; internal helpers are not callable with
`grip::` and are not a supported public interface. Signatures below are historical
design summaries, not current usage documentation. Consult the installed help
for current arguments, or use the public weighted-graph constructors, whose
results include coordinates as well as edges and weights.

Scope note:
- This note focuses on family constructors, weighted-geometry wrappers, mask constructors, and benchmark-oriented helpers.
- Older topology-only generators that already existed before the thread, such as `edges.mesh()`, `edges.torus()`, `edges.sphere()`, `edges.cylinder()`, and `edges.kary.tree()`, are mentioned only when a new family wraps or extends them.

## Common conventions

- `*.surface.embedding(...)` returns coordinates only.
- `*.surface.graph(...)` returns a reusable weighted-graph bundle, typically including `edges`, `edge_weights`, `coords_surface`, `coords_param` or equivalent canonical coordinates, and family metadata.
- `*.solid.embedding(...)` / `*.solid.graph(...)` follow the same pattern for volumetric families.
- `normalize` controls post-hoc edge-weight normalization:
  - `"median"`: rescale so median edge weight is `1`.
  - `"mean"`: rescale so mean edge weight is `1`.
  - `"none"`: keep raw induced lengths.
- For all recursive families, `level` is the recursion depth.

## 1. Curved lifts of regular lattice families

### 1.1 Mesh surface family

| Function | Signature | Implementation | API status (0.2.0) |
|---|---|---|---|
| `mesh.surface.embedding` | `mesh.surface.embedding(h, w = h, surface = c("saddle", "paraboloid", "ripple"), amplitude = 0.75, freq_u = 1, freq_v = 1, x_scale = 1, y_scale = 1)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `mesh.surface.graph` | `mesh.surface.graph(h, w = h, surface = c("saddle", "paraboloid", "ripple"), amplitude = 0.75, freq_u = 1, freq_v = 1, x_scale = 1, y_scale = 1, normalize = c("median", "mean", "none"))` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |

Parameters:
- `h`, `w`: mesh height and width in grid cells.
- `surface`: lift profile; `saddle` uses a signed quadratic, `paraboloid` uses a radial quadratic, and `ripple` uses oscillatory height variation.
- `amplitude`: deformation strength in the lifted coordinate.
- `freq_u`, `freq_v`: oscillation frequencies along the two grid directions; mainly relevant for `ripple`.
- `x_scale`, `y_scale`: anisotropic scaling of the parameter grid before lifting.
- `normalize`: edge-weight normalization mode for the weighted-graph wrapper.

### 1.1a Irregular rectangle surface family

| Function | Signature | Implementation | API status (0.2.0) |
|---|---|---|---|
| `irregular.rectangle.param.coords` | `irregular.rectangle.param.coords(h, w = h, x_scale = 1, y_scale = 1, row_irregularity = 0.20, col_irregularity = 0.20, row_phase = 0.35, col_phase = 0.65, interior_warp = 0.08, shear = 0, min_step_ratio = 0.30)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `irregular.rectangle.surface.embedding` | `irregular.rectangle.surface.embedding(h, w = h, surface = c("flat", "saddle", "paraboloid", "ripple"), amplitude = 0.75, freq_u = 1, freq_v = 1, x_scale = 1, y_scale = 1, row_irregularity = 0.20, col_irregularity = 0.20, row_phase = 0.35, col_phase = 0.65, interior_warp = 0.08, shear = 0, min_step_ratio = 0.30)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `irregular.rectangle.surface.graph` | `irregular.rectangle.surface.graph(h, w = h, surface = c("flat", "saddle", "paraboloid", "ripple"), amplitude = 0.75, freq_u = 1, freq_v = 1, x_scale = 1, y_scale = 1, row_irregularity = 0.20, col_irregularity = 0.20, row_phase = 0.35, col_phase = 0.65, interior_warp = 0.08, shear = 0, min_step_ratio = 0.30, normalize = c("median", "mean", "none"))` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |

Parameters:
- `h`, `w`: rectangle height and width in grid cells; squares are the special case `h = w`.
- `row_irregularity`, `col_irregularity`: deterministic irregularity levels for row and column spacing.
- `row_phase`, `col_phase`: phase offsets used in the spacing perturbations.
- `interior_warp`: boundary-vanishing interior warp strength.
- `shear`: optional affine shear applied after the interior warp.
- `min_step_ratio`: lower bound on perturbed row/column interval lengths.
- `surface`, `amplitude`, `freq_u`, `freq_v`, `x_scale`, `y_scale`, `normalize`: same role as in the regular mesh family, except that `surface` also allows `flat`.

### 1.2 Cylinder surface family

| Function | Signature | Implementation | API status (0.2.0) |
|---|---|---|---|
| `cylinder.surface.embedding` | `cylinder.surface.embedding(h, w = h, surface = c("barrel", "hourglass", "wavy"), radius = 1, height = 2, amplitude = 0.3, freq_theta = 2, freq_z = 1, twist = 0.25)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `cylinder.surface.graph` | `cylinder.surface.graph(h, w = h, surface = c("barrel", "hourglass", "wavy"), radius = 1, height = 2, amplitude = 0.3, freq_theta = 2, freq_z = 1, twist = 0.25, normalize = c("median", "mean", "none"))` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |

Parameters:
- `h`, `w`: counts along the axial and wrapped directions.
- `surface`: radial profile; `barrel` bulges outward, `hourglass` narrows in the middle, and `wavy` adds oscillatory variation.
- `radius`, `height`: base cylinder dimensions.
- `amplitude`: strength of radial modulation.
- `freq_theta`, `freq_z`: frequencies around the cylinder and along its axis.
- `twist`: helical phase drift along the axis.
- `normalize`: edge-weight normalization mode.

### 1.3 Torus surface family

| Function | Signature | Implementation | API status (0.2.0) |
|---|---|---|---|
| `torus.surface.embedding` | `torus.surface.embedding(h, w = h, surface = c("standard", "pinched", "wavy"), major_radius = 2, minor_radius = 0.75, amplitude = 0.2, freq_major = 2, freq_minor = 1, twist = 0.25)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `torus.surface.graph` | `torus.surface.graph(h, w = h, surface = c("standard", "pinched", "wavy"), major_radius = 2, minor_radius = 0.75, amplitude = 0.2, freq_major = 2, freq_minor = 1, twist = 0.25, normalize = c("median", "mean", "none"))` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |

Parameters:
- `h`, `w`: counts along the major and minor torus cycles.
- `surface`: torus deformation mode; `standard` is unmodified, `pinched` constricts the tube, `wavy` oscillates the tube radius.
- `major_radius`, `minor_radius`: torus radii.
- `amplitude`: strength of radius modulation.
- `freq_major`, `freq_minor`: oscillation counts along the two torus directions.
- `twist`: phase twist between the major and minor cycles.
- `normalize`: edge-weight normalization mode.

### 1.4 Sphere surface family

| Function | Signature | Implementation | API status (0.2.0) |
|---|---|---|---|
| `sphere.surface.embedding` | `sphere.surface.embedding(h, w = h, surface = c("standard", "ellipsoid", "wavy"), radius = 1, amplitude = 0.2, freq_theta = 3, freq_lat = 2, twist = 0.25)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `sphere.surface.graph` | `sphere.surface.graph(h, w = h, surface = c("standard", "ellipsoid", "wavy"), radius = 1, amplitude = 0.2, freq_theta = 3, freq_lat = 2, twist = 0.25, normalize = c("median", "mean", "none"))` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |

Parameters:
- `h`, `w`: numbers of latitude bands and samples per band in the regular spherical family.
- `surface`: spherical deformation mode; `standard` is round, `ellipsoid` introduces axis distortion, `wavy` perturbs the radius.
- `radius`: base radius.
- `amplitude`: deformation strength.
- `freq_theta`, `freq_lat`: oscillation counts along longitude and latitude.
- `twist`: phase drift by latitude.
- `normalize`: edge-weight normalization mode.

## 2. Recursive square-mask families

### 2.1 Mask constructors

| Function | Signature | Implementation | API status (0.2.0) |
|---|---|---|---|
| `mask.cross` | `mask.cross(k = 5, arm_width = 1)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |
| `mask.border` | `mask.border(k = 5, thickness = 1)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |
| `mask.corner` | `mask.corner(k = 5, width = 2, corner = c("top_left", "top_right", "bottom_left", "bottom_right"))` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |
| `mask.asymmetric.holes` | `mask.asymmetric.holes(k = 5, hole_size = 1)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |

Parameters:
- `k`: side length of the square recursive mask.
- `arm_width`: thickness of the retained cross arms.
- `thickness`: retained border thickness.
- `width`: retained size of the selected corner block.
- `corner`: which corner to keep for `mask.corner`.
- `hole_size`: size of the asymmetric interior holes.

### 2.2 Generic recursive square-mask graph family

| Function | Signature | Implementation | API status (0.2.0) |
|---|---|---|---|
| `edges.recursive.mask.grid` | `edges.recursive.mask.grid(mask, level = 2)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `recursive.mask.grid.surface.embedding` | `recursive.mask.grid.surface.embedding(mask, level = 2, surface = c("saddle", "paraboloid", "ripple"), amplitude = 0.75, freq_u = 1, freq_v = 1, x_scale = 1, y_scale = 1)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `recursive.mask.grid.surface.graph` | `recursive.mask.grid.surface.graph(mask, level = 2, surface = c("saddle", "paraboloid", "ripple"), amplitude = 0.75, freq_u = 1, freq_v = 1, x_scale = 1, y_scale = 1, normalize = c("median", "mean", "none"))` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |

Parameters:
- `mask`: logical square keep-mask used at each recursion step.
- `level`: recursion depth.
- `surface`, `amplitude`, `freq_u`, `freq_v`, `x_scale`, `y_scale`: same meaning as in the mesh surface family, applied to the occupied recursive cells.
- `normalize`: edge-weight normalization mode.

### 2.3 Named square-mask wrappers

#### Sierpinski carpet

| Function | Signature | Implementation | API status (0.2.0) |
|---|---|---|---|
| `edges.sierpinski.carpet` | `edges.sierpinski.carpet(level = 2)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |
| `sierpinski.carpet.surface.embedding` | `sierpinski.carpet.surface.embedding(level = 2, surface = c("saddle", "paraboloid", "ripple"), amplitude = 0.75, freq_u = 1, freq_v = 1, x_scale = 1, y_scale = 1)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `sierpinski.carpet.surface.graph` | `sierpinski.carpet.surface.graph(level = 2, surface = c("saddle", "paraboloid", "ripple"), amplitude = 0.75, freq_u = 1, freq_v = 1, x_scale = 1, y_scale = 1, normalize = c("median", "mean", "none"))` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |

Parameters:
- `level`: carpet recursion depth.
- Remaining geometry parameters match `recursive.mask.grid.surface.*`.

#### Vicsek fractal

| Function | Signature | Implementation | API status (0.2.0) |
|---|---|---|---|
| `edges.vicsek` | `edges.vicsek(level = 2)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `vicsek.surface.embedding` | `vicsek.surface.embedding(level = 2, surface = c("saddle", "paraboloid", "ripple"), amplitude = 0.75, freq_u = 1, freq_v = 1, x_scale = 1, y_scale = 1)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `vicsek.surface.graph` | `vicsek.surface.graph(level = 2, surface = c("saddle", "paraboloid", "ripple"), amplitude = 0.75, freq_u = 1, freq_v = 1, x_scale = 1, y_scale = 1, normalize = c("median", "mean", "none"))` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |

Parameters:
- `level`: Vicsek recursion depth.
- Remaining geometry parameters match `recursive.mask.grid.surface.*`.

## 3. Deterministic perforated meshes

### 3.1 Occupancy-mask helpers

| Function | Signature | Implementation | API status (0.2.0) |
|---|---|---|---|
| `keep.periodic.holes` | `keep.periodic.holes(h, w = h, hole_period = 4, hole_height = 1, hole_width = hole_height, row_offset = 2, col_offset = 2)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |
| `keep.staggered.windows` | `keep.staggered.windows(h, w = h, window_height = 1, window_width = 2, row_period = 4, col_period = 5, row_offset = 2, col_offset = 2)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |
| `keep.slit.channels` | `keep.slit.channels(h, w = h, orientation = c("vertical", "horizontal"), slit_period = 5, slit_width = 1, bridge_spacing = 4, bridge_size = 1, offset = 2)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |
| `keep.asymmetric.notches` | `keep.asymmetric.notches(h, w = h, notch_depth = 3, notch_width = 2)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |

Parameters:
- `h`, `w`: occupancy-grid dimensions.
- `hole_period`, `row_period`, `col_period`, `slit_period`: periodic spacing between removed regions.
- `hole_height`, `hole_width`, `window_height`, `window_width`, `slit_width`, `bridge_size`, `notch_depth`, `notch_width`: sizes of removed or retained motif pieces.
- `row_offset`, `col_offset`, `offset`: phase offsets of the repeating motif.
- `orientation`: whether slits run vertically or horizontally.
- `bridge_spacing`: spacing between retained bridges across slit channels.

### 3.2 Finite occupied-mesh family

| Function | Signature | Implementation | API status (0.2.0) |
|---|---|---|---|
| `edges.occupied.mesh` | `edges.occupied.mesh(keep)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `occupied.mesh.surface.embedding` | `occupied.mesh.surface.embedding(keep, surface = c("saddle", "paraboloid", "ripple"), amplitude = 0.75, freq_u = 1, freq_v = 1, x_scale = 1, y_scale = 1)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `occupied.mesh.surface.graph` | `occupied.mesh.surface.graph(keep, surface = c("saddle", "paraboloid", "ripple"), amplitude = 0.75, freq_u = 1, freq_v = 1, x_scale = 1, y_scale = 1, normalize = c("median", "mean", "none"))` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |

Parameters:
- `keep`: logical occupancy matrix; `TRUE` cells are retained as mesh vertices/cells.
- `surface`, `amplitude`, `freq_u`, `freq_v`, `x_scale`, `y_scale`: same lifted-geometry controls as in the mesh family.
- `normalize`: edge-weight normalization mode.

## 4. Recursive triangular gasket families

### 4.1 Triangle mask constructors

| Function | Signature | Implementation | API status (0.2.0) |
|---|---|---|---|
| `mask.triangle.classic` | `mask.triangle.classic()` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |
| `mask.triangle.bridge` | `mask.triangle.bridge(missing = c("top", "left", "right"))` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |

Parameters:
- `missing`: which corner subtriangle is removed in the bridge variant.

### 4.2 Generic recursive triangle-mask family

| Function | Signature | Implementation | API status (0.2.0) |
|---|---|---|---|
| `edges.recursive.triangle.mask` | `edges.recursive.triangle.mask(mask = mask.triangle.classic(), level = 2)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `recursive.triangle.mask.surface.embedding` | `recursive.triangle.mask.surface.embedding(mask = mask.triangle.classic(), level = 2, surface = c("flat", "saddle", "paraboloid", "ripple", "folded"), amplitude = 0.75, freq_u = 1, freq_v = 1, x_scale = 1, y_scale = 1)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `recursive.triangle.mask.surface.graph` | `recursive.triangle.mask.surface.graph(mask = mask.triangle.classic(), level = 2, surface = c("flat", "saddle", "paraboloid", "ripple", "folded"), amplitude = 0.75, freq_u = 1, freq_v = 1, x_scale = 1, y_scale = 1, normalize = c("median", "mean", "none"))` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |

Parameters:
- `mask`: length-4 keep-mask for the recursive subtriangles.
- `level`: recursion depth.
- `surface`: lifted-triangle geometry; `flat` leaves the triangle planar, `folded` introduces a crease-like lift, and the remaining options are analogous to the square-grid lifts.
- `amplitude`, `freq_u`, `freq_v`, `x_scale`, `y_scale`: deformation controls in triangle parameter coordinates.
- `normalize`: edge-weight normalization mode.

### 4.3 Named Sierpinski triangle wrapper

| Function | Signature | Implementation | API status (0.2.0) |
|---|---|---|---|
| `edges.sierpinski.triangle` | `edges.sierpinski.triangle(level = 2)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |
| `sierpinski.triangle.surface.embedding` | `sierpinski.triangle.surface.embedding(level = 2, surface = c("flat", "saddle", "paraboloid", "ripple", "folded"), amplitude = 0.75, freq_u = 1, freq_v = 1, x_scale = 1, y_scale = 1)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `sierpinski.triangle.surface.graph` | `sierpinski.triangle.surface.graph(level = 2, surface = c("flat", "saddle", "paraboloid", "ripple", "folded"), amplitude = 0.75, freq_u = 1, freq_v = 1, x_scale = 1, y_scale = 1, normalize = c("median", "mean", "none"))` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |

Parameters:
- `level`: recursion depth for the classic Sierpinski triangle.
- Remaining geometry parameters match `recursive.triangle.mask.surface.*`.

## 5. Recursive tetrahedral gasket families

### 5.1 Tetrahedron mask constructors

| Function | Signature | Implementation | API status (0.2.0) |
|---|---|---|---|
| `mask.tetrahedron.classic` | `mask.tetrahedron.classic()` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |
| `mask.tetrahedron.corner.missing` | `mask.tetrahedron.corner.missing(omit = c("apex", "base_left", "base_right", "base_back"))` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |

Parameters:
- `omit`: which corner subtetrahedron to remove from the retained corner set.

### 5.2 Generic recursive tetrahedron-mask family

| Function | Signature | Implementation | API status (0.2.0) |
|---|---|---|---|
| `edges.recursive.tetrahedron.mask` | `edges.recursive.tetrahedron.mask(mask = mask.tetrahedron.classic(), level = 2)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `recursive.tetrahedron.mask.surface.embedding` | `recursive.tetrahedron.mask.surface.embedding(mask = mask.tetrahedron.classic(), level = 2, surface = c("standard", "squashed", "twisted", "wavy"), amplitude = 0.3, freq = 2, twist = 0.6)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `recursive.tetrahedron.mask.surface.graph` | `recursive.tetrahedron.mask.surface.graph(mask = mask.tetrahedron.classic(), level = 2, surface = c("standard", "squashed", "twisted", "wavy"), amplitude = 0.3, freq = 2, twist = 0.6, normalize = c("median", "mean", "none"))` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |

Parameters:
- `mask`: length-4 keep-mask over the corner subtetrahedra.
- `level`: recursion depth.
- `surface`: tetrahedral deformation mode; `squashed` changes aspect, `twisted` rotates layers, `wavy` modulates coordinates.
- `amplitude`: deformation strength.
- `freq`: oscillation count used by the `wavy` profile.
- `twist`: twist intensity for twisted/twisted-like variants.
- `normalize`: edge-weight normalization mode.

### 5.3 Named Sierpinski tetrahedron wrapper

| Function | Signature | Implementation | API status (0.2.0) |
|---|---|---|---|
| `edges.sierpinski.tetrahedron` | `edges.sierpinski.tetrahedron(level = 2)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |
| `sierpinski.tetrahedron.surface.embedding` | `sierpinski.tetrahedron.surface.embedding(level = 2, surface = c("standard", "squashed", "twisted", "wavy"), amplitude = 0.3, freq = 2, twist = 0.6)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `sierpinski.tetrahedron.surface.graph` | `sierpinski.tetrahedron.surface.graph(level = 2, surface = c("standard", "squashed", "twisted", "wavy"), amplitude = 0.3, freq = 2, twist = 0.6, normalize = c("median", "mean", "none"))` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |

Parameters:
- `level`: recursion depth for the classic Sierpinski tetrahedron.
- Remaining geometry parameters match `recursive.tetrahedron.mask.surface.*`.

## 6. Recursive cube, Menger sponge, and porous cube families

### 6.1 Generic recursive cube-mask family

| Function | Signature | Implementation | API status (0.2.0) |
|---|---|---|---|
| `edges.recursive.cube.mask` | `edges.recursive.cube.mask(mask, level = 2)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `recursive.cube.mask.surface.embedding` | `recursive.cube.mask.surface.embedding(mask, level = 2, surface = c("standard", "bulged", "twisted", "wavy"), amplitude = 0.2, freq = 2, twist = 0.6, x_scale = 1, y_scale = 1, z_scale = 1)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `recursive.cube.mask.surface.graph` | `recursive.cube.mask.surface.graph(mask, level = 2, surface = c("standard", "bulged", "twisted", "wavy"), amplitude = 0.2, freq = 2, twist = 0.6, x_scale = 1, y_scale = 1, z_scale = 1, normalize = c("median", "mean", "none"))` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |

Parameters:
- `mask`: logical cubic keep-array.
- `level`: recursion depth.
- `surface`: deformation mode for voxel-center coordinates; `bulged`, `twisted`, and `wavy` introduce different 3D distortions.
- `amplitude`: deformation strength.
- `freq`: oscillation count for `wavy`.
- `twist`: twist intensity.
- `x_scale`, `y_scale`, `z_scale`: anisotropic scaling of the cubic parameter space.
- `normalize`: edge-weight normalization mode.

### 6.2 Menger sponge wrapper

| Function | Signature | Implementation | API status (0.2.0) |
|---|---|---|---|
| `edges.menger.sponge` | `edges.menger.sponge(level = 2)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `menger.sponge.surface.embedding` | `menger.sponge.surface.embedding(level = 2, surface = c("standard", "bulged", "twisted", "wavy"), amplitude = 0.2, freq = 2, twist = 0.6, x_scale = 1, y_scale = 1, z_scale = 1)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `menger.sponge.surface.graph` | `menger.sponge.surface.graph(level = 2, surface = c("standard", "bulged", "twisted", "wavy"), amplitude = 0.2, freq = 2, twist = 0.6, x_scale = 1, y_scale = 1, z_scale = 1, normalize = c("median", "mean", "none"))` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |

Parameters:
- `level`: Menger-sponge recursion depth.
- Remaining geometry parameters match `recursive.cube.mask.surface.*`.

### 6.3 Porous cube mask constructors

| Function | Signature | Implementation | API status (0.2.0) |
|---|---|---|---|
| `mask.cube.periodic.tunnels` | `mask.cube.periodic.tunnels(side = 5, tunnel_width = 1, tunnel_period = 2, tunnel_offset = 2)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |
| `mask.cube.asymmetric.cavities` | `mask.cube.asymmetric.cavities(side = 5, cavity_size = 2, pocket_size = max(1L, cavity_size - 1L))` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |
| `mask.cube.channel.network` | `mask.cube.channel.network(side = 5, channel_width = 1, branch_offset = 2)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |

Parameters:
- `side`: side length of the cubic mask.
- `tunnel_width`: width of removed tunnel corridors.
- `tunnel_period`: spacing between parallel tunnels.
- `tunnel_offset`: phase offset of the tunnel pattern.
- `cavity_size`: size of the main asymmetric cavity.
- `pocket_size`: size of the secondary pocket attached to the main cavity.
- `channel_width`: width of the carved channels.
- `branch_offset`: offset controlling where side branches attach in the channel network.

### 6.4 Porous cube named families

#### Periodic tunnels

| Function | Signature | Implementation | API status (0.2.0) |
|---|---|---|---|
| `edges.cube.periodic.tunnels` | `edges.cube.periodic.tunnels(level = 2, side = 5, tunnel_width = 1, tunnel_period = 2, tunnel_offset = 2)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `cube.periodic.tunnels.surface.embedding` | `cube.periodic.tunnels.surface.embedding(level = 2, side = 5, tunnel_width = 1, tunnel_period = 2, tunnel_offset = 2, surface = c("standard", "bulged", "twisted", "wavy"), amplitude = 0.2, freq = 2, twist = 0.6, x_scale = 1, y_scale = 1, z_scale = 1)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `cube.periodic.tunnels.surface.graph` | `cube.periodic.tunnels.surface.graph(level = 2, side = 5, tunnel_width = 1, tunnel_period = 2, tunnel_offset = 2, surface = c("standard", "bulged", "twisted", "wavy"), amplitude = 0.2, freq = 2, twist = 0.6, x_scale = 1, y_scale = 1, z_scale = 1, normalize = c("median", "mean", "none"))` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |

Parameters:
- `level`: recursion depth.
- `side`, `tunnel_width`, `tunnel_period`, `tunnel_offset`: same meaning as in `mask.cube.periodic.tunnels`.
- Remaining geometry parameters match `recursive.cube.mask.surface.*`.

#### Asymmetric cavities

| Function | Signature | Implementation | API status (0.2.0) |
|---|---|---|---|
| `edges.cube.asymmetric.cavities` | `edges.cube.asymmetric.cavities(level = 2, side = 5, cavity_size = 2, pocket_size = max(1L, cavity_size - 1L))` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `cube.asymmetric.cavities.surface.embedding` | `cube.asymmetric.cavities.surface.embedding(level = 2, side = 5, cavity_size = 2, pocket_size = max(1L, cavity_size - 1L), surface = c("standard", "bulged", "twisted", "wavy"), amplitude = 0.2, freq = 2, twist = 0.6, x_scale = 1, y_scale = 1, z_scale = 1)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `cube.asymmetric.cavities.surface.graph` | `cube.asymmetric.cavities.surface.graph(level = 2, side = 5, cavity_size = 2, pocket_size = max(1L, cavity_size - 1L), surface = c("standard", "bulged", "twisted", "wavy"), amplitude = 0.2, freq = 2, twist = 0.6, x_scale = 1, y_scale = 1, z_scale = 1, normalize = c("median", "mean", "none"))` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |

Parameters:
- `level`: recursion depth.
- `side`, `cavity_size`, `pocket_size`: same meaning as in `mask.cube.asymmetric.cavities`.
- Remaining geometry parameters match `recursive.cube.mask.surface.*`.

#### Channel network

| Function | Signature | Implementation | API status (0.2.0) |
|---|---|---|---|
| `edges.cube.channel.network` | `edges.cube.channel.network(level = 2, side = 5, channel_width = 1, branch_offset = 2)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `cube.channel.network.surface.embedding` | `cube.channel.network.surface.embedding(level = 2, side = 5, channel_width = 1, branch_offset = 2, surface = c("standard", "bulged", "twisted", "wavy"), amplitude = 0.2, freq = 2, twist = 0.6, x_scale = 1, y_scale = 1, z_scale = 1)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `cube.channel.network.surface.graph` | `cube.channel.network.surface.graph(level = 2, side = 5, channel_width = 1, branch_offset = 2, surface = c("standard", "bulged", "twisted", "wavy"), amplitude = 0.2, freq = 2, twist = 0.6, x_scale = 1, y_scale = 1, z_scale = 1, normalize = c("median", "mean", "none"))` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |

Parameters:
- `level`: recursion depth.
- `side`, `channel_width`, `branch_offset`: same meaning as in `mask.cube.channel.network`.
- Remaining geometry parameters match `recursive.cube.mask.surface.*`.

## 7. Triangulated manifold families

### 7.1 Closed triangulated polyhedra

| Function | Signature | Implementation | API status (0.2.0) |
|---|---|---|---|
| `edges.triangulated.polyhedron` | `edges.triangulated.polyhedron(base = c("tetrahedron", "octahedron", "icosahedron"), level = 1)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `triangulated.polyhedron.surface.embedding` | `triangulated.polyhedron.surface.embedding(base = c("tetrahedron", "octahedron", "icosahedron"), level = 1, surface = c("standard", "inflated", "twisted", "wavy"), amplitude = 0.25, freq = 2, twist = 0.6)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `triangulated.polyhedron.surface.graph` | `triangulated.polyhedron.surface.graph(base = c("tetrahedron", "octahedron", "icosahedron"), level = 1, surface = c("standard", "inflated", "twisted", "wavy"), amplitude = 0.25, freq = 2, twist = 0.6, normalize = c("median", "mean", "none"))` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |

Parameters:
- `base`: starting Platonic surface to subdivide.
- `level`: subdivision depth; each face is split recursively.
- `surface`: deformation mode on the subdivided closed surface.
- `amplitude`: deformation strength.
- `freq`: oscillation count for `wavy`.
- `twist`: twist intensity for twisted variants.
- `normalize`: edge-weight normalization mode.

### 7.2 Triangulated annulus

| Function | Signature | Implementation | API status (0.2.0) |
|---|---|---|---|
| `edges.triangulated.annulus` | `edges.triangulated.annulus(resolution = 12, outer_radius = 1, inner_radius = 0.45)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `triangulated.annulus.surface.embedding` | `triangulated.annulus.surface.embedding(resolution = 12, outer_radius = 1, inner_radius = 0.45, surface = c("flat", "saddle", "paraboloid", "ripple", "folded"), amplitude = 0.6, freq_u = 1, freq_v = 1)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `triangulated.annulus.surface.graph` | `triangulated.annulus.surface.graph(resolution = 12, outer_radius = 1, inner_radius = 0.45, surface = c("flat", "saddle", "paraboloid", "ripple", "folded"), amplitude = 0.6, freq_u = 1, freq_v = 1, normalize = c("median", "mean", "none"))` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |

Parameters:
- `resolution`: sampling density of the clipped triangular lattice.
- `outer_radius`, `inner_radius`: outer and inner annulus radii in the parameter domain.
- `surface`: surface-lift choice; `flat` stays planar, `folded` introduces a crease, the others add smooth curvature.
- `amplitude`: deformation strength.
- `freq_u`, `freq_v`: oscillation counts in the 2D parameter domain.
- `normalize`: edge-weight normalization mode.

### 7.3 Triangulated pair of pants

| Function | Signature | Implementation | API status (0.2.0) |
|---|---|---|---|
| `edges.triangulated.pair.of.pants` | `edges.triangulated.pair.of.pants(resolution = 12, outer_radius = 1.1, hole_radius = 0.24, hole_offset = 0.38, hole_height = 0.18)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `triangulated.pair.of.pants.surface.embedding` | `triangulated.pair.of.pants.surface.embedding(resolution = 12, outer_radius = 1.1, hole_radius = 0.24, hole_offset = 0.38, hole_height = 0.18, surface = c("flat", "saddle", "paraboloid", "ripple", "folded"), amplitude = 0.6, freq_u = 1, freq_v = 1)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `triangulated.pair.of.pants.surface.graph` | `triangulated.pair.of.pants.surface.graph(resolution = 12, outer_radius = 1.1, hole_radius = 0.24, hole_offset = 0.38, hole_height = 0.18, surface = c("flat", "saddle", "paraboloid", "ripple", "folded"), amplitude = 0.6, freq_u = 1, freq_v = 1, normalize = c("median", "mean", "none"))` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |

Parameters:
- `resolution`: sampling density of the clipped triangular lattice.
- `outer_radius`: radius of the outer boundary.
- `hole_radius`: radius of each leg opening.
- `hole_offset`: horizontal offset of the leg holes.
- `hole_height`: vertical placement of the upper holes relative to the waist.
- `surface`, `amplitude`, `freq_u`, `freq_v`, `normalize`: same meaning as in `triangulated.annulus.surface.*`.

## 8. Irregular point-sampled manifold families

### 8.1 Irregular annulus

| Function | Signature | Implementation | API status (0.2.0) |
|---|---|---|---|
| `edges.irregular.annulus` | `edges.irregular.annulus(rings = 6, outer_count = 28, outer_radius = 1, inner_radius = 0.45, count_irregularity = 0.2, radial_irregularity = 0.35, phase_twist = 0.35)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `irregular.annulus.surface.embedding` | `irregular.annulus.surface.embedding(rings = 6, outer_count = 28, outer_radius = 1, inner_radius = 0.45, count_irregularity = 0.2, radial_irregularity = 0.35, phase_twist = 0.35, surface = c("flat", "saddle", "paraboloid", "ripple", "folded"), amplitude = 0.6, freq_u = 1, freq_v = 1)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `irregular.annulus.surface.graph` | `irregular.annulus.surface.graph(rings = 6, outer_count = 28, outer_radius = 1, inner_radius = 0.45, count_irregularity = 0.2, radial_irregularity = 0.35, phase_twist = 0.35, surface = c("flat", "saddle", "paraboloid", "ripple", "folded"), amplitude = 0.6, freq_u = 1, freq_v = 1, normalize = c("median", "mean", "none"))` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |

Parameters:
- `rings`: number of concentric irregular rings.
- `outer_count`: nominal sample count on the outside ring.
- `outer_radius`, `inner_radius`: outer and inner annulus radii.
- `count_irregularity`: amount of variation in samples-per-ring.
- `radial_irregularity`: amount of radial perturbation.
- `phase_twist`: angular phase offset between successive rings.
- `surface`, `amplitude`, `freq_u`, `freq_v`, `normalize`: same meaning as in the triangulated annulus wrappers.

### 8.2 Irregular sphere

| Function | Signature | Implementation | API status (0.2.0) |
|---|---|---|---|
| `edges.irregular.sphere` | `edges.irregular.sphere(bands = 6, equator_count = 28, count_irregularity = 0.2, lat_irregularity = 0.35, phase_twist = 0.35)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `irregular.sphere.surface.embedding` | `irregular.sphere.surface.embedding(bands = 6, equator_count = 28, count_irregularity = 0.2, lat_irregularity = 0.35, phase_twist = 0.35, surface = c("standard", "ellipsoid", "wavy"), radius = 1, amplitude = 0.2, freq_theta = 3, freq_lat = 2, twist = 0.25)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `irregular.sphere.surface.graph` | `irregular.sphere.surface.graph(bands = 6, equator_count = 28, count_irregularity = 0.2, lat_irregularity = 0.35, phase_twist = 0.35, surface = c("standard", "ellipsoid", "wavy"), radius = 1, amplitude = 0.2, freq_theta = 3, freq_lat = 2, twist = 0.25, normalize = c("median", "mean", "none"))` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |

Parameters:
- `bands`: number of irregular latitude bands.
- `equator_count`: nominal sample count near the equator.
- `count_irregularity`: variation in band sample counts.
- `lat_irregularity`: irregularity in band placement.
- `phase_twist`: azimuthal offset between adjacent bands.
- `surface`, `radius`, `amplitude`, `freq_theta`, `freq_lat`, `twist`, `normalize`: same meaning as in the regular sphere surface family.

### 8.3 Irregular pair of pants

| Function | Signature | Implementation | API status (0.2.0) |
|---|---|---|---|
| `edges.irregular.pair.of.pants` | `edges.irregular.pair.of.pants(slices = 11, outer_count = 28, outer_radius = 1.1, hole_radius = 0.24, hole_offset = 0.38, hole_height = 0.18, count_irregularity = 0.2, vertical_irregularity = 0.35, phase_twist = 0.35)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `irregular.pair.of.pants.surface.embedding` | `irregular.pair.of.pants.surface.embedding(slices = 11, outer_count = 28, outer_radius = 1.1, hole_radius = 0.24, hole_offset = 0.38, hole_height = 0.18, count_irregularity = 0.2, vertical_irregularity = 0.35, phase_twist = 0.35, surface = c("flat", "saddle", "paraboloid", "ripple", "folded"), amplitude = 0.6, freq_u = 1, freq_v = 1)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `irregular.pair.of.pants.surface.graph` | `irregular.pair.of.pants.surface.graph(slices = 11, outer_count = 28, outer_radius = 1.1, hole_radius = 0.24, hole_offset = 0.38, hole_height = 0.18, count_irregularity = 0.2, vertical_irregularity = 0.35, phase_twist = 0.35, surface = c("flat", "saddle", "paraboloid", "ripple", "folded"), amplitude = 0.6, freq_u = 1, freq_v = 1, normalize = c("median", "mean", "none"))` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |

Parameters:
- `slices`: number of horizontal slice levels.
- `outer_count`: nominal outer-ring sample count.
- `outer_radius`, `hole_radius`, `hole_offset`, `hole_height`: geometric controls of the three-legged parameter domain.
- `count_irregularity`: variation in per-slice sample counts.
- `vertical_irregularity`: irregularity in slice spacing.
- `phase_twist`: angular offset between neighboring slices.
- `surface`, `amplitude`, `freq_u`, `freq_v`, `normalize`: same meaning as in the triangulated pair-of-pants wrappers.

### 8.4 Irregular torus

| Function | Signature | Implementation | API status (0.2.0) |
|---|---|---|---|
| `edges.irregular.torus` | `edges.irregular.torus(major_rings = 8, tube_count = 16, count_irregularity = 0.2, major_irregularity = 0.25, phase_twist = 0.35)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `irregular.torus.surface.embedding` | `irregular.torus.surface.embedding(major_rings = 8, tube_count = 16, count_irregularity = 0.2, major_irregularity = 0.25, phase_twist = 0.35, surface = c("standard", "pinched", "wavy"), major_radius = 2, minor_radius = 0.75, amplitude = 0.2, freq_major = 2, freq_minor = 1, twist = 0.25)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `irregular.torus.surface.graph` | `irregular.torus.surface.graph(major_rings = 8, tube_count = 16, count_irregularity = 0.2, major_irregularity = 0.25, phase_twist = 0.35, surface = c("standard", "pinched", "wavy"), major_radius = 2, minor_radius = 0.75, amplitude = 0.2, freq_major = 2, freq_minor = 1, twist = 0.25, normalize = c("median", "mean", "none"))` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |

Parameters:
- `major_rings`: number of irregular major-cycle rings.
- `tube_count`: nominal samples around each tube cross-section.
- `count_irregularity`: variation in per-ring sample counts.
- `major_irregularity`: irregularity of the major-angle spacing.
- `phase_twist`: angular phase shift between adjacent major rings.
- `surface`, `major_radius`, `minor_radius`, `amplitude`, `freq_major`, `freq_minor`, `twist`, `normalize`: same meaning as in the regular torus surface family.

### 8.5 Irregular double torus

| Function | Signature | Implementation | API status (0.2.0) |
|---|---|---|---|
| `edges.irregular.double.torus` | `edges.irregular.double.torus(slices = 11, tube_count = 14, branch_length = 0.85, branch_offset = 0.72, tube_radius = 0.28, transition_width = 0.42, count_irregularity = 0.2, axial_irregularity = 0.3, phase_twist = 0.35)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `irregular.double.torus.surface.embedding` | `irregular.double.torus.surface.embedding(slices = 11, tube_count = 14, branch_length = 0.85, branch_offset = 0.72, tube_radius = 0.28, transition_width = 0.42, count_irregularity = 0.2, axial_irregularity = 0.3, phase_twist = 0.35, surface = c("standard", "bulged", "twisted", "wavy"), amplitude = 0.25, freq_x = 2, freq_theta = 2, twist = 0.6)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `irregular.double.torus.surface.graph` | `irregular.double.torus.surface.graph(slices = 11, tube_count = 14, branch_length = 0.85, branch_offset = 0.72, tube_radius = 0.28, transition_width = 0.42, count_irregularity = 0.2, axial_irregularity = 0.3, phase_twist = 0.35, surface = c("standard", "bulged", "twisted", "wavy"), amplitude = 0.25, freq_x = 2, freq_theta = 2, twist = 0.6, normalize = c("median", "mean", "none"))` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |

Parameters:
- `slices`: number of slice levels in the cyclic genus-2 construction.
- `tube_count`: nominal samples around each tube component.
- `branch_length`: axial extent of the two torus lobes.
- `branch_offset`: separation of the two lobes from the center.
- `tube_radius`: radius of each handle tube.
- `transition_width`: width of the middle transition region between the two handles.
- `count_irregularity`: variation in slice/component sample counts.
- `axial_irregularity`: irregularity in slice spacing along the main axis.
- `phase_twist`: angular phase drift from slice to slice.
- `surface`: 3D deformation mode; `bulged`, `twisted`, and `wavy` modulate the base genus-2 embedding.
- `amplitude`: deformation strength.
- `freq_x`, `freq_theta`: oscillation counts along the main axis and around the tube.
- `twist`: twist intensity.
- `normalize`: edge-weight normalization mode.

## 9. Irregular tetrahedralized 3D solids

### 9.1 Irregular ball

| Function | Signature | Implementation | API status (0.2.0) |
|---|---|---|---|
| `edges.irregular.ball` | `edges.irregular.ball(base = c("tetrahedron", "octahedron", "icosahedron"), level = 1, layers = 3, outer_radius = 1, radial_irregularity = 0.25, layer_twist = 0.35)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `irregular.ball.solid.embedding` | `irregular.ball.solid.embedding(base = c("tetrahedron", "octahedron", "icosahedron"), level = 1, layers = 3, outer_radius = 1, radial_irregularity = 0.25, layer_twist = 0.35, surface = c("standard", "bulged", "twisted", "wavy"), amplitude = 0.2, freq_theta = 2, freq_phi = 2, twist = 0.6)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `irregular.ball.solid.graph` | `irregular.ball.solid.graph(base = c("tetrahedron", "octahedron", "icosahedron"), level = 1, layers = 3, outer_radius = 1, radial_irregularity = 0.25, layer_twist = 0.35, surface = c("standard", "bulged", "twisted", "wavy"), amplitude = 0.2, freq_theta = 2, freq_phi = 2, twist = 0.6, normalize = c("median", "mean", "none"))` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |

Parameters:
- `base`: outer shell template used for the layered tetrahedralized solid.
- `level`: subdivision depth of each shell.
- `layers`: number of radial shell layers; the ball also adds a center vertex.
- `outer_radius`: outermost radius.
- `radial_irregularity`: nonuniformity of intermediate shell radii.
- `layer_twist`: rotational offset between radial layers.
- `surface`: volumetric deformation mode applied to the shell-center coordinates.
- `amplitude`: deformation strength.
- `freq_theta`, `freq_phi`: angular oscillation counts.
- `twist`: twist intensity.
- `normalize`: edge-weight normalization mode.

### 9.2 Irregular shell

| Function | Signature | Implementation | API status (0.2.0) |
|---|---|---|---|
| `edges.irregular.shell` | `edges.irregular.shell(base = c("tetrahedron", "octahedron", "icosahedron"), level = 1, layers = 3, inner_radius = 0.45, outer_radius = 1, radial_irregularity = 0.25, layer_twist = 0.35)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `irregular.shell.solid.embedding` | `irregular.shell.solid.embedding(base = c("tetrahedron", "octahedron", "icosahedron"), level = 1, layers = 3, inner_radius = 0.45, outer_radius = 1, radial_irregularity = 0.25, layer_twist = 0.35, surface = c("standard", "bulged", "twisted", "wavy"), amplitude = 0.2, freq_theta = 2, freq_phi = 2, twist = 0.6)` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | internal |
| `irregular.shell.solid.graph` | `irregular.shell.solid.graph(base = c("tetrahedron", "octahedron", "icosahedron"), level = 1, layers = 3, inner_radius = 0.45, outer_radius = 1, radial_irregularity = 0.25, layer_twist = 0.35, surface = c("standard", "bulged", "twisted", "wavy"), amplitude = 0.2, freq_theta = 2, freq_phi = 2, twist = 0.6, normalize = c("median", "mean", "none"))` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |

Parameters:
- `base`, `level`, `layers`: same meaning as in `edges.irregular.ball`.
- `inner_radius`, `outer_radius`: inner and outer shell radii.
- `radial_irregularity`: nonuniformity of intermediate shell radii.
- `layer_twist`: rotational offset between layers.
- `surface`, `amplitude`, `freq_theta`, `freq_phi`, `twist`, `normalize`: same meaning as in `irregular.ball.solid.*`.

## 10. Intrinsic weighted trees

| Function | Signature | Implementation | API status (0.2.0) |
|---|---|---|---|
| `kary.tree.weighted.graph` | `kary.tree.weighted.graph(k = 2, depth = 2, base_length = 1, depth_rule = c("geometric", "constant", "custom"), depth_decay = 0.85, depth_factors = NULL, branch_rule = c("linear", "uniform", "custom"), branch_spread = 0.3, branch_factors = NULL, normalize = c("median", "mean", "none"))` | [graph_helpers.R](https://github.com/pgajer/grip/blob/v0.2.0/R/graph_helpers.R) | exported |

Parameters:
- `k`: branching factor.
- `depth`: tree depth.
- `base_length`: base multiplicative scale for all edge lengths.
- `depth_rule`: how edge length changes with depth:
  - `"geometric"` uses `depth_decay^depth`,
  - `"constant"` uses no depth taper,
  - `"custom"` uses `depth_factors`.
- `depth_decay`: geometric taper factor used when `depth_rule = "geometric"`.
- `depth_factors`: explicit per-depth multipliers for `depth_rule = "custom"`.
- `branch_rule`: how sibling branches vary:
  - `"linear"` spreads lengths across child slots,
  - `"uniform"` keeps all siblings equal,
  - `"custom"` uses `branch_factors`.
- `branch_spread`: spread magnitude for the `"linear"` sibling rule.
- `branch_factors`: explicit per-child-slot multipliers for `branch_rule = "custom"`.
- `normalize`: edge-weight normalization mode.

## 11. Quick index by family

- Regular lifted lattices: `mesh.surface.*`, `cylinder.surface.*`, `torus.surface.*`, `sphere.surface.*`
- Recursive square masks: `mask.*`, `edges.recursive.mask.grid`, `recursive.mask.grid.surface.*`, `edges.sierpinski.carpet`, `sierpinski.carpet.surface.*`, `edges.vicsek`, `vicsek.surface.*`
- Perforated meshes: `keep.*`, `edges.occupied.mesh`, `occupied.mesh.surface.*`
- Recursive triangle gaskets: `mask.triangle.*`, `edges.recursive.triangle.mask`, `recursive.triangle.mask.surface.*`, `edges.sierpinski.triangle`, `sierpinski.triangle.surface.*`
- Recursive tetrahedron gaskets: `mask.tetrahedron.*`, `edges.recursive.tetrahedron.mask`, `recursive.tetrahedron.mask.surface.*`, `edges.sierpinski.tetrahedron`, `sierpinski.tetrahedron.surface.*`
- Recursive cubes and porous cubes: `edges.recursive.cube.mask`, `recursive.cube.mask.surface.*`, `edges.menger.sponge`, `menger.sponge.surface.*`, `mask.cube.*`, `edges.cube.*`, `cube.*.surface.*`
- Triangulated manifolds: `edges.triangulated.*`, `triangulated.*.surface.*`
- Irregular sampled manifolds: `edges.irregular.*`, `irregular.*.surface.*`
- Volumetric solids: `edges.irregular.ball`, `irregular.ball.solid.*`, `edges.irregular.shell`, `irregular.shell.solid.*`
- Intrinsic weighted trees: `kary.tree.weighted.graph`
