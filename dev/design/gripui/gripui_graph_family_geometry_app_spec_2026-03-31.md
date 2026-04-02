# Spec: Shiny App for Exploring Synthetic Graph Family Geometries

This document specifies a Shiny app for interactively exploring the synthetic graph families implemented in [R/graph_helpers.R](/Users/pgajer/current_projects/grip/R/graph_helpers.R). The app should let a user browse all families, adjust family-specific topology and geometry parameters, and inspect the resulting graph geometry through embedded `rglwidget` 3D views.

The recommendation is to build this as a geometry-first sibling of the existing `gripui` app rather than forcing it into the current layout-catalog workflow. The existing `gripui` infrastructure is still highly reusable:

- app entry pattern in [R/gripui-app.R](/Users/pgajer/current_projects/grip/R/gripui-app.R)
- theme/layout structure in [R/gripui-ui.R](/Users/pgajer/current_projects/grip/R/gripui-ui.R)
- server organization in [R/gripui-server.R](/Users/pgajer/current_projects/grip/R/gripui-server.R)
- `rglwidget` rendering in [R/gripui-render.R](/Users/pgajer/current_projects/grip/R/gripui-render.R)
- CSS in [inst/app/www/gripui.css](/Users/pgajer/current_projects/grip/inst/app/www/gripui.css)

## 1. Goal

Provide a Shiny app that allows a user to:

- browse every synthetic graph-family branch added in the thread
- adjust all relevant graph-topology and geometry parameters
- see only the controls that make sense for the selected family
- inspect the generated graph in 3D with an embedded `rglwidget`
- optionally inspect 2D parameter coordinates, graph summaries, and edge-weight summaries
- compare presets or parameter variants within and across families
- export reproducible R code for the currently selected example

## 2. Recommended app identity

Recommended exported API:

- `gripui_family_app(catalog = gripui_graph_family_catalog(), title = "Graph Family Geometry Explorer")`
- `run_gripui_family_app(...)`

Recommended positioning:

- Keep existing `gripui_*` for layout-catalog inspection.
- Add a new family-explorer app for synthetic graph generation.
- Reuse the same visual language, but keep the data model separate.

Reason:

- current `gripui` assumes a static catalog of already-generated layouts
- the new app must generate graphs reactively from parameterized constructors
- forcing both concerns into one project object would make the existing UI harder to maintain

## 3. Families the app must cover

The app should cover all family branches added in the thread.

### 3.1 Regular lifted lattice families

- `mesh.surface.graph()`
- `cylinder.surface.graph()`
- `torus.surface.graph()`
- `sphere.surface.graph()`

### 3.2 Recursive square-mask families

- `recursive.mask.grid.surface.graph()`
- `sierpinski.carpet.surface.graph()`
- `vicsek.surface.graph()`
- mask constructors:
  - `mask.cross()`
  - `mask.border()`
  - `mask.corner()`
  - `mask.asymmetric.holes()`

### 3.3 Deterministic perforated meshes

- `occupied.mesh.surface.graph()`
- occupancy-pattern helpers:
  - `keep.periodic.holes()`
  - `keep.staggered.windows()`
  - `keep.slit.channels()`
  - `keep.asymmetric.notches()`

### 3.4 Recursive triangle-mask families

- `recursive.triangle.mask.surface.graph()`
- `sierpinski.triangle.surface.graph()`
- mask helpers:
  - `mask.triangle.classic()`
  - `mask.triangle.bridge()`

### 3.5 Recursive tetrahedron-mask families

- `recursive.tetrahedron.mask.surface.graph()`
- `sierpinski.tetrahedron.surface.graph()`
- mask helpers:
  - `mask.tetrahedron.classic()`
  - `mask.tetrahedron.corner.missing()`

### 3.6 Recursive cube and porous cube families

- `recursive.cube.mask.surface.graph()`
- `menger.sponge.surface.graph()`
- `cube.periodic.tunnels.surface.graph()`
- `cube.asymmetric.cavities.surface.graph()`
- `cube.channel.network.surface.graph()`
- cube mask helpers:
  - `mask.cube.periodic.tunnels()`
  - `mask.cube.asymmetric.cavities()`
  - `mask.cube.channel.network()`

### 3.7 Triangulated manifold families

- `triangulated.polyhedron.surface.graph()`
- `triangulated.annulus.surface.graph()`
- `triangulated.pair.of.pants.surface.graph()`

### 3.8 Irregular point-sampled manifold families

- `irregular.annulus.surface.graph()`
- `irregular.sphere.surface.graph()`
- `irregular.pair.of.pants.surface.graph()`
- `irregular.torus.surface.graph()`
- `irregular.double.torus.surface.graph()`

### 3.9 Volumetric families

- `irregular.ball.solid.graph()`
- `irregular.shell.solid.graph()`

### 3.10 Intrinsic weighted trees

- `kary.tree.weighted.graph()`

## 4. Core UX

### 4.1 Main layout

Recommended page layout:

- left sidebar:
  - category selector
  - family selector
  - preset selector
  - dynamic parameter panel
  - render/apply button
  - reset-to-defaults button
- main body:
  - top row:
    - 3D geometry viewer
    - parameter-domain or auxiliary 2D view
  - middle row:
    - graph summary card
    - edge-weight summary card
    - function call / reproducible-code card
  - bottom row:
    - compare/preset grid or family notes

### 4.2 Suggested tabs

Recommended tabs:

- `Explore`
- `Compare`
- `Gallery`
- `About`

Purpose:

- `Explore`: single-family interactive builder
- `Compare`: side-by-side variants of one family or multiple families
- `Gallery`: curated presets across all families, effectively an in-app version of the R Markdown gallery
- `About`: family descriptions, implementation references, and tips

## 5. Data model

The app should be driven by a registry, not by hard-coded conditional UI logic scattered across the server.

Recommended internal constructor:

- `gripui_graph_family_catalog()`

This should return a list of family descriptors, one per UI-selectable family.

### 5.1 Family descriptor schema

Each family descriptor should include:

- `id`
- `label`
- `category`
- `subcategory`
- `kind`
  - one of: `surface_graph`, `solid_graph`, `weighted_graph`
- `builder`
  - function that returns the graph bundle shown in the app
- `topology_builder`
  - optional topology-only builder if separate from geometry builder
- `params`
  - list of parameter specs
- `presets`
  - named preset list
- `summary`
  - one-paragraph description
- `implementation_ref`
  - source location in [R/graph_helpers.R](/Users/pgajer/current_projects/grip/R/graph_helpers.R)
- `output_adapter`
  - function that normalizes builder output into a standard app payload
- `supports_param_view`
  - whether a useful `coords_param` view exists
- `supports_weight_hist`
  - whether edge weights are meaningful to display
- `tags`
  - optional vector like `c("fractal", "surface", "closed", "genus1", "porous")`

### 5.2 Standard app payload

Every family should be normalized into a common payload object:

```r
list(
  family_id = "irregular.torus",
  family_label = "Irregular torus",
  edges = <integer matrix>,
  n = <integer>,
  edge_weights = <numeric vector or NULL>,
  coords_surface = <numeric matrix or NULL>,
  coords_param = <numeric matrix or NULL>,
  coords_display = <numeric matrix>,
  label = <character>,
  note = <character or NULL>,
  summary = list(
    n = ...,
    m = ...,
    weight_min = ...,
    weight_max = ...,
    weight_cv = ...,
    bbox = ...
  ),
  code = "reproducible R call as text"
)
```

Rules:

- `coords_display` is the matrix used by the 3D viewer.
- for most families, `coords_display = coords_surface`
- for `kary.tree.weighted.graph()`, `coords_display` must be a synthetic 3D embedding derived from the tree topology and edge weights

## 6. Dynamic parameter system

The app’s main feature is that the parameter panel adapts to the selected family. That should be metadata-driven.

### 6.1 Parameter descriptor schema

Each parameter spec should include:

- `id`
- `label`
- `type`
  - `integer`
  - `double`
  - `choice`
  - `logical`
  - `choice_set`
  - `vector_numeric`
  - `mask_square`
  - `mask_triangle`
  - `mask_tetrahedron`
  - `mask_cube`
- `default`
- `choices`
  - for `choice` inputs
- `min`
- `max`
- `step`
- `help`
- `group`
  - such as `Topology`, `Geometry`, `Weights`, `Mask`, `Rendering`
- `visible_if`
  - function of current parameter state
- `validate`
  - validator function returning `NULL` or an error string
- `coerce`
  - coercion function from UI value to builder-ready value
- `advanced`
  - whether hidden by default behind an “Advanced” toggle

### 6.2 Parameter groups

Recommended groups:

- `Topology`
- `Mask`
- `Geometry`
- `Weights`
- `Rendering`
- `Presets`

### 6.3 Dependent visibility rules

Examples:

- `freq_u` and `freq_v` should be shown only when the chosen `surface` uses them meaningfully
- `twist` should be shown only for families/surfaces that expose twist
- `normalize` should be hidden in pure-embedding preview mode if the user is viewing only coordinates
- `depth_factors` should be shown only if `depth_rule == "custom"`
- `branch_factors` should be shown only if `branch_rule == "custom"`
- mask editors should appear only for generic recursive families

## 7. Specialized editors for non-scalar parameters

Several families cannot be handled well by simple sliders alone.

### 7.1 Square mask editor

Used by:

- `recursive.mask.grid.surface.graph()`

UI:

- `k` selector
- `mode` selector:
  - `named`
  - `custom`
- named choices:
  - `cross`
  - `border`
  - `corner`
  - `asymmetric_holes`
- for `custom`, show a `k x k` clickable checkbox grid

Validation:

- nonempty
- connected

### 7.2 Triangle mask editor

Used by:

- `recursive.triangle.mask.surface.graph()`

UI:

- four toggles:
  - `left`
  - `right`
  - `top`
  - `center`
- preset shortcuts:
  - `classic`
  - `bridge_top`
  - `bridge_left`
  - `bridge_right`

Validation:

- at least one retained subtriangle

### 7.3 Tetrahedron mask editor

Used by:

- `recursive.tetrahedron.mask.surface.graph()`

UI:

- four toggles:
  - `apex`
  - `base_left`
  - `base_right`
  - `base_back`
- preset shortcuts:
  - `classic`
  - `missing_apex`
  - `missing_base_left`
  - `missing_base_right`
  - `missing_base_back`

Validation:

- at least one retained corner

### 7.4 Cube mask editor

Used by:

- `recursive.cube.mask.surface.graph()`

UI:

- `k` selector
- z-slice tabs
- clickable `k x k` voxel occupancy grid per slice
- preset shortcuts for:
  - periodic tunnels
  - asymmetric cavities
  - channel network
  - Menger sponge mask

Validation:

- nonempty
- connected

Implementation note:

- the generic cube editor may be Phase 2 if Phase 1 ships with named porous families only

### 7.5 Occupied-mesh editor

Used by:

- `occupied.mesh.surface.graph()`

Recommendation:

- do not expose a raw `keep` matrix editor in Phase 1
- instead expose deterministic pattern families:
  - periodic holes
  - staggered windows
  - slit channels
  - asymmetric notches

This is simpler and still covers the full implemented perforated-family branch.

## 8. Rendering and views

### 8.1 Primary 3D view

Use a refactored version of `gripui.render.rglwidget()` from [R/gripui-render.R](/Users/pgajer/current_projects/grip/R/gripui-render.R).

Recommended changes:

- accept a standardized graph-family payload instead of layout coordinates only
- support edge coloring by:
  - uniform
  - edge weight
  - branch index for trees
- support vertex coloring by:
  - z-value
  - degree
  - depth or band index when available
- support edge toggle
- support point size and edge-width sliders

Recommended helper:

- `gripui.render.geometry.rglwidget(payload, color_by = "z", show_edges = TRUE, show_vertices = TRUE, edge_color_by = "weight")`

### 8.2 Secondary view

Show one of:

- `coords_param` scatter plot when available
- 2D slice view for volumetric families
- occupancy/mask preview for recursive mask families
- schematic tree view for intrinsic weighted trees

This secondary view is important because many families distinguish between:

- intrinsic graph/weight structure
- parameter-domain structure
- extrinsic 3D embedding used to induce weights

### 8.3 Summary cards

Recommended cards:

- `Graph summary`
  - family id
  - vertex count
  - edge count
  - connectedness
  - min/mean/max degree
- `Weight summary`
  - min, median, mean, max
  - coefficient of variation
  - normalization mode
- `Geometry summary`
  - bounding box
  - dimensionality
  - whether `coords_param` exists
- `Reproducible code`
  - exact call used to create the current object

## 9. Compare mode

The app should support comparison, but it should remain lighter than the main layout-comparison `gripui`.

Recommended compare modes:

- compare presets within one family
- compare the same family under multiple geometry settings
- compare multiple families with one representative preset each

Comparison display:

- 2x2 or 3x2 grid of `rglwidget` panels
- shared summary table below

Useful compare presets:

- same topology, different geometry
- same family, different recursion depth
- regular vs irregular variant
- open vs closed surface

## 10. Gallery mode

The app should contain an in-app gallery sourced from a fixed list of presets.

This gallery should mirror the generated R Markdown/HTML gallery:

- [graph_geometry_gallery_2026-03-31.Rmd](/Users/pgajer/current_projects/grip/dev/design/graph_geometry_gallery_2026-03-31.Rmd)
- [graph_geometry_gallery_2026-03-31.html](/Users/pgajer/current_projects/grip/dev/design/graph_geometry_gallery_2026-03-31.html)

Recommendation:

- build the gallery preset list once in code
- reuse it in both:
  - the Shiny app `Gallery` tab
  - the R Markdown gallery document

That avoids duplicated preset definitions drifting over time.

## 11. Implementation architecture

Recommended new files:

- [R/gripui-family-catalog.R](/Users/pgajer/current_projects/grip/R/gripui-family-catalog.R)
  - family registry
  - parameter metadata
  - preset definitions
- [R/gripui-family-normalize.R](/Users/pgajer/current_projects/grip/R/gripui-family-normalize.R)
  - normalize family outputs to the common payload
  - tree display embedding helper
- [R/gripui-family-ui.R](/Users/pgajer/current_projects/grip/R/gripui-family-ui.R)
  - app UI
  - dynamic parameter module UI
- [R/gripui-family-server.R](/Users/pgajer/current_projects/grip/R/gripui-family-server.R)
  - server logic
  - reactive build pipeline
- [R/gripui-family-app.R](/Users/pgajer/current_projects/grip/R/gripui-family-app.R)
  - exported app entry points
- [R/gripui-family-render.R](/Users/pgajer/current_projects/grip/R/gripui-family-render.R)
  - 3D viewer
  - parameter-domain plot
  - summary-card helpers
- [inst/app/www/gripui-family.css](/Users/pgajer/current_projects/grip/inst/app/www/gripui-family.css)
  - app-specific styling

### 11.1 Reuse from existing `gripui`

Direct reuse candidates:

- `gripui.enable.rgl.null.device()`
- `gripui.require.app.packages()`
- theme choices from `gripui_ui()`
- rendering patterns from `gripui.render.rglwidget()`

Possible refactor:

- move generic `rglwidget` rendering pieces into shared helpers
- keep layout-specific and family-specific server logic separate

## 12. Reactive flow

Recommended server pipeline:

1. user selects category and family
2. app loads the family descriptor
3. UI module renders only parameters for that family
4. user edits controls
5. app validates/coerces parameters
6. app builds the graph bundle
7. output adapter normalizes the result
8. viewer cards update
9. compare/gallery tabs consume the same normalized payloads

Important implementation detail:

- use an explicit `Render` button or debounced reactivity

Reason:

- some families can get moderately expensive as recursion depth rises
- immediate rebuild on every slider move will feel noisy and slow

Recommended behavior:

- lightweight fields can update reactively with a 250–500 ms debounce
- expensive families or high-depth settings should require explicit apply

## 13. Caching

The app should cache built graph payloads by a hash of:

- family id
- parameter values
- display mode if relevant

Recommended cache key:

- `digest::digest(list(family_id, params))`

Cache target:

- normalized payload, not only raw family output

Reason:

- same object may be reused across Explore, Compare, and Gallery tabs

## 14. Validation rules

The app must fail gently and explain what is wrong in plain language.

Examples:

- recursion depth too high for current memory budget
- custom mask is empty
- custom mask is disconnected
- `inner_radius >= outer_radius`
- `layers < 2` for shell families
- `length(depth_factors) < depth`
- `length(branch_factors) != k`

Recommended UI behavior:

- inline validation message next to the relevant control
- summary alert at top of viewer area
- preserve last successful render until the new parameter set becomes valid

## 15. Suggested presets

Each family descriptor should include 2–5 named presets.

Examples:

- mesh:
  - `saddle_small`
  - `paraboloid_wide`
  - `ripple_dense`
- recursive mask grid:
  - `cross`
  - `border`
  - `corner_top_left`
  - `asymmetric`
- perforated mesh:
  - `periodic_holes`
  - `staggered_windows`
  - `vertical_slits`
  - `asymmetric_notches`
- irregular torus:
  - `pinched`
  - `wavy`
- irregular double torus:
  - `twisted`
  - `bulged`
- weighted tree:
  - `geometric_linear`
  - `custom_asymmetric`

## 16. Export features

Recommended export actions:

- copy exact R call
- download graph payload as `.rds`
- download edge list with weights as `.csv`
- download current viewer screenshot as `.png`

Optional later:

- export a mini HTML widget for the current family

## 17. Testing strategy

Recommended tests:

- registry coverage test:
  - every family descriptor builds successfully
- parameter visibility tests:
  - family-specific controls appear/disappear correctly
- normalization tests:
  - all builders produce a valid common payload
- tree-display tests:
  - weighted-tree display preserves edge-length ordering
- Shiny snapshot tests:
  - one snapshot per major family category

Recommended files:

- [tests/testthat/test-gripui-family-catalog.R](/Users/pgajer/current_projects/grip/tests/testthat/test-gripui-family-catalog.R)
- [tests/testthat/test-gripui-family-normalize.R](/Users/pgajer/current_projects/grip/tests/testthat/test-gripui-family-normalize.R)
- [tests/testthat/test-gripui-family-app.R](/Users/pgajer/current_projects/grip/tests/testthat/test-gripui-family-app.R)

## 18. Phased implementation plan

### Phase 1

- create family registry
- implement Explore tab
- cover all named families
- reuse `rglwidget` renderer
- add parameter-domain/summary panels

### Phase 2

- add generic mask editors
- add Compare tab
- add caching
- add reproducible-code export

### Phase 3

- add Gallery tab backed by shared preset registry
- add download/export helpers
- add more polished family notes and cross-links to docs

## 19. Recommended design decisions

- Build a sibling app, not a forced extension of the current layout-catalog `gripui`.
- Drive the UI from a registry with explicit parameter metadata.
- Use specialized editors for mask families instead of overloading plain text inputs.
- Normalize all family outputs to one internal payload format.
- Reuse `rglwidget` rendering and the current `gripui` look-and-feel.
- Use debounced or explicit rebuilds, not continuous recomputation on every keystroke.
- Share preset definitions between the app and the generated gallery document.

## 20. Most important implementation risk

The biggest risk is not rendering, but drift between:

- family constructors
- family parameter metadata
- gallery presets
- compare presets

The best mitigation is to make the family registry the single source of truth. The app, compare mode, and static gallery should all read from that same registry.
