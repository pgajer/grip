# gripui Design for `grip`

Date: 2026-03-24

## Summary

`gripui` should live inside the `grip` package as an optional Shiny layer for
exploring the layout landscape of a single graph.

The core idea is:

- `grip` remains the package that computes layouts, scores them, and stores
  graph data.
- `gripui` becomes the package-facing exploration surface for comparing many
  realized layouts of the same graph.
- the UI must work in two modes:
  - in-memory mode for direct use after `grip.compare.layouts()`
  - file-backed mode for larger projects such as HMP/U01
- file-backed ingestion should converge on one canonical bundle contract;
  project-specific table/manifest merges are compatibility adapters, not the
  target steady state

This avoids splitting the project into a second package while still keeping the
UI code isolated from the core layout code.

## Why it belongs in `grip`

- The main abstractions already exist in `grip`:
  `grip.layout()`, `grip.layout.trace()`, `grip.score.layout()`, and
  `grip.compare.layouts()`.
- The HMP/U01 coarsened graph is already packaged in `grip`.
- The UI is not a general graph-analysis workflow like `gflowui`; it is a thin
  exploration layer over `grip` outputs.
- Keeping it in `grip` makes it easier to stabilize one artifact contract for:
  graph data, layout metrics, realized coordinates, and saved visual assets.

## Product goals

- explore many layouts for one graph without re-running expensive jobs
- filter and rank layouts by parameters and quality metrics
- load a selected layout quickly in an interactive 2D or 3D viewer
- compare a small number of layouts side by side
- gracefully support projects where some candidates have metrics only and only a
  subset have saved coordinates or HTML artifacts

## Non-goals

- general graph construction
- endpoint detection, arm workflows, or data-ingest workflows from `gflowui`
- distributed job management
- automatic parsing of free-form handoff notes

## Package and dependency plan

Keep core `grip` dependencies unchanged as much as possible.

Add UI dependencies to `Suggests`:

- `shiny`
- `bslib`
- `DT`
- `htmltools`
- `rgl`

Optional later:

- `plotly` if we want a plotly-based metric scatter with linked brushing

The exported UI entry points should fail with a clear message if suggested UI
packages are missing.

## Public API

### Phase 1 exports

These are the functions I would add first.

#### `gripui_app(project, ...)`

Build a Shiny app object from a normalized `gripui_project`.

#### `run_gripui(project, host = "127.0.0.1", port = getOption("shiny.port"), launch.browser = interactive(), ...)`

Run the app.

#### `gripui_project(graph, layouts, title = NULL, subtitle = NULL, notes = NULL)`

Create a normalized in-memory project object.

Expected fields:

- `graph`
  a list with `adj_list`, optional `weight_list`, optional `vertex_data`,
  optional `graph_info`
- `layouts`
  a data frame with one row per realized layout
- optional list-columns or columns for artifacts such as `coords`, `coords_path`,
  `html_path`, `gif_path`, `thumbnail_path`

#### `gripui_project_from_compare(compare_obj, graph, vertex_data = NULL, graph_info = NULL, title = NULL)`

Create a `gripui_project` directly from `grip.compare.layouts(..., return.layouts = TRUE)`.

This is the easiest package-facing workflow for small and medium examples.

#### `gripui_project_from_dir(root, graph = NULL, title = NULL, subtitle = NULL)`

Create a `gripui_project` from a file-backed project directory or from a set of
 known manifests and run tables under `root`.

This is the main entry point for large saved searches such as HMP/U01.

#### `gripui_validate_project(project)`

Validate structure and required columns before app launch.

### Phase 2 exports

These are useful once the app works.

#### `gripui_write_bundle(project, path)`

Write a portable on-disk bundle.

#### `gripui_build_catalog(...)`

Build a layout catalog data frame from one or more run tables plus optional
artifact manifests.

#### `gripui_export_layout_artifacts(...)`

Save embeddings, thumbnails, and optional interactive HTML for one or more
layouts after a search run.

This is the missing link for the current HMP/U01 workflow.

## `gripui_project` object shape

Use a plain list with a lightweight class:

```r
list(
  graph = list(
    adj_list = ...,
    weight_list = ...,
    vertex_data = ...,
    graph_info = ...
  ),
  layouts = data.frame(...),
  meta = list(
    title = ...,
    subtitle = ...,
    notes = ...
  )
)
```

Class:

```r
class(x) <- "gripui_project"
```

Required `layouts` columns:

- `layout_id`
- `candidate`
- `stage`
- `seed`
- `status`
- `viewable`

Recommended parameter columns:

- `placement`
- `rounds`
- `final_rounds`
- `num_init`
- `num_nbrs`
- `r`
- `s`
- `repulsion_factor`
- `tinit_factor`

Recommended metric columns:

- `sampled_stress`
- `edge_length_cv`
- `sampled_nonedge_sep_ratio`
- `cluster_separation`
- `cst_cluster_separation`
- `subcst_cluster_separation`
- `dcst_depth1_absorb_cluster_separation`
- `dcst_depth2_absorb_cluster_separation`
- `stability_procrustes_mean`
- `score_composite`
- `score_composite_extended`
- `run_score_extended`
- `elapsed_sec`

Recommended artifact columns:

- `coords`
- `coords_path`
- `html_path`
- `gif_path`
- `thumbnail_path`
- `color_view_default`

The app should treat `coords` and `coords_path` as interchangeable sources.

## File-backed project contract

For large external runs, I would standardize on a simple directory contract.

```text
project_root/
  graph/
    graph.rds
  catalog/
    layout_catalog.tsv
  artifacts/
    layouts/<layout_id>/embedding.tsv
    html/<layout_id>.html
    gif/<layout_id>.gif
    png/<layout_id>.png
  metadata/
    notes.md
    source_files.tsv
```

### `graph/graph.rds`

An RDS containing:

- `adj_list`
- optional `weight_list`
- optional `vertex_data`
- optional `graph_info`

### `catalog/layout_catalog.tsv`

One row per realized layout.

Required columns:

- `layout_id`
- `candidate`
- `stage`
- `seed`
- `status`
- `viewable`

Recommended columns:

- all tuning parameters
- all summary metrics
- artifact paths relative to `project_root`

### Why a catalog matters

The app should never scrape directories or parse prose to infer state. It
should read one normalized table and then load optional assets lazily.

## Ingestion strategy

The intended ingestion architecture is:

1. one canonical internal model:
   `gripui_project`
2. one canonical on-disk bundle contract:
   `graph/graph.rds` plus `catalog/layout_catalog.tsv` plus artifact files
3. one or more compatibility adapters for older project layouts

Under this plan:

- the current HMP/U01 coarse/full/repulsion merge is a useful first adapter
- the app itself should only consume the normalized `gripui_project`
- future exporters should write the canonical bundle directly so `gripui` does
  not need project-specific merge logic

This distinction is important. The HMP adapter is a development bootstrap and a
real-world test case, but it should not become the preferred producer contract.

## Proposed files inside `grip`

### `R/gripui-package.R`

Small package-level helpers and roxygen block for the `gripui` feature.

### `R/gripui-project.R`

Functions:

- `gripui_project()`
- `gripui_project_from_compare()`
- `gripui_project_from_dir()`
- `gripui_validate_project()`

This file may contain compatibility loaders for legacy project layouts, but the
goal is for those loaders to normalize everything into the same
`gripui_project`.

### `R/gripui-catalog.R`

Functions:

- `gripui_build_catalog()`
- `gripui_merge_run_tables()`
- `gripui_attach_artifact_manifest()`
- `gripui_normalize_catalog_names()`
- `gripui_make_layout_id()`

This file contains generic catalog normalization helpers. Any HMP/U01-specific
joining logic should stay thin and should adapt into this generic layer rather
than define a second contract.

### `R/gripui-io.R`

Functions:

- `gripui_read_embedding_tsv()`
- `gripui_read_graph_rds()`
- `gripui_resolve_project_path()`
- `gripui_load_layout_coords()`

### `R/gripui-render.R`

Functions:

- `gripui_render_rglwidget()`
- `gripui_render_saved_html()`
- `gripui_project_coords_2d()`
- `gripui_make_vertex_colors()`
- `gripui_make_edge_segments()`

This file should reuse `grip` graph helpers where practical.

### `R/gripui-app.R`

Functions:

- `gripui_app()`
- `run_gripui()`

### `R/gripui-ui.R`

Functions:

- `gripui_ui()`
- `gripui_sidebar_ui()`
- `gripui_main_ui()`

### `R/gripui-server.R`

Functions:

- `gripui_server()`
- `gripui_catalog_server()`
- `gripui_viewer_server()`
- `gripui_compare_server()`

### `inst/app/www/gripui.css`

Minimal styling only. Keep the UI clean and analytical rather than ornamental.

### Tests

- `tests/testthat/test-gripui-project.R`
- `tests/testthat/test-gripui-catalog.R`
- `tests/testthat/test-gripui-io.R`
- `tests/testthat/test-gripui-app.R`

## App layout

Use a simpler layout than `gflowui`.

### Left sidebar

- project title and graph summary
- stage filter
- candidate search box
- parameter filters
- metric filters
- color-by selector
- edge toggle
- seed selector

### Main panel

Top row:

- landscape plot
- selected layout viewer

Bottom row:

- sortable catalog table
- selected-layout metadata card

### Compare drawer

Allow the user to pin up to 4 layouts and compare them side by side.

Key features:

- shared color mapping
- shared camera
- optional Procrustes alignment to a reference layout
- quick metric deltas against the reference

## Viewer behavior

### Default rendering

- use interactive `rglwidget` for coordinates available in memory or on disk
- use saved HTML only as a fallback or convenience link

Reason:

- saved HTML files are useful artifacts, but the app should not depend on them
  for core viewing
- coordinate-driven rendering lets us keep camera linking, recoloring, and
  side-by-side comparison inside the app

### Color views

Support:

- `plain`
- categorical views from `vertex_data`, for example `cst`, `subcst`,
  `dcst.depth1.absorb`, `dcst.depth2.absorb`
- continuous views such as `ph` or `log10_reads`

### Layout availability states

Each row should clearly show one of:

- `interactive`
  coordinates available, can render inside app
- `artifact-only`
  saved HTML or image exists but raw coordinates are absent
- `summary-only`
  metrics exist but no layout asset is available

This distinction is important for the current HMP/U01 outputs.

## Landscape plot

The landscape plot is central. It should not just be a table.

### MVP

A 2D scatter where the user chooses:

- x axis from metrics or tuning parameters
- y axis from metrics or tuning parameters
- point color by stage or one selected metric
- point shape by stage or seed

Clicking a point selects the layout row and opens it in the viewer.

### Later

- pareto-front highlighting
- brushing to compare selected subsets
- density overlays for larger searches

## HMP/U01 mapping

The current HMP/U01 artifacts naturally map into the proposed contract.

### Graph source

Use the packaged graph from `grip` or the upstream coarsening files:

- `data/hmp.u01.gc.coarse.rda`
- `inst/extdata/hmp_u01_gc_coarse/*`

### Catalog sources

The catalog can be built from:

- `coarse_stage/coarse_stage_runs.tsv`
- `full_stage1/full_stage1_runs.tsv`
- `full_stage2/full_stage2_runs.tsv`
- `repulsion_sweep_runs.tsv`

### Artifact manifests

Join optional paths from:

- `top_layout_visual_manifest.tsv`
- `repulsion_sweep_visual_manifest.tsv`

### Status of this loader

This HMP/U01 loader should be treated as a compatibility adapter for an older
output shape.

Its role is:

- prove that the canonical `gripui_project` model is expressive enough for a
  real project
- let us build and validate the MVP app on real artifacts
- bridge the gap until search/export scripts can emit canonical bundles

It should not be treated as the desired long-term producer contract.

### Important current limitation

For the 2026-03-23 HMP/U01 outputs:

- all candidates have summary metrics
- only the top 3 finalists have saved embeddings and visual assets in the main
  coarsened-search output
- only the 6 repulsion-sweep layouts have saved embeddings in the sweep output

So the first HMP/U01 app should expose:

- full landscape browsing across all summary rows
- full interactive viewing for the 9 layouts with saved embeddings
- disabled viewer controls for summary-only rows

## Canonical producer path after the MVP

After the MVP app is working on the HMP adapter, the next ingestion milestone
should be to make canonical bundle writing the preferred path.

That means:

- keep `gripui_project_from_dir()` able to read the HMP adapter path
- add bundle-writing helpers so new searches can write the canonical layout
  catalog and artifact tree directly
- update future analysis/export scripts to produce canonical bundles first and
  compatibility adapters only when needed

In other words, the HMP merge path comes first chronologically, but canonical
bundle writing is the next desired steady-state ingestion target.

## Changes needed in future search/export scripts

To make `gripui` fully useful for large searches, future exporters should save
coordinates for every successful run, not only for finalists.

### Minimum artifact contract per successful run

- `embedding.tsv`
- one thumbnail PNG
- one catalog row

Optional:

- one saved interactive HTML
- one GIF

### Recommended helper to add in `grip`

`gripui_export_layout_artifacts()` should take:

- a layout matrix
- vertex ids
- optional graph and metadata
- output root
- `layout_id`

and write:

- `embedding.tsv`
- thumbnail PNG
- optional HTML widget

That lets analysis scripts produce app-ready outputs consistently.

## Relationship to existing `grip` code

`gripui` should sit on top of existing code rather than duplicating logic.

Reuse directly:

- `grip.compare.layouts()`
- `grip.score.layout()`
- `grip.layout.trace()`
- packaged graph objects such as `hmp.u01.gc.coarse`

Possible future helper additions to core `grip`:

- a small exported helper to build an edge matrix from an adjacency list
- a small exported helper to normalize coordinate matrices for comparison views

## Phased implementation plan

### Phase 1: data model and HMP reader

Files:

- `R/gripui-project.R`
- `R/gripui-catalog.R`
- `R/gripui-io.R`

Deliverables:

- `gripui_project()`
- `gripui_project_from_dir()`
- `gripui_validate_project()`
- HMP/U01 project loader
- tests for table normalization and manifest joins

### Phase 2: MVP app

Files:

- `R/gripui-app.R`
- `R/gripui-ui.R`
- `R/gripui-server.R`
- `inst/app/www/gripui.css`

Deliverables:

- stage filters
- catalog table
- single-layout viewer
- landscape plot
- status badges for interactive vs summary-only rows

### Phase 3: compare mode

Deliverables:

- pinned layouts
- shared camera
- optional Procrustes alignment
- metric deltas

### Phase 4: artifact export helpers

Deliverables:

- `gripui_export_layout_artifacts()`
- `gripui_write_bundle()`
- documentation for analysis scripts

This phase should be treated as the first canonical producer phase, not as an
optional cleanup pass. Once it lands, new projects should prefer the canonical
bundle contract over custom table/manifest merges.

## Testing plan

### Unit tests

- project validation
- path resolution
- catalog joins
- summary-only versus interactive row handling

### Integration tests

- app constructs from a tiny synthetic project
- app constructs from a `grip.compare.layouts(..., return.layouts = TRUE)` object
- HMP/U01 loader returns expected counts for rows and viewable layouts

### Manual checks

- start app on packaged HMP/U01 example
- verify that `coarse.search.088` and the repulsion sweep layouts render
- verify that summary-only rows remain selectable in the table but do not break
  the viewer

## Immediate implementation recommendation

Start with one narrow vertical slice:

1. implement `gripui_project_from_dir()` for the HMP/U01 output structure
2. build an MVP app that can:
   - show the catalog
   - filter by stage
   - view any row with `embedding.tsv`
3. once that works, add `gripui_project_from_compare()` so the same UI also
   works for smaller package examples

That path gets the real HMP/U01 use case working first while still leading to a
clean reusable design.
