#!/usr/bin/env Rscript

repo_root <- normalizePath("/Users/pgajer/current_projects/grip", winslash = "/", mustWork = TRUE)
setwd(repo_root)

options(rgl.useNULL = TRUE)

benchmark_root <- file.path(repo_root, "output", "benchmarks", "weighted-grip-phase5-family-panel-2026-04-02")
output_dir <- file.path(repo_root, "output", "html")
output_html <- file.path(output_dir, "weighted_grip_phase5_layout_gallery_2026-04-02.html")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(repo_root, quiet = TRUE, helpers = FALSE)
} else if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(repo_root, quiet = TRUE, export_all = FALSE, helpers = FALSE)
} else if (requireNamespace("grip", quietly = TRUE)) {
  library(grip)
} else {
  stop("Install 'pkgload', 'devtools', or the 'grip' package to run this renderer.")
}

for (pkg in c("htmltools", "htmlwidgets", "rgl", "igraph")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required to render the gallery.", pkg))
  }
}

benchmark_cfg <- list(
  dims = 3L,
  lgkk_max_iter = 4L,
  gkk_max_iter = 4L,
  lgkk_local_nbrs = 8L,
  lgkk_landmark_count = 8L,
  weighted_core_lgkk_rounds = 1L,
  weighted_core_lgkk_scope = "all",
  weighted_core_lgkk_active_limit = 512L,
  weighted_lgkk_polish_rounds = 4L,
  stress_sample = 1000L,
  jitter_xy = 0.08,
  jitter_z = 0.08
)

`%||%` <- function(x, y) if (is.null(x)) y else x

family_configs <- list(
  list(
    id = "mesh",
    label = "Mesh saddle",
    preset = "mesh",
    builder = function() mesh.surface.graph(6, 6, surface = "saddle", amplitude = 0.75)
  ),
  list(
    id = "cylinder",
    label = "Cylinder hourglass",
    preset = "cylinder",
    builder = function() cylinder.surface.graph(8, 12, surface = "hourglass", amplitude = 0.30)
  ),
  list(
    id = "torus",
    label = "Torus pinched",
    preset = "torus",
    builder = function() torus.surface.graph(6, 6, surface = "pinched", amplitude = 0.22)
  ),
  list(
    id = "sphere",
    label = "Sphere wavy",
    preset = "sphere",
    builder = function() sphere.surface.graph(8, 10, surface = "wavy", amplitude = 0.18)
  ),
  list(
    id = "sierpinski_carpet",
    label = "Sierpinski carpet ripple",
    preset = "carpet",
    builder = function() sierpinski.carpet.surface.graph(
      level = 2,
      surface = "ripple",
      amplitude = 0.70,
      freq_u = 2,
      freq_v = 1
    )
  ),
  list(
    id = "cube_channel_network",
    label = "Cube channel network twisted",
    preset = "irregular",
    builder = function() cube.channel.network.surface.graph(
      level = 1,
      surface = "twisted",
      amplitude = 0.16,
      twist = 0.55
    )
  ),
  list(
    id = "irregular_annulus",
    label = "Irregular annulus folded",
    preset = "irregular",
    builder = function() irregular.annulus.surface.graph(
      rings = 6,
      outer_count = 24,
      surface = "folded",
      amplitude = 0.45
    )
  ),
  list(
    id = "irregular_torus",
    label = "Irregular torus pinched",
    preset = "irregular",
    builder = function() irregular.torus.surface.graph(
      major_rings = 6,
      tube_count = 10,
      surface = "pinched",
      amplitude = 0.18
    )
  )
)

family_configs <- family_configs[vapply(
  family_configs,
  function(cfg) cfg$id %in% c("mesh", "torus", "sierpinski_carpet", "cube_channel_network", "irregular_torus"),
  logical(1L)
)]

grip_ns <- asNamespace("grip")
pkg_internal <- function(name) get(name, envir = grip_ns)

fmt_num <- function(x, digits = 4L) {
  ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "NA")
}

normalize_coords <- function(coords) {
  coords <- as.matrix(coords)
  centered <- sweep(coords, 2L, colMeans(coords), "-", check.margin = FALSE)
  radius <- max(sqrt(rowSums(centered^2)))
  if (!is.finite(radius) || radius <= 0) {
    return(centered)
  }
  centered / radius
}

normalize_coords3d <- function(coords) {
  coords <- as.matrix(coords)
  center <- colMeans(coords)
  coords <- sweep(coords, 2L, center, "-", check.margin = FALSE)
  radius <- max(abs(coords))
  if (is.finite(radius) && radius > 0) {
    coords <- coords / radius
  }
  coords
}

base_param_coords <- function(spec) {
  coords <- as.matrix(spec$coords_param)
  if (ncol(coords) < 2L) {
    coords <- cbind(coords, rep(0, nrow(coords)))
  }
  coords[, 1:2, drop = FALSE]
}

build_initial_layout <- function(spec, dim, seed) {
  dim <- as.integer(dim)
  base2d <- normalize_coords(base_param_coords(spec))
  fam_index <- match(spec$family_id, vapply(family_configs, `[[`, "", "id"))
  draw_seed <- as.integer(100000L + 1000L * seed + 100L * dim + 10L * fam_index)
  set.seed(draw_seed)
  if (dim == 2L) {
    coords <- base2d + matrix(stats::rnorm(spec$n * 2L, sd = benchmark_cfg$jitter_xy), ncol = 2L)
  } else {
    coords <- cbind(base2d, 0)
    jitter <- matrix(stats::rnorm(spec$n * 3L, sd = benchmark_cfg$jitter_z), ncol = 3L)
    jitter[, 1:2] <- 0.45 * jitter[, 1:2, drop = FALSE]
    coords <- coords + jitter
  }
  storage.mode(coords) <- "double"
  sweep(coords, 2L, colMeans(coords), "-", check.margin = FALSE)
}

resolve_weighted_preset_args <- function(preset, dim) {
  resolver <- pkg_internal("grip.resolve.weighted.preset")
  resolver(
    preset = preset,
    dim = dim,
    placement = "barycenter",
    placement_missing = TRUE,
    rounds = 160L,
    rounds_missing = TRUE,
    final_rounds = 384L,
    final_rounds_missing = TRUE,
    num_init = 24L,
    num_init_missing = TRUE,
    num_nbrs = 20L,
    num_nbrs_missing = TRUE,
    r = 0.03,
    r_missing = TRUE,
    s = 7.5,
    s_missing = TRUE,
    repulsion_factor = 2.5,
    repulsion_factor_missing = TRUE
  )
}

prepare_family_spec <- function(cfg) {
  spec <- cfg$builder()
  spec$family_id <- cfg$id
  spec$family_label <- cfg$label
  spec$preset <- cfg$preset
  spec$graph <- igraph::graph_from_edgelist(as.matrix(spec$edges), directed = FALSE)
  spec$prepared_gkk <- grip.prepare.geodesic.kk(
    edges = spec$edges,
    n = spec$n,
    edge_weights = spec$edge_weights
  )
  spec$prepared_lgkk <- grip.prepare.landmark.geodesic.kk(
    edges = spec$edges,
    n = spec$n,
    edge_weights = spec$edge_weights,
    local_nbrs = benchmark_cfg$lgkk_local_nbrs,
    landmark_count = benchmark_cfg$lgkk_landmark_count
  )
  spec
}

score_layout_metrics <- function(spec, coords, runtime_sec) {
  gkk <- grip.score.geodesic.kk(coords = coords, prepared = spec$prepared_gkk)
  procrustes_rmse <- pkg_internal("grip.align.to.target.nd")(coords, spec$coords_surface)$rmse
  list(
    runtime_sec = runtime_sec,
    gkk_rel_rmse = gkk$gkk.weighted.rel.rmse[[1L]],
    procrustes_rmse = procrustes_rmse
  )
}

compute_layouts_3d <- function(spec, seed) {
  dim <- 3L
  tuning <- resolve_weighted_preset_args(spec$preset, dim = dim)
  initial <- build_initial_layout(spec, dim = dim, seed = seed)
  seed_method <- function(offset) as.integer(1000L * seed + 100L * dim + offset)

  methods <- list()
  methods$target <- list(coords = spec$coords_surface, metrics = NULL)
  methods$start <- list(coords = initial, metrics = score_layout_metrics(spec, initial, runtime_sec = NA_real_))

  timed <- system.time({
    coords_grip <- grip.layout.globalrep(
      edges = spec$edges,
      edge_weights = spec$edge_weights,
      n = spec$n,
      dim = dim,
      placement = tuning$placement,
      rounds = tuning$rounds,
      final_rounds = tuning$final_rounds,
      num_init = tuning$num_init,
      num_nbrs = tuning$num_nbrs,
      r = tuning$r,
      s = tuning$s,
      repulsion_factor = tuning$repulsion_factor,
      seed = seed_method(1L)
    )
  })
  methods$grip <- list(coords = coords_grip, metrics = score_layout_metrics(spec, coords_grip, unname(timed[["elapsed"]])))

  timed <- system.time({
    coords_wgrip <- grip.layout.globalrep.weighted(
      edges = spec$edges,
      edge_weights = spec$edge_weights,
      n = spec$n,
      dim = dim,
      preset = spec$preset,
      seed = seed_method(2L)
    )
  })
  methods$wgrip <- list(coords = coords_wgrip, metrics = score_layout_metrics(spec, coords_wgrip, unname(timed[["elapsed"]])))

  timed <- system.time({
    coords_wgrip_core_lgkk <- grip.layout.globalrep.weighted(
      edges = spec$edges,
      edge_weights = spec$edge_weights,
      n = spec$n,
      dim = dim,
      preset = spec$preset,
      lgkk_multiscale_rounds = benchmark_cfg$weighted_core_lgkk_rounds,
      lgkk_local_nbrs = benchmark_cfg$lgkk_local_nbrs,
      lgkk_landmark_count = benchmark_cfg$lgkk_landmark_count,
      lgkk_multiscale_scope = benchmark_cfg$weighted_core_lgkk_scope,
      lgkk_active_limit = benchmark_cfg$weighted_core_lgkk_active_limit,
      seed = seed_method(3L)
    )
  })
  methods$wgrip_core_lgkk <- list(coords = coords_wgrip_core_lgkk, metrics = score_layout_metrics(spec, coords_wgrip_core_lgkk, unname(timed[["elapsed"]])))

  timed <- system.time({
    coords_wgrip_polish_lgkk <- grip.layout.globalrep.weighted(
      edges = spec$edges,
      edge_weights = spec$edge_weights,
      n = spec$n,
      dim = dim,
      preset = spec$preset,
      lgkk_polish_rounds = benchmark_cfg$weighted_lgkk_polish_rounds,
      lgkk_local_nbrs = benchmark_cfg$lgkk_local_nbrs,
      lgkk_landmark_count = benchmark_cfg$lgkk_landmark_count,
      seed = seed_method(4L)
    )
  })
  methods$wgrip_polish_lgkk <- list(coords = coords_wgrip_polish_lgkk, metrics = score_layout_metrics(spec, coords_wgrip_polish_lgkk, unname(timed[["elapsed"]])))

  timed <- system.time({
    coords_kk <- igraph::layout_with_kk(
      spec$graph,
      coords = initial,
      dim = dim,
      weights = spec$edge_weights
    )
  })
  methods$kk <- list(coords = coords_kk, metrics = score_layout_metrics(spec, coords_kk, unname(timed[["elapsed"]])))

  timed <- system.time({
    coords_gkk <- grip.optimize.geodesic.kk(
      coords = coords_kk,
      prepared = spec$prepared_gkk,
      max_iter = benchmark_cfg$gkk_max_iter,
      scale_mode = "profiled",
      return_trace = TRUE
    )$coords
  })
  methods$gkk <- list(coords = coords_gkk, metrics = score_layout_metrics(spec, coords_gkk, unname(timed[["elapsed"]])))

  timed <- system.time({
    coords_lgkk <- grip.optimize.landmark.geodesic.kk(
      coords = coords_kk,
      prepared = spec$prepared_lgkk,
      max_iter = benchmark_cfg$lgkk_max_iter,
      return_trace = TRUE
    )$coords
  })
  methods$lgkk <- list(coords = coords_lgkk, metrics = score_layout_metrics(spec, coords_lgkk, unname(timed[["elapsed"]])))

  methods
}

edge_color_map <- function(values) {
  pal <- grDevices::colorRampPalette(c("#355070", "#588157", "#e9c46a", "#f4a261", "#e76f51"))(256)
  idx <- 1L + floor((values - min(values)) / max(1e-12, diff(range(values))) * 255)
  pal[pmax(1L, pmin(256L, idx))]
}

vertex_color_map <- function(coords) {
  pal <- grDevices::colorRampPalette(c("#264653", "#2a9d8f", "#e9c46a", "#f4a261"))(256)
  z <- coords[, 3L]
  idx <- 1L + floor((z - min(z)) / max(1e-12, diff(range(z))) * 255)
  pal[pmax(1L, pmin(256L, idx))]
}

scene_widget <- function(coords, edges, edge_weights, width = 320, height = 260, theta = 30, phi = 20) {
  coords <- normalize_coords3d(coords)
  edges <- as.matrix(edges)
  p1 <- coords[edges[, 1L], , drop = FALSE]
  p2 <- coords[edges[, 2L], , drop = FALSE]
  seg <- matrix(NA_real_, nrow = 3L * nrow(edges), ncol = 3L)
  seg[seq(1L, nrow(seg), by = 3L), ] <- p1
  seg[seq(2L, nrow(seg), by = 3L), ] <- p2

  edge_cols <- if (is.null(edge_weights)) rep("#7b8794", nrow(edges)) else edge_color_map(edge_weights)
  edge_cols <- rep(edge_cols, each = 3L)
  point_cols <- vertex_color_map(coords)
  point_size <- if (nrow(coords) <= 80L) 6 else if (nrow(coords) <= 180L) 5 else 4
  edge_lwd <- if (nrow(edges) <= 160L) 2 else 1.4

  rgl::open3d(useNULL = TRUE)
  rgl::bg3d(color = "#ffffff")
  rgl::material3d(specular = "#444444")
  rgl::lines3d(seg, color = edge_cols, lwd = edge_lwd)
  rgl::points3d(coords, color = point_cols, size = point_size, alpha = 0.95)
  rgl::aspect3d(1, 1, 1)
  rgl::view3d(theta = theta, phi = phi, zoom = 0.80, fov = 30)
  widget <- rgl::rglwidget(width = width, height = height)
  rgl::close3d()
  widget
}

metric_line <- function(metrics) {
  if (is.null(metrics)) {
    return(NULL)
  }
  sprintf(
    "runtime = %ss; GKK rel. RMSE = %s; Procrustes RMSE = %s",
    fmt_num(metrics$runtime_sec, 3L),
    fmt_num(metrics$gkk_rel_rmse, 4L),
    fmt_num(metrics$procrustes_rmse, 4L)
  )
}

method_label <- function(method) {
  switch(
    method,
    target = "Target geometry",
    start = "Shared start",
    grip = "GRIP",
    wgrip = "Weighted GRIP",
    wgrip_core_lgkk = "Weighted GRIP + core LGKK",
    wgrip_polish_lgkk = "Weighted GRIP + polish LGKK",
    kk = "KK",
    gkk = "KK->GKK",
    lgkk = "KK->LGKK",
    method
  )
}

gallery_card <- function(spec, method, coords, note = NULL) {
  htmltools::tags$div(
    class = "gallery-card",
    htmltools::tags$h4(method_label(method)),
    if (!is.null(note)) htmltools::tags$p(class = "gallery-metric", note),
    scene_widget(coords = coords, edges = spec$edges, edge_weights = spec$edge_weights)
  )
}

family_section <- function(cfg) {
  spec <- prepare_family_spec(cfg)
  layouts <- compute_layouts_3d(spec, seed = 1L)
  methods <- c("target", "start", "grip", "wgrip", "wgrip_core_lgkk", "wgrip_polish_lgkk", "kk", "gkk", "lgkk")
  cards <- lapply(methods, function(method) {
    gallery_card(
      spec = spec,
      method = method,
      coords = layouts[[method]]$coords,
      note = if (method == "target") {
        "Intrinsic 3D target geometry for the weighted family."
      } else {
        metric_line(layouts[[method]]$metrics)
      }
    )
  })

  htmltools::tags$section(
    class = "family-section",
    htmltools::tags$h2(spec$family_label),
    htmltools::tags$p(
      class = "family-meta",
      sprintf("preset = %s; n = %d; m = %d; family id = %s; representative benchmark seed = 1", spec$preset, spec$n, nrow(spec$edges), spec$family_id)
    ),
    htmltools::tags$div(class = "gallery-grid", cards)
  )
}

page <- htmltools::tagList(
  htmltools::tags$html(
    htmltools::tags$head(
      htmltools::tags$title("Weighted GRIP Phase 5 Layout Gallery"),
      htmltools::tags$style(htmltools::HTML('
body {
  font-family: "Avenir Next", "Helvetica Neue", Helvetica, Arial, sans-serif;
  margin: 0;
  color: #1f2933;
  background: linear-gradient(180deg, #f7f3ea 0%, #ffffff 24%);
}
.page {
  max-width: 1560px;
  margin: 0 auto;
  padding: 28px 24px 48px 24px;
}
h1 {
  margin-top: 0;
  color: #14213d;
}
.lede {
  max-width: 980px;
  font-size: 1.02rem;
  line-height: 1.55;
}
.family-section {
  margin-top: 2rem;
  padding-top: 1rem;
  border-top: 1px solid #d8dee9;
}
.family-meta {
  color: #52606d;
}
.gallery-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(340px, 1fr));
  gap: 18px;
  margin: 1rem 0 2rem 0;
}
.gallery-card {
  border: 1px solid #d9dee7;
  border-radius: 14px;
  padding: 14px;
  background: linear-gradient(180deg, #ffffff 0%, #f8fbff 100%);
  box-shadow: 0 8px 24px rgba(40, 52, 71, 0.08);
}
.gallery-card h4 {
  margin: 0 0 0.45rem 0;
  color: #1f2a44;
}
.gallery-metric {
  margin: 0 0 0.65rem 0;
  color: #52606d;
  font-size: 0.92rem;
}
.seed-details {
  margin-top: 1rem;
}
.seed-details summary {
  cursor: pointer;
  font-weight: 600;
  color: #183a5a;
  margin-bottom: 0.8rem;
}
      '))
    ),
    htmltools::tags$body(
      htmltools::tags$div(
        class = "page",
        htmltools::tags$h1("Weighted GRIP Phase 5: Interactive 3D Layout Gallery"),
        htmltools::tags$p(
          class = "lede",
          "This HTML companion covers the primary 3D layouts from the completed Phase 5 smoke-validation panel. ",
          "For each family, it shows the target geometry and every benchmarked layout method for the representative benchmark seed 1. ",
          "The PDF report carries the quantitative 2D and 3D figures and tables; this gallery is the visual 3D companion."
        ),
        lapply(family_configs, family_section)
      )
    )
  )
)

htmltools::save_html(page, file = output_html, background = "white")
message(sprintf("Wrote %s", output_html))
