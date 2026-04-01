#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
if (length(file_arg) == 0L) {
  stop("Unable to determine script path from commandArgs().")
}
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
setwd(repo_root)

options(rgl.useNULL = TRUE)

manual_root <- file.path(repo_root, "dev", "manual")
interactive_dir <- file.path(manual_root, "interactive-prototypes")
tmp_dir <- file.path(manual_root, "tmp", "paraboloid-gmds-bending-2026-03-31")
output_html <- file.path(
  interactive_dir,
  "paraboloid_gmds_bending_figure1_diagnostics_2026-03-31.html"
)

dir.create(interactive_dir, recursive = TRUE, showWarnings = FALSE)

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(repo_root, quiet = TRUE, helpers = FALSE)
} else if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(repo_root, quiet = TRUE, export_all = FALSE, helpers = FALSE)
} else if (requireNamespace("grip", quietly = TRUE)) {
  library(grip)
} else {
  stop("Install 'pkgload', 'devtools', or the 'grip' package to run this renderer.")
}

for (pkg in c("htmltools", "rgl")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required to render the diagnostic HTML.", pkg))
  }
}

fmt_num <- function(x, digits = 4L) {
  ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "NA")
}

fmt_time <- function(x) {
  ifelse(
    is.finite(x),
    ifelse(
      x < 1,
      formatC(x, format = "f", digits = 3L),
      formatC(x, format = "f", digits = 2L)
    ),
    "NA"
  )
}

grid_mesh_triangles <- function(h, w) {
  index <- matrix(seq_len(h * w), nrow = h, ncol = w, byrow = TRUE)
  triangles <- vector("list", 2L * (h - 1L) * (w - 1L))
  k <- 1L
  for (r in seq_len(h - 1L)) {
    for (c in seq_len(w - 1L)) {
      triangles[[k]] <- c(index[r, c], index[r + 1L, c], index[r, c + 1L])
      k <- k + 1L
      triangles[[k]] <- c(index[r + 1L, c], index[r + 1L, c + 1L], index[r, c + 1L])
      k <- k + 1L
    }
  }
  do.call(rbind, triangles)
}

build_case <- function(side = 12L, iter_budget = 25L) {
  bundle <- mesh.surface.graph(
    side,
    side,
    surface = "paraboloid",
    amplitude = 0.35,
    connectivity = "orthogonal",
    normalize = "median"
  )
  prepared <- grip.prepare.geodesic.kk(
    edges = bundle$edges,
    n = bundle$n,
    edge_weights = bundle$edge_weights,
    tie_mode = "average"
  )
  started <- proc.time()[["elapsed"]]
  cmd <- grip:::grip.classical.mds.embedding(prepared, dim = 3L, eig = TRUE)
  cmd_elapsed <- proc.time()[["elapsed"]] - started

  list(
    side = as.integer(side),
    iter_budget = as.integer(iter_budget),
    bundle = bundle,
    truth = bundle$coords_surface,
    edges = bundle$edges,
    prepared = prepared,
    cmd = cmd,
    cmd_elapsed = cmd_elapsed,
    triangles = grid_mesh_triangles(side, side),
    bending_stencils = grip:::grip.rectangular.grid.bending.stencils(side, side)
  )
}

run_fit <- function(case, ...) {
  started <- proc.time()[["elapsed"]]
  fit <- grip.optimize.geodesic.mds(
    coords = case$cmd$coords,
    prepared = case$prepared,
    max_iter = case$iter_budget,
    engine = "cpp",
    return_trace = FALSE,
    n_threads = 1L,
    ...
  )
  elapsed <- proc.time()[["elapsed"]] - started
  aligned <- grip:::grip.align.to.target.nd(
    fit$coords,
    case$truth,
    allow.reflection = TRUE
  )
  list(
    fit = fit,
    elapsed_sec = elapsed,
    aligned = aligned$aligned,
    rmse = aligned$rmse
  )
}

make_widget <- function(coords,
                        edges,
                        triangles,
                        vertex_colors,
                        face_color = "#cbd5e1",
                        face_alpha = 0.34,
                        highlight_idx = integer(0),
                        label_text = character(0),
                        width = 360,
                        height = 280,
                        theta = 35,
                        phi = 24,
                        zoom = 0.84,
                        fov = 30) {
  coords <- as.matrix(coords)
  edges <- as.matrix(edges)
  triangles <- as.matrix(triangles)

  tri_xyz <- coords[as.vector(t(triangles)), , drop = FALSE]
  edge_xyz <- matrix(NA_real_, nrow = 3L * nrow(edges), ncol = 3L)
  edge_xyz[seq(1L, nrow(edge_xyz), by = 3L), ] <- coords[edges[, 1L], , drop = FALSE]
  edge_xyz[seq(2L, nrow(edge_xyz), by = 3L), ] <- coords[edges[, 2L], , drop = FALSE]

  rgl::open3d(useNULL = TRUE)
  on.exit(try(rgl::close3d(), silent = TRUE), add = TRUE)
  rgl::bg3d(color = "#ffffff")
  rgl::material3d(specular = "#525252")
  rgl::light3d(theta = 35, phi = 20, viewpoint.rel = TRUE)
  rgl::light3d(theta = -60, phi = 28, viewpoint.rel = TRUE)

  rgl::triangles3d(
    x = tri_xyz[, 1L],
    y = tri_xyz[, 2L],
    z = tri_xyz[, 3L],
    color = grDevices::adjustcolor(face_color, alpha.f = face_alpha)
  )
  rgl::lines3d(
    edge_xyz,
    color = grDevices::adjustcolor("#7b8794", alpha.f = 0.75),
    lwd = 1.3
  )
  rgl::points3d(
    coords[, 1L],
    coords[, 2L],
    coords[, 3L],
    color = vertex_colors,
    size = 6,
    alpha = 0.96
  )

  if (length(highlight_idx) > 0L) {
    rgl::points3d(
      coords[highlight_idx, 1L],
      coords[highlight_idx, 2L],
      coords[highlight_idx, 3L],
      color = "#b91c1c",
      size = 9,
      alpha = 1
    )
    rgl::text3d(
      coords[highlight_idx, 1L],
      coords[highlight_idx, 2L],
      coords[highlight_idx, 3L],
      texts = label_text,
      color = "#111827",
      cex = 0.9,
      adj = c(0.5, -0.2)
    )
  }

  rgl::aspect3d(1, 1, 1)
  rgl::view3d(theta = theta, phi = phi, zoom = zoom, fov = fov)
  rgl::rglwidget(width = width, height = height)
}

diagnostic_card <- function(title_text, subtitle_text, widget) {
  htmltools::tags$div(
    class = "diag-card",
    htmltools::tags$div(
      class = "diag-card-head",
      htmltools::tags$div(class = "diag-card-title", title_text),
      htmltools::tags$div(class = "diag-card-subtitle", subtitle_text)
    ),
    htmltools::tags$div(class = "diag-widget", widget)
  )
}

method_section <- function(method_title, metrics_line, cards) {
  htmltools::tags$section(
    class = "method-section",
    htmltools::tags$div(
      class = "method-head",
      htmltools::tags$h2(method_title),
      htmltools::tags$p(metrics_line)
    ),
    htmltools::tags$div(class = "diag-grid", cards)
  )
}

selected_csv <- file.path(tmp_dir, "paraboloid_bending_selected.csv")
if (!file.exists(selected_csv)) {
  stop("Expected benchmark metadata at: ", selected_csv)
}
selected_meta <- utils::read.csv(selected_csv, stringsAsFactors = FALSE)

row_by_id <- function(id) {
  hit <- selected_meta[selected_meta$setting_id == id, , drop = FALSE]
  if (nrow(hit) != 1L) {
    stop("Could not find a unique selected row for setting id: ", id)
  }
  hit
}

case <- build_case(side = 12L, iter_budget = 25L)
truth_display <- grip:::grip.normalize.coords(case$truth)
cmd_aligned <- grip:::grip.align.to.target.nd(
  case$cmd$coords,
  case$truth,
  allow.reflection = TRUE
)$aligned

gmds_avg <- run_fit(case)
anchor_only <- run_fit(
  case,
  anchor_mode = "cmdscale",
  anchor_weight = 0.1,
  anchor_weight_end = 0.1,
  continuation = "constant"
)
bending_only <- run_fit(
  case,
  bending_stencils = case$bending_stencils,
  bending_weight = 1,
  anchor_mode = "none"
)
hybrid <- run_fit(
  case,
  bending_stencils = case$bending_stencils,
  bending_weight = 0.2,
  anchor_mode = "cmdscale",
  anchor_weight = 0.1,
  anchor_weight_end = 0.1,
  continuation = "constant"
)

index <- matrix(seq_len(case$side * case$side), nrow = case$side, ncol = case$side, byrow = TRUE)
row_idx <- rep(seq_len(case$side), each = case$side)
col_idx <- rep(seq_len(case$side), times = case$side)
row_palette <- grDevices::hcl.colors(case$side, palette = "Temps")
col_palette <- grDevices::hcl.colors(case$side, palette = "YlGnBu")
row_colors <- row_palette[row_idx]
col_colors <- col_palette[col_idx]
center_cells <- expand.grid(r = 5:8, c = 5:8, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
highlight_idx <- index[as.matrix(center_cells)]
highlight_labels <- sprintf("%d,%d", center_cells$r, center_cells$c)

methods <- list(
  list(
    id = "truth",
    title = "Reference surface",
    metrics = "Ground-truth paraboloid mesh. Row and column colors should appear as clean, non-overlapping bands.",
    coords = truth_display
  ),
  list(
    id = "cmdscale",
    title = "cmdscale",
    metrics = with(row_by_id("cmdscale"), sprintf(
      "sigma %s, rho %s, eta %s, alpha_0.05 %s, t %ss",
      fmt_num(gmds_stress[[1L]]),
      fmt_num(procrustes_rmse[[1L]]),
      fmt_num(roughness[[1L]]),
      fmt_num(area_q05_ratio[[1L]]),
      fmt_time(elapsed_sec[[1L]])
    )),
    coords = cmd_aligned
  ),
  list(
    id = "gmds_average",
    title = "GMDS avg",
    metrics = with(row_by_id("gmds_average"), sprintf(
      "sigma %s, rho %s, eta %s, alpha_0.05 %s, t %ss",
      fmt_num(gmds_stress[[1L]]),
      fmt_num(procrustes_rmse[[1L]]),
      fmt_num(roughness[[1L]]),
      fmt_num(area_q05_ratio[[1L]]),
      fmt_time(elapsed_sec[[1L]])
    )),
    coords = gmds_avg$aligned
  ),
  list(
    id = "gmds_anchor_0.1",
    title = "anchor only",
    metrics = with(row_by_id("gmds_anchor_0.1"), sprintf(
      "sigma %s, rho %s, eta %s, alpha_0.05 %s, t %ss",
      fmt_num(gmds_stress[[1L]]),
      fmt_num(procrustes_rmse[[1L]]),
      fmt_num(roughness[[1L]]),
      fmt_num(area_q05_ratio[[1L]]),
      fmt_time(elapsed_sec[[1L]])
    )),
    coords = anchor_only$aligned
  ),
  list(
    id = "bend_beta_1",
    title = "beta = 1.000",
    metrics = with(row_by_id("bend_beta_1"), sprintf(
      "sigma %s, rho %s, eta %s, alpha_0.05 %s, t %ss",
      fmt_num(gmds_stress[[1L]]),
      fmt_num(procrustes_rmse[[1L]]),
      fmt_num(roughness[[1L]]),
      fmt_num(area_q05_ratio[[1L]]),
      fmt_time(elapsed_sec[[1L]])
    )),
    coords = bending_only$aligned
  ),
  list(
    id = "hybrid_lambda_0.1_beta_0.2",
    title = "lambda = 0.100, beta = 0.200",
    metrics = with(row_by_id("hybrid_lambda_0.1_beta_0.2"), sprintf(
      "sigma %s, rho %s, eta %s, alpha_0.05 %s, t %ss",
      fmt_num(gmds_stress[[1L]]),
      fmt_num(procrustes_rmse[[1L]]),
      fmt_num(roughness[[1L]]),
      fmt_num(area_q05_ratio[[1L]]),
      fmt_time(elapsed_sec[[1L]])
    )),
    coords = hybrid$aligned
  )
)

method_sections <- lapply(methods, function(method) {
  cards <- list(
    diagnostic_card(
      "Rows colored",
      "If topology is intact, row bands should stay identifiable even when they overlap in 3D.",
      make_widget(
        coords = method$coords,
        edges = case$edges,
        triangles = case$triangles,
        vertex_colors = row_colors,
        face_color = "#d6ccc2",
        face_alpha = 0.36
      )
    ),
    diagnostic_card(
      "Columns colored",
      "Column colors help show whether the apparent middle spine is really multiple folded strips from different columns.",
      make_widget(
        coords = method$coords,
        edges = case$edges,
        triangles = case$triangles,
        vertex_colors = col_colors,
        face_color = "#d8e2dc",
        face_alpha = 0.36
      )
    ),
    diagnostic_card(
      "Center labels",
      "The red points mark the 4x4 center block, labeled by grid coordinates (row,column).",
      make_widget(
        coords = method$coords,
        edges = case$edges,
        triangles = case$triangles,
        vertex_colors = rep("#475569", length(row_colors)),
        face_color = "#cbd5e1",
        face_alpha = 0.28,
        highlight_idx = highlight_idx,
        label_text = highlight_labels
      )
    )
  )
  method_section(method$title, method$metrics, cards)
})

page <- htmltools::browsable(
  htmltools::tagList(
    htmltools::tags$head(
      htmltools::tags$title("Paraboloid Figure 1 Diagnostics"),
      htmltools::tags$style(htmltools::HTML(
        "
        body {
          margin: 0;
          background: #f4efe8;
          color: #1f2933;
          font-family: 'Iowan Old Style', 'Palatino Linotype', 'Book Antiqua', Palatino, serif;
        }
        .page {
          max-width: 1540px;
          margin: 0 auto;
          padding: 28px 28px 42px;
        }
        .hero {
          background: linear-gradient(135deg, #fffdfa 0%, #f6ecdf 100%);
          border: 1px solid #dcc8ad;
          border-radius: 18px;
          box-shadow: 0 10px 28px rgba(72, 43, 18, 0.08);
          padding: 22px 24px;
          margin-bottom: 22px;
        }
        .hero h1 {
          margin: 0 0 10px;
          font-size: 31px;
          line-height: 1.2;
          color: #6d4024;
        }
        .hero p {
          margin: 10px 0 0;
          font-size: 17px;
          line-height: 1.55;
          max-width: 1100px;
        }
        .hero code {
          background: rgba(109, 64, 36, 0.08);
          border-radius: 5px;
          padding: 1px 5px;
        }
        .method-section {
          margin-top: 20px;
          padding: 18px 18px 20px;
          background: #fffdfb;
          border: 1px solid #ddd3c7;
          border-radius: 18px;
          box-shadow: 0 8px 22px rgba(51, 41, 28, 0.07);
        }
        .method-head h2 {
          margin: 0;
          font-size: 24px;
          color: #4b3524;
        }
        .method-head p {
          margin: 8px 0 0;
          font-family: 'SFMono-Regular', Menlo, Consolas, monospace;
          font-size: 13px;
          color: #59636e;
          line-height: 1.5;
        }
        .diag-grid {
          display: grid;
          grid-template-columns: repeat(auto-fit, minmax(360px, 1fr));
          gap: 16px;
          margin-top: 16px;
        }
        .diag-card {
          background: #fffdfa;
          border: 1px solid #ebe2d5;
          border-radius: 16px;
          overflow: hidden;
        }
        .diag-card-head {
          min-height: 74px;
          padding: 14px 16px 10px;
          border-bottom: 1px solid #efe7db;
          background: linear-gradient(180deg, #fffdfa 0%, #faf5ee 100%);
        }
        .diag-card-title {
          font-size: 19px;
          font-weight: 700;
          color: #4c392b;
          margin-bottom: 4px;
        }
        .diag-card-subtitle {
          font-size: 13px;
          line-height: 1.45;
          color: #5e5c59;
        }
        .diag-widget {
          padding: 10px 12px 6px;
        }
        .footer-note {
          margin-top: 20px;
          padding: 14px 18px;
          background: #fff8ef;
          border: 1px solid #e3d7c3;
          border-radius: 14px;
          font-size: 14px;
          line-height: 1.55;
        }
        "
      ))
    ),
    htmltools::tags$div(
      class = "page",
      htmltools::tags$section(
        class = "hero",
        htmltools::tags$h1("Paraboloid Figure 1 Diagnostics"),
        htmltools::tags$p(
          "This page is meant to answer one narrow question: are the strange central structures in the GMDS panels coming from a bad graph topology, or from a valid rectangular mesh that has folded over itself in 3D?"
        ),
        htmltools::tags$p(
          "The diagnostics are all built from the same six methods used in Figure 1 of ",
          htmltools::tags$code("paraboloid_gmds_bending_report_2026-03-31.tex"),
          ". Each section shows the same embedding three ways: vertices colored by row index, vertices colored by column index, and the center 4x4 block labeled by grid coordinates."
        ),
        htmltools::tags$p(
          "The translucent faces are especially useful here: if the graph topology were wrong, the row and column color bands would break combinatorially. If the topology is intact but the embedding has folded over, the same row and column bands remain coherent while multiple strips pile up spatially in the center."
        )
      ),
      method_sections,
      htmltools::tags$div(
        class = "footer-note",
        "The underlying mesh graph remains the standard 12x12 orthogonal rectangular grid. The non-reference layouts were recomputed from the same benchmark setup, then globally aligned to the reference surface for comparison. The alignment is only a rotation/reflection/rescaling and does not alter adjacency or create the central fold."
      )
    )
  )
)

htmltools::save_html(
  html = page,
  file = output_html,
  libdir = paste0(tools::file_path_sans_ext(basename(output_html)), "_files"),
  background = "#f4efe8"
)

cat(output_html, "\n")
