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
  "paraboloid_gmds_bending_figure1_interactive_2026-03-31.html"
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

for (pkg in c("htmltools", "htmlwidgets", "rgl")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required to render the interactive HTML.", pkg))
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
    bending_stencils = grip:::grip.rectangular.grid.bending.stencils(side, side)
  )
}

run_fit <- function(case,
                    label,
                    ...) {
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
  score <- grip.score.geodesic.mds(
    fit$coords,
    prepared = case$prepared,
    anchor_coords = fit$anchor_coords,
    anchor_weight = fit$final_anchor_weight,
    bending_stencils = case$bending_stencils,
    bending_weight = fit$final_bending_weight
  )
  aligned <- grip:::grip.align.to.target.nd(
    fit$coords,
    case$truth,
    allow.reflection = TRUE
  )

  list(
    label = label,
    coords = fit$coords,
    aligned = aligned$aligned,
    elapsed_sec = elapsed,
    fit = fit,
    score = score,
    rmse = aligned$rmse
  )
}

build_panel_title <- function(kind, row = NULL) {
  if (identical(kind, "truth")) {
    return("Reference surface")
  }
  if (is.null(row) || !is.data.frame(row) || nrow(row) != 1L) {
    stop("Expected a one-row data frame for non-truth panel titles.")
  }

  label <- row$setting_label[[1L]]
  if (identical(row$setting_id[[1L]], "gmds_anchor_0.1")) {
    label <- "anchor only"
  }
  if (identical(row$setting_id[[1L]], "gmds_average")) {
    label <- "GMDS avg"
  }
  sprintf(
    "%s\nsigma %s, rho %s\neta %s, alpha_0.05 %s\nt %ss",
    label,
    fmt_num(row$gmds_stress[[1L]]),
    fmt_num(row$procrustes_rmse[[1L]]),
    fmt_num(row$roughness[[1L]]),
    fmt_num(row$area_q05_ratio[[1L]]),
    fmt_time(row$elapsed_sec[[1L]])
  )
}

make_widget <- function(coords,
                        edges,
                        vertex_col,
                        edge_col = "#adb5bd",
                        width = 430,
                        height = 330,
                        theta = 35,
                        phi = 24,
                        zoom = 0.82,
                        fov = 30) {
  coords <- as.matrix(coords)
  edges <- as.matrix(edges)

  keep <- apply(coords, 1L, function(x) all(is.finite(x)))
  good_edges <- keep[edges[, 1L]] & keep[edges[, 2L]]
  edges <- edges[good_edges, , drop = FALSE]

  seg <- matrix(NA_real_, nrow = 3L * nrow(edges), ncol = 3L)
  if (nrow(edges) > 0L) {
    seg[seq(1L, nrow(seg), by = 3L), ] <- coords[edges[, 1L], , drop = FALSE]
    seg[seq(2L, nrow(seg), by = 3L), ] <- coords[edges[, 2L], , drop = FALSE]
  }

  rgl::open3d(useNULL = TRUE)
  on.exit(try(rgl::close3d(), silent = TRUE), add = TRUE)
  rgl::bg3d(color = "#ffffff")
  rgl::material3d(specular = "#555555")
  rgl::light3d(theta = 35, phi = 20, viewpoint.rel = TRUE)
  rgl::light3d(theta = -55, phi = 25, viewpoint.rel = TRUE)

  if (nrow(edges) > 0L) {
    rgl::lines3d(seg, color = grDevices::adjustcolor(edge_col, alpha.f = 0.75), lwd = 1.5)
  }
  if (any(keep)) {
    rgl::points3d(
      coords[keep, 1L],
      coords[keep, 2L],
      coords[keep, 3L],
      color = vertex_col,
      size = 6,
      alpha = 0.95
    )
  }

  rgl::aspect3d(1, 1, 1)
  rgl::view3d(theta = theta, phi = phi, zoom = zoom, fov = fov)
  rgl::rglwidget(width = width, height = height)
}

card_div <- function(title_text, widget, note = NULL) {
  lines <- strsplit(title_text, "\n", fixed = TRUE)[[1L]]
  htmltools::tags$div(
    class = "panel-card",
    htmltools::tags$div(
      class = "panel-card-title",
      htmltools::tags$div(class = "panel-card-line panel-card-line-main", lines[[1L]]),
      if (length(lines) > 1L) {
        lapply(lines[-1L], function(x) {
          htmltools::tags$div(class = "panel-card-line panel-card-line-metric", x)
        })
      }
    ),
    htmltools::tags$div(class = "panel-widget", widget),
    if (!is.null(note)) htmltools::tags$p(class = "panel-card-note", note)
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

gmds_avg <- run_fit(case, "GMDS avg")
anchor_only <- run_fit(
  case,
  "anchor only",
  anchor_mode = "cmdscale",
  anchor_weight = 0.1,
  anchor_weight_end = 0.1,
  continuation = "constant"
)
bending_only <- run_fit(
  case,
  "beta = 1.000",
  bending_stencils = case$bending_stencils,
  bending_weight = 1,
  anchor_mode = "none"
)
hybrid <- run_fit(
  case,
  "lambda = 0.100, beta = 0.200",
  bending_stencils = case$bending_stencils,
  bending_weight = 0.2,
  anchor_mode = "cmdscale",
  anchor_weight = 0.1,
  anchor_weight_end = 0.1,
  continuation = "constant"
)

panels <- list(
  list(
    title = build_panel_title("truth"),
    coords = truth_display,
    vertex_col = "#bc6c25",
    note = "Reference paraboloid mesh used in Figure 1."
  ),
  list(
    title = build_panel_title("cmdscale", row_by_id("cmdscale")),
    coords = cmd_aligned,
    vertex_col = "#3a5a40",
    note = "Classical MDS initialization on the graph-geodesic distance matrix."
  ),
  list(
    title = build_panel_title("gmds_average", row_by_id("gmds_average")),
    coords = gmds_avg$aligned,
    vertex_col = "#3a5a40",
    note = "Untethered tie-averaged GMDS after 25 extra correction steps."
  ),
  list(
    title = build_panel_title("gmds_anchor_0.1", row_by_id("gmds_anchor_0.1")),
    coords = anchor_only$aligned,
    vertex_col = "#3a5a40",
    note = "GMDS with a fixed quadratic tether to the MDS embedding."
  ),
  list(
    title = build_panel_title("bend_beta_1", row_by_id("bend_beta_1")),
    coords = bending_only$aligned,
    vertex_col = "#3a5a40",
    note = "Best bending-only run from the selected Figure 1 comparison."
  ),
  list(
    title = build_panel_title("hybrid_lambda_0.1_beta_0.2", row_by_id("hybrid_lambda_0.1_beta_0.2")),
    coords = hybrid$aligned,
    vertex_col = "#3a5a40",
    note = "Best anchor-plus-bending run from the selected Figure 1 comparison."
  )
)

panel_cards <- lapply(panels, function(panel) {
  widget <- make_widget(
    coords = panel$coords,
    edges = case$edges,
    vertex_col = panel$vertex_col
  )
  card_div(panel$title, widget, note = panel$note)
})

page <- htmltools::browsable(
  htmltools::tagList(
    htmltools::tags$head(
      htmltools::tags$title("Interactive Figure 1: Paraboloid GMDS Bending Report"),
      htmltools::tags$style(htmltools::HTML(
        "
        body {
          margin: 0;
          background: #f5f1ea;
          color: #1f2933;
          font-family: 'Iowan Old Style', 'Palatino Linotype', 'Book Antiqua', Palatino, serif;
        }
        .page {
          max-width: 1500px;
          margin: 0 auto;
          padding: 28px 28px 36px;
        }
        .hero {
          background: linear-gradient(135deg, #fffdf9 0%, #f3ece1 100%);
          border: 1px solid #d8c9b7;
          border-radius: 18px;
          padding: 22px 24px;
          box-shadow: 0 10px 28px rgba(58, 36, 14, 0.08);
          margin-bottom: 22px;
        }
        .hero h1 {
          margin: 0 0 10px;
          font-size: 30px;
          line-height: 1.2;
          color: #6b3e26;
        }
        .hero p {
          margin: 10px 0 0;
          font-size: 17px;
          line-height: 1.55;
          max-width: 980px;
        }
        .hero code {
          background: rgba(107, 62, 38, 0.08);
          border-radius: 5px;
          padding: 1px 5px;
        }
        .panel-grid {
          display: grid;
          grid-template-columns: repeat(auto-fit, minmax(430px, 1fr));
          gap: 18px;
        }
        .panel-card {
          background: #fffdfb;
          border: 1px solid #d9d3c7;
          border-radius: 16px;
          box-shadow: 0 8px 22px rgba(51, 41, 28, 0.08);
          overflow: hidden;
        }
        .panel-card-title {
          padding: 16px 18px 8px;
          border-bottom: 1px solid #ede6da;
          background: linear-gradient(180deg, #fffdfa 0%, #faf5ed 100%);
          min-height: 92px;
        }
        .panel-card-line {
          margin: 0;
          white-space: pre-wrap;
        }
        .panel-card-line-main {
          font-size: 20px;
          font-weight: 700;
          color: #4a3728;
          margin-bottom: 4px;
        }
        .panel-card-line-metric {
          font-family: 'SFMono-Regular', Menlo, Consolas, monospace;
          font-size: 13px;
          color: #5b6570;
          line-height: 1.45;
        }
        .panel-widget {
          padding: 12px 14px 0;
        }
        .panel-card-note {
          margin: 0;
          padding: 12px 18px 16px;
          color: #5f5445;
          font-size: 14px;
          line-height: 1.45;
        }
        .footer-note {
          margin-top: 18px;
          padding: 14px 18px;
          background: #fff8ef;
          border: 1px solid #e4d7c3;
          border-radius: 14px;
          font-size: 14px;
          line-height: 1.5;
        }
        "
      ))
    ),
    htmltools::tags$div(
      class = "page",
      htmltools::tags$section(
        class = "hero",
        htmltools::tags$h1("Interactive Figure 1: Paraboloid GMDS Bending Report"),
        htmltools::tags$p(
          "This page reproduces the six layouts from Figure 1 of ",
          htmltools::tags$code("paraboloid_gmds_bending_report_2026-03-31.tex"),
          " as interactive ",
          htmltools::tags$code("rglwidget"),
          " panels. Drag to rotate, scroll to zoom, and right-drag to pan."
        ),
        htmltools::tags$p(
          "Each non-reference layout is Procrustes-aligned to the reference paraboloid exactly as in the static report, so the panels stay visually comparable. The titles reproduce the method labels and diagnostics used in the report."
        )
      ),
      htmltools::tags$div(class = "panel-grid", panel_cards),
      htmltools::tags$div(
        class = "footer-note",
        "Source settings were read from ",
        htmltools::tags$code(selected_csv),
        ". Geometry was recomputed deterministically from the same 12x12 paraboloid benchmark setup and then aligned to the reference surface before rendering."
      )
    )
  )
)

htmltools::save_html(
  html = page,
  file = output_html,
  libdir = paste0(tools::file_path_sans_ext(basename(output_html)), "_files"),
  background = "#f5f1ea"
)

cat(output_html, "\n")
