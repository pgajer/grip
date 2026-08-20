#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
if (length(file_arg) == 0L) {
  stop("Unable to determine script path from commandArgs().")
}
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."), winslash = "/", mustWork = TRUE)
setwd(repo_root)

options(rgl.useNULL = TRUE)

manual_root <- file.path(repo_root, "output", "geodesic_mds_paper")
tmp_dir <- file.path(manual_root, "tmp", "gmds-spring-repulsion-phase1-2026-04-01")
interactive_dir <- file.path(manual_root, "html")
output_html <- file.path(
  interactive_dir,
  "gmds_spring_repulsion_phase1_3d_panels_2026-04-01.html"
)
rds_path <- file.path(tmp_dir, "gmds_spring_repulsion_phase1_results.rds")

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

if (!file.exists(rds_path)) {
  stop("Expected benchmark bundle at: ", rds_path)
}

bundle <- readRDS(rds_path)

fmt_num <- function(x, digits = 4L) {
  ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "NA")
}

fmt_time <- function(x) {
  ifelse(
    is.finite(x),
    ifelse(x < 1, formatC(x, format = "f", digits = 3L), formatC(x, format = "f", digits = 2L)),
    "--"
  )
}

build_panel_title <- function(method_result) {
  row <- method_result$metrics[1L, , drop = FALSE]
  if (identical(row$method_id[[1L]], "reference")) {
    return(row$method_label[[1L]])
  }
  sprintf(
    "%s\nsigma %s, rho %s\nalpha_0.05 %s, t %ss",
    row$method_label[[1L]],
    fmt_num(row$gmds_stress[[1L]]),
    fmt_num(row$procrustes_rmse[[1L]]),
    fmt_num(row$area_q05_ratio[[1L]], 3L),
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

method_note <- function(method_id) {
  switch(
    method_id,
    reference = "Reference geometry used for side-by-side comparison.",
    cmdscale = "Classical MDS on the graph-geodesic distance matrix.",
    gmds_average = "Tie-averaged GMDS baseline from the compiled solver.",
    gmds_anchor = "GMDS with a linear MDS tether schedule.",
    gmds_spring_repulsion_anchor = "Anchor-GMDS plus the new generic edge-spring and graph-aware repulsion terms.",
    NULL
  )
}

case_sections <- lapply(bundle$case_results, function(case_result) {
  case <- case_result$case
  if (case$dim != 3L) {
    return(NULL)
  }
  panel_cards <- lapply(case_result$methods, function(method_result) {
    widget <- make_widget(
      coords = method_result$display_coords,
      edges = case$edges,
      vertex_col = if (identical(method_result$metrics$method_id[[1L]], "reference")) "#bc6c25" else "#3a5a40"
    )
    card_div(
      title_text = build_panel_title(method_result),
      widget = widget,
      note = method_note(method_result$metrics$method_id[[1L]])
    )
  })

  htmltools::tags$section(
    class = "case-section",
    htmltools::tags$h2(case$label),
    htmltools::tags$p(
      class = "case-lead",
      sprintf(
        "All non-reference layouts are rigidly aligned to the reference %s before rendering. The alignment changes only translation, uniform scale, rotation, and optional reflection.",
        if (identical(case$surface, "saddle")) "surface" else tolower(case$surface)
      )
    ),
    htmltools::tags$div(class = "panel-grid", panel_cards)
  )
})
case_sections <- Filter(Negate(is.null), case_sections)

page <- htmltools::browsable(
  htmltools::tagList(
    htmltools::tags$head(
      htmltools::tags$title("Interactive 3D Panels: GMDS Spring-Repulsion Phase 1"),
      htmltools::tags$style(htmltools::HTML(
        "
        body {
          margin: 0;
          background: #f5f1ea;
          color: #1f2933;
          font-family: 'Iowan Old Style', 'Palatino Linotype', 'Book Antiqua', Palatino, serif;
        }
        .page {
          max-width: 1560px;
          margin: 0 auto;
          padding: 28px 28px 36px;
        }
        .hero {
          background: linear-gradient(135deg, #fffdf9 0%, #f3ece1 100%);
          border: 1px solid #d8c9b7;
          border-radius: 18px;
          padding: 22px 24px;
          box-shadow: 0 10px 28px rgba(58, 36, 14, 0.08);
          margin-bottom: 24px;
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
          max-width: 1080px;
        }
        .hero code {
          background: rgba(107, 62, 38, 0.08);
          border-radius: 5px;
          padding: 1px 5px;
        }
        .case-section {
          margin-top: 26px;
        }
        .case-section h2 {
          margin: 0 0 8px;
          font-size: 28px;
          color: #4a3728;
        }
        .case-lead {
          margin: 0 0 14px;
          color: #5f5445;
          font-size: 15px;
          line-height: 1.5;
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
          margin-top: 22px;
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
        htmltools::tags$h1("Interactive 3D Panels: GMDS Spring-Repulsion Phase 1"),
        htmltools::tags$p(
          "This page shows every 3D layout included in the spring-repulsion Phase 1 report as an interactive ",
          htmltools::tags$code("rglwidget"),
          " panel. Drag to rotate, scroll to zoom, and right-drag to pan."
        ),
        htmltools::tags$p(
          "The panels are loaded from the saved benchmark object at ",
          htmltools::tags$code(rds_path),
          ", so they match the report figures rather than recomputing a potentially different run."
        )
      ),
      case_sections,
      htmltools::tags$div(
        class = "footer-note",
        "The displayed layouts are the same benchmark outputs used in the LaTeX report. Non-reference panels are aligned only by translation, uniform scaling, rotation, and optional reflection so that the shapes are comparable; the alignment does not bend or smooth the mesh."
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
