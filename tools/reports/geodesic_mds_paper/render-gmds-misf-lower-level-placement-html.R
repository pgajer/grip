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

design_root <- file.path(repo_root, "output", "geodesic_mds_paper")
tmp_dir <- file.path(design_root, "tmp", "gmds-misf-lower-level-placement-2026-04-02")
interactive_dir <- file.path(design_root, "html")
output_html <- file.path(
  interactive_dir,
  "gmds_misf_lower_level_placement_3d_2026-04-02.html"
)
rds_path <- file.path(tmp_dir, "gmds_misf_lower_level_placement_results.rds")

dir.create(interactive_dir, recursive = TRUE, showWarnings = FALSE)

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

make_widget <- function(coords,
                        edges,
                        vertex_col = "#355070",
                        edge_col = "#adb5bd",
                        width = 420,
                        height = 320,
                        theta = 35,
                        phi = 24,
                        zoom = 0.82,
                        fov = 30) {
  coords <- as.matrix(coords)
  edges <- as.matrix(edges)

  keep <- apply(coords, 1L, function(x) all(is.finite(x)))
  good_edges <- if (nrow(edges)) keep[edges[, 1L]] & keep[edges[, 2L]] else logical(0L)
  edges <- if (nrow(edges)) edges[good_edges, , drop = FALSE] else edges

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
    rgl::lines3d(seg, color = grDevices::adjustcolor(edge_col, alpha.f = 0.78), lwd = 1.5)
  }
  if (any(keep)) {
    rgl::points3d(
      coords[keep, 1L],
      coords[keep, 2L],
      coords[keep, 3L],
      color = vertex_col,
      size = 6,
      alpha = 0.96
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

build_panel_title <- function(method_result) {
  row <- method_result$metrics[1L, , drop = FALSE]
  sprintf(
    "%s\nsigma %s, rho %s\ntotal %ss",
    row$method_label[[1L]],
    fmt_num(row$gmds_stress[[1L]], 4L),
    fmt_num(row$procrustes_rmse[[1L]], 4L),
    fmt_time(row$elapsed_sec[[1L]])
  )
}

stage_card_title <- function(method_label, snapshot) {
  sprintf("%s\n%s", method_label, snapshot$stage_label)
}

case_sections <- lapply(bundle$case_results, function(case_result) {
  ref_card <- card_div(
    title_text = "Reference surface",
    widget = make_widget(
      coords = case_result$case$truth,
      edges = case_result$case$edges,
      vertex_col = "#bc6c25"
    ),
    note = sprintf(
      "%s. Fixed top-level seed reused in this case: %s (%s).",
      case_result$case$label,
      case_result$top_seed$seed_label,
      case_result$top_seed$seed_kind
    )
  )

  final_cards <- c(
    list(ref_card),
    lapply(case_result$methods, function(method_result) {
      row <- method_result$metrics[1L, , drop = FALSE]
      card_div(
        title_text = build_panel_title(method_result),
        widget = make_widget(
          coords = method_result$display_coords,
          edges = case_result$case$edges,
          vertex_col = "#3a5a40"
        ),
        note = sprintf(
          "Top seed: %s (%s). t_place = %ss, t_ref = %ss, t_final = %ss.",
          row$top_seed_source[[1L]],
          row$top_seed_kind[[1L]],
          fmt_time(row$placement_elapsed_sec[[1L]]),
          fmt_time(row$refinement_elapsed_sec[[1L]]),
          fmt_time(row$final_elapsed_sec[[1L]])
        )
      )
    })
  )

  stage_sections <- lapply(case_result$methods, function(method_result) {
    row <- method_result$metrics[1L, , drop = FALSE]
    cards <- lapply(method_result$stage_display, function(snapshot) {
      card_div(
        title_text = stage_card_title(row$method_label[[1L]], snapshot),
        widget = make_widget(
          coords = snapshot$display_coords,
          edges = case_result$case$edges,
          vertex_col = "#355070"
        ),
        note = sprintf("Level %d, active n = %d.", snapshot$level, snapshot$active_n)
      )
    })
    htmltools::tags$div(
      class = "trace-block",
      htmltools::tags$h3(row$method_label[[1L]]),
      htmltools::tags$div(class = "panel-grid", cards)
    )
  })

  htmltools::tags$section(
    class = "case-section",
    htmltools::tags$h2(case_result$case$label),
    htmltools::tags$p(
      class = "case-lead",
      "All non-reference layouts are Procrustes-aligned to the reference paraboloid before rendering. The first gallery compares final layouts after full refinement and polish; the second shows the full lower-level placement trajectory for every method."
    ),
    htmltools::tags$h3("Final Layouts"),
    htmltools::tags$div(class = "panel-grid", final_cards),
    htmltools::tags$h3("Lower-Level Placement Trajectories"),
    stage_sections
  )
})

page <- htmltools::browsable(
  htmltools::tagList(
    htmltools::tags$head(
      htmltools::tags$title("Interactive 3D Panels: Phase D Lower-Level Placement"),
      htmltools::tags$style(htmltools::HTML(
        "
        body {
          margin: 0;
          font-family: 'Avenir Next', 'Segoe UI', sans-serif;
          background: #f7f5f2;
          color: #1f2933;
        }
        .page {
          max-width: 1500px;
          margin: 0 auto;
          padding: 28px 28px 56px;
        }
        h1, h2, h3 {
          margin: 0 0 12px;
          font-weight: 600;
          line-height: 1.2;
        }
        h1 { font-size: 30px; }
        h2 { font-size: 24px; margin-top: 34px; }
        h3 { font-size: 19px; margin-top: 22px; }
        p {
          line-height: 1.6;
          margin: 0 0 14px;
          max-width: 980px;
        }
        .case-section + .case-section {
          margin-top: 38px;
          padding-top: 24px;
          border-top: 1px solid #d7d0c8;
        }
        .panel-grid {
          display: grid;
          grid-template-columns: repeat(auto-fit, minmax(360px, 1fr));
          gap: 18px;
          align-items: start;
        }
        .panel-card {
          background: #ffffff;
          border: 1px solid #ddd6ce;
          border-radius: 16px;
          padding: 14px 14px 16px;
          box-shadow: 0 12px 30px rgba(31, 41, 51, 0.08);
        }
        .panel-card-title {
          margin-bottom: 10px;
        }
        .panel-card-line-main {
          font-size: 15px;
          font-weight: 600;
        }
        .panel-card-line-metric {
          font-size: 13px;
          color: #52606d;
          margin-top: 2px;
        }
        .panel-card-note {
          font-size: 13px;
          color: #52606d;
          margin-top: 10px;
          margin-bottom: 0;
        }
        .trace-block + .trace-block {
          margin-top: 18px;
        }
        "
      ))
    ),
    htmltools::tags$div(
      class = "page",
      htmltools::tags$h1("Phase D: MISF-GMDS Lower-Level Placement"),
      htmltools::tags$p(
        "This gallery fixes the top MISF scaffold to the best Phase C seed for each case, then compares how lower-level vertices are placed before the same sparse refinement and final pure-GMDS polish. It is the interactive companion to the Phase D report."
      ),
      case_sections
    )
  )
)

htmltools::save_html(
  page,
  file = output_html,
  libdir = paste0(tools::file_path_sans_ext(basename(output_html)), "_files")
)

message("Wrote interactive HTML: ", output_html)
