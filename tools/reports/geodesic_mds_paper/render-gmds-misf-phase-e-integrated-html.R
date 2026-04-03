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

design_root <- file.path(repo_root, "dev", "design")
tmp_dir <- file.path(design_root, "tmp", "gmds-misf-phase-e-integrated-2026-04-02")
interactive_dir <- file.path(design_root, "interactive-prototypes")
output_html <- file.path(interactive_dir, "gmds_misf_phase_e_integrated_3d_2026-04-02.html")
rds_path <- file.path(tmp_dir, "gmds_misf_phase_e_integrated_results.rds")

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

make_widget <- function(coords,
                        edges,
                        vertex_col = "#355070",
                        edge_col = "#adb5bd",
                        width = 400,
                        height = 300,
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
    rgl::points3d(coords[keep, 1L], coords[keep, 2L], coords[keep, 3L], color = vertex_col, size = 6, alpha = 0.96)
  }
  rgl::aspect3d(1, 1, 1)
  rgl::view3d(theta = theta, phi = phi, zoom = zoom, fov = fov)
  rgl::rglwidget(width = width, height = height)
}

card_div <- function(title_text, widget, note = NULL) {
  htmltools::tags$div(
    class = "panel-card",
    htmltools::tags$div(class = "panel-card-title", title_text),
    htmltools::tags$div(class = "panel-widget", widget),
    if (!is.null(note)) htmltools::tags$p(class = "panel-card-note", note)
  )
}

stage_label <- function(name) {
  switch(
    name,
    top_level = "top level",
    after_insertion = "after insertion",
    after_refinement = "after refinement",
    final = "final",
    name
  )
}

case_sections <- lapply(bundle$case_results, function(case_result) {
  candidate_cards <- lapply(case_result$candidates, function(candidate) {
    look <- candidate$lookahead_metrics[1L, , drop = FALSE]
    card_div(
      title_text = sprintf(
        "%s\nlook sigma %s, rho %s\nfinal sigma %s, rho %s",
        candidate$method_label,
        fmt_num(look$lookahead_sigma[[1L]], 4L),
        fmt_num(look$lookahead_rho[[1L]], 4L),
        fmt_num(candidate$final_metrics$gmds_stress[[1L]], 4L),
        fmt_num(candidate$final_metrics$procrustes_rmse[[1L]], 4L)
      ),
      widget = make_widget(candidate$final_display, case_result$case$edges),
      note = "Candidate final pipeline"
    )
  })

  top_cards <- c(
    list(card_div("Reference surface", make_widget(case_result$case$truth, case_result$case$edges, vertex_col = "#bc6c25"))),
    lapply(case_result$candidates, function(candidate) {
      card_div(candidate$method_label, make_widget(candidate$top_display, case_result$case$edges), "Corrected top-level seed")
    })
  )

  lookahead_cards <- c(
    list(card_div("Reference surface", make_widget(case_result$case$truth, case_result$case$edges, vertex_col = "#bc6c25"))),
    lapply(case_result$candidates, function(candidate) {
      look <- candidate$lookahead_metrics[1L, , drop = FALSE]
      card_div(
        sprintf("%s\nsigma %s, rho %s", candidate$method_label, fmt_num(look$lookahead_sigma[[1L]], 4L), fmt_num(look$lookahead_rho[[1L]], 4L)),
        make_widget(candidate$lookahead_display, case_result$case$edges),
        "One-level weighted-KK lookahead"
      )
    })
  )

  comparison_cards <- lapply(case_result$comparison_entries, function(entry) {
    card_div(
      sprintf("%s\nsigma %s, rho %s", entry$method_label, fmt_num(entry$gmds_stress, 4L), fmt_num(entry$procrustes_rmse, 4L)),
      make_widget(entry$display_coords, case_result$case$edges, vertex_col = if (entry$method_id == "reference") "#bc6c25" else "#355070"),
      "Final comparison"
    )
  })

  stage_sections <- lapply(case_result$selected_stage_methods, function(method) {
    cards <- lapply(names(method$stage_display), function(stage_name) {
      card_div(
        sprintf("%s\n%s", method$method_label, stage_label(stage_name)),
        make_widget(method$stage_display[[stage_name]], case_result$case$edges),
        NULL
      )
    })
    htmltools::tags$div(
      class = "trace-block",
      htmltools::tags$h3(method$method_label),
      htmltools::tags$div(class = "panel-grid", cards)
    )
  })

  htmltools::tags$section(
    class = "case-section",
    htmltools::tags$h2(case_result$case$label),
    htmltools::tags$p(
      class = "case-lead",
      "All non-reference layouts are Procrustes-aligned to the reference paraboloid before rendering. The sections below show the corrected top-level seeds, the one-level weighted-KK lookahead used for candidate ranking, the final integrated pipelines, and stagewise views for the selected pipelines."
    ),
    htmltools::tags$h3("Corrected Top-Level Seeds"),
    htmltools::tags$div(class = "panel-grid", top_cards),
    htmltools::tags$h3("Weighted-KK Lookahead"),
    htmltools::tags$div(class = "panel-grid", lookahead_cards),
    htmltools::tags$h3("Final Candidate Pipelines"),
    htmltools::tags$div(class = "panel-grid", candidate_cards),
    htmltools::tags$h3("Final Comparison"),
    htmltools::tags$div(class = "panel-grid", comparison_cards),
    if (length(stage_sections)) {
      htmltools::tagList(
        htmltools::tags$h3("Selected Stage Layouts"),
        stage_sections
      )
    }
  )
})

page <- htmltools::browsable(
  htmltools::tagList(
    htmltools::tags$head(
      htmltools::tags$title("Interactive 3D Panels: Phase E MISF-GMDS"),
      htmltools::tags$style(htmltools::HTML(
        "
        body { margin: 0; background: #f5f1ea; color: #1f2933; font-family: 'Iowan Old Style', 'Palatino Linotype', 'Book Antiqua', Palatino, serif; }
        .page { max-width: 1640px; margin: 0 auto; padding: 28px 28px 36px; }
        .hero { background: linear-gradient(135deg, #fffdf9 0%, #f3ece1 100%); border: 1px solid #d8c9b7; border-radius: 18px; padding: 22px 24px; box-shadow: 0 10px 28px rgba(58, 36, 14, 0.08); margin-bottom: 24px; }
        .hero h1 { margin: 0 0 10px; font-size: 30px; line-height: 1.2; color: #6b3e26; }
        .hero p { margin: 10px 0 0; font-size: 17px; line-height: 1.55; max-width: 1120px; }
        .case-section { margin-top: 26px; }
        .case-lead { max-width: 1200px; line-height: 1.55; }
        .panel-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(340px, 1fr)); gap: 16px; margin: 14px 0 22px; }
        .panel-card { background: #fffdfa; border: 1px solid #d9ccb9; border-radius: 16px; padding: 12px 12px 14px; box-shadow: 0 10px 24px rgba(58, 36, 14, 0.07); }
        .panel-card-title { white-space: pre-line; font-size: 15px; line-height: 1.35; font-weight: 600; color: #5e3b24; margin-bottom: 8px; }
        .panel-card-note { font-size: 13px; line-height: 1.45; margin: 8px 0 0; color: #5c6670; }
        .trace-block { margin-top: 12px; }
        h2, h3 { color: #5e3b24; }
        "
      ))
    ),
    htmltools::tags$div(
      class = "page",
      htmltools::tags$section(
        class = "hero",
        htmltools::tags$h1("Phase E: Integrated MISF-GMDS Seed Selection"),
        htmltools::tags$p(
          "This gallery compares corrected top-level seed candidates, their one-level weighted-KK lookahead layouts, their final full weighted-KK pipelines, and the selected stagewise trajectories. The intention is to see whether a simple lookahead-based seed selector can improve on the fixed-top pipelines from Phases C and D."
        )
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
message("Wrote Phase E interactive HTML: ", output_html)
