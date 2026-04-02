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

design_root <- file.path(repo_root, "dev", "design")
tmp_dir <- file.path(design_root, "tmp", "gmds-misf-seeded-pipeline-2026-04-02")
interactive_dir <- file.path(design_root, "interactive-prototypes")
output_html <- file.path(
  interactive_dir,
  "gmds_misf_seeded_pipeline_3d_2026-04-02.html"
)
rds_path <- file.path(tmp_dir, "gmds_misf_seeded_pipeline_results.rds")

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
  if (identical(row$method_id[[1L]], "reference")) {
    return("Reference surface")
  }
  sprintf(
    "%s\nsigma %s, rho %s\ntotal %ss",
    row$method_label[[1L]],
    fmt_num(row$gmds_stress[[1L]], 4L),
    fmt_num(row$procrustes_rmse[[1L]], 4L),
    fmt_time(row$elapsed_sec[[1L]])
  )
}

method_note <- function(case_result, row) {
  if (identical(row$method_id[[1L]], "reference")) {
    return(sprintf("%s. n = %d vertices.", case_result$case$label, row$n[[1L]]))
  }
  if (grepl("^misf_", row$method_id[[1L]])) {
    return(sprintf(
      "%s. seed source: %s. t_seed = %ss, t_pipe = %ss.",
      row$method_label[[1L]],
      ifelse(nzchar(row$seed_source[[1L]]), row$seed_source[[1L]], "internal top solve"),
      fmt_time(row$seed_elapsed_sec[[1L]]),
      fmt_time(row$pipeline_elapsed_sec[[1L]])
    ))
  }
  sprintf("%s. n = %d vertices.", row$method_label[[1L]], row$n[[1L]])
}

stage_card_title <- function(method_label, stage_name) {
  label <- switch(stage_name,
    top_level = "top level",
    after_insertion = "after insertion",
    after_refinement = "after refinement",
    final = "final",
    stage_name
  )
  sprintf("%s\n%s", method_label, label)
}

case_sections <- lapply(bundle$case_results, function(case_result) {
  final_cards <- lapply(case_result$methods, function(method_result) {
    row <- method_result$metrics[1L, , drop = FALSE]
    widget <- make_widget(
      coords = method_result$display_coords,
      edges = case_result$case$edges,
      vertex_col = if (identical(row$method_id[[1L]], "reference")) "#bc6c25" else "#3a5a40"
    )
    card_div(
      title_text = build_panel_title(method_result),
      widget = widget,
      note = method_note(case_result, row)
    )
  })

  stage_methods <- Filter(function(method) !is.null(method$stage_display), case_result$methods)
  stage_sections <- lapply(stage_methods, function(method_result) {
    row <- method_result$metrics[1L, , drop = FALSE]
    stage_names <- c("top_level", "after_insertion", "after_refinement", "final")
    cards <- lapply(stage_names, function(stage_name) {
      card_div(
        title_text = stage_card_title(row$method_label[[1L]], stage_name),
        widget = make_widget(
          coords = method_result$stage_display[[stage_name]],
          edges = case_result$case$edges,
          vertex_col = "#355070"
        ),
        note = if (identical(stage_name, "top_level")) {
          ifelse(nzchar(row$seed_source[[1L]]), paste("Top-level seed source:", row$seed_source[[1L]]), "Internal MISF top-level solve")
        } else {
          NULL
        }
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
      "All non-reference layouts are Procrustes-aligned to the reference paraboloid before rendering, so the panels differ by intrinsic shape rather than by rigid viewing transforms. The stage section below shows how each MISF pipeline evolves from its top-level scaffold to the final polished layout."
    ),
    htmltools::tags$h3("Final Layouts"),
    htmltools::tags$div(class = "panel-grid", final_cards),
    if (length(stage_sections)) {
      htmltools::tagList(
        htmltools::tags$h3("MISF Stage Layouts"),
        stage_sections
      )
    }
  )
})

page <- htmltools::browsable(
  htmltools::tagList(
    htmltools::tags$head(
      htmltools::tags$title("Interactive 3D Panels: Phase C Seeded MISF-GMDS"),
      htmltools::tags$style(htmltools::HTML(
        "
        body {
          margin: 0;
          background: #f5f1ea;
          color: #1f2933;
          font-family: 'Iowan Old Style', 'Palatino Linotype', 'Book Antiqua', Palatino, serif;
        }
        .page {
          max-width: 1620px;
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
          max-width: 1120px;
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
        .case-section h3 {
          margin: 16px 0 10px;
          font-size: 22px;
          color: #5a4636;
        }
        .case-lead {
          margin: 0 0 14px;
          color: #5f5445;
          font-size: 15px;
          line-height: 1.5;
        }
        .panel-grid {
          display: grid;
          grid-template-columns: repeat(auto-fit, minmax(420px, 1fr));
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
          color: #6b3e26;
        }
        .panel-card-line-metric {
          margin-top: 4px;
          font-size: 14px;
          color: #4b5563;
        }
        .panel-widget {
          padding: 14px 14px 10px;
          display: flex;
          justify-content: center;
          background: #fff;
        }
        .panel-card-note {
          margin: 0;
          padding: 0 18px 16px;
          color: #5f5445;
          font-size: 14px;
          line-height: 1.5;
        }
        .trace-block {
          margin-top: 12px;
          padding-top: 12px;
          border-top: 1px solid #d9d3c7;
        }
        "
      ))
    ),
    htmltools::tags$div(
      class = "page",
      htmltools::tags$section(
        class = "hero",
        htmltools::tags$h1("Phase C: Seeded MISF-GMDS Full-Pipeline Gallery"),
        htmltools::tags$p(
          "This interactive companion shows the full Phase C comparison. Each case includes the final layouts for the direct cmdscale control, the current MISF-GMDS baseline, and the two seeded MISF variants that reuse the best Phase B top-level corrected coarse embeddings."
        ),
        htmltools::tags$p(
          "The stage panels make the main question easy to inspect in 3D: once we replace the coarse top-level scaffold, do insertion, sparse refinement, and final polish stay in a better basin all the way to the full-graph embedding?"
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

message("Wrote Phase C interactive HTML: ", output_html)
