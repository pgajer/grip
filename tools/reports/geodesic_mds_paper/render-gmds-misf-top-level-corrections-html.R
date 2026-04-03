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
tmp_dir <- file.path(design_root, "tmp", "gmds-misf-top-level-corrections-2026-04-02")
interactive_dir <- file.path(design_root, "interactive-prototypes")
output_html <- file.path(
  interactive_dir,
  "gmds_misf_top_level_corrections_3d_2026-04-02.html"
)
rds_path <- file.path(tmp_dir, "gmds_misf_top_level_corrections_results.rds")

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
                        triangles = NULL,
                        vertex_col = "#355070",
                        edge_col = "#adb5bd",
                        face_col = "#84a98c",
                        width = 420,
                        height = 320,
                        theta = 35,
                        phi = 24,
                        zoom = 0.84,
                        fov = 30) {
  coords <- as.matrix(coords)
  edges <- as.matrix(edges)
  triangles <- if (is.null(triangles)) matrix(integer(), ncol = 3L) else as.matrix(triangles)

  keep <- apply(coords, 1L, function(x) all(is.finite(x)))
  good_edges <- if (nrow(edges)) keep[edges[, 1L]] & keep[edges[, 2L]] else logical(0L)
  edges <- if (nrow(edges)) edges[good_edges, , drop = FALSE] else edges
  good_tri <- if (nrow(triangles)) {
    keep[triangles[, 1L]] & keep[triangles[, 2L]] & keep[triangles[, 3L]]
  } else {
    logical(0L)
  }
  triangles <- if (nrow(triangles)) triangles[good_tri, , drop = FALSE] else triangles

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

  if (nrow(triangles) > 0L) {
    tri_coords <- coords[t(triangles), , drop = FALSE]
    rgl::triangles3d(
      tri_coords[, 1L],
      tri_coords[, 2L],
      tri_coords[, 3L],
      color = grDevices::adjustcolor(face_col, alpha.f = 0.18),
      front = "fill",
      back = "lines"
    )
  }
  if (nrow(edges) > 0L) {
    rgl::lines3d(seg, color = grDevices::adjustcolor(edge_col, alpha.f = 0.82), lwd = 1.5)
  }
  if (any(keep)) {
    rgl::points3d(
      coords[keep, 1L],
      coords[keep, 2L],
      coords[keep, 3L],
      color = vertex_col,
      size = 7,
      alpha = 0.98
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

build_stage_title <- function(entry, stage_name = NULL) {
  metrics <- entry$metrics
  main <- if (is.null(stage_name)) {
    metrics$initializer_label[[1L]]
  } else {
    sprintf("%s (%s)", metrics$initializer_label[[1L]], stage_name)
  }
  sprintf(
    "%s\nsigma %s, rho %s\nt %ss",
    main,
    fmt_num(metrics$sigma_geo[[1L]], 4L),
    fmt_num(metrics$rho[[1L]], 4L),
    fmt_time(metrics$elapsed_sec[[1L]])
  )
}

build_trace_title <- function(method_label, stage_name, trace_frame) {
  sprintf(
    "%s (%s)\niter %s, sigma %s\nanchor %s",
    method_label,
    stage_name,
    ifelse(is.na(trace_frame$iteration), "NA", as.character(trace_frame$iteration)),
    fmt_num(trace_frame$gmds_energy, 4L),
    fmt_num(trace_frame$anchor_weight, 3L)
  )
}

reference_card <- function(case) {
  card_div(
    title_text = "Reference sample\nactive-set target geometry",
    widget = make_widget(
      coords = case$truth,
      edges = case$display_edges,
      triangles = case$display_triangles,
      vertex_col = "#bc6c25",
      face_col = "#dda15e"
    ),
    note = sprintf("Reference active-set sample. L = %d, |V_L| = %d.", case$top_level, case$top_n)
  )
}

make_stage_cards <- function(case_result, stage, stage_name, vertex_col, face_col) {
  c(
    list(reference_card(case_result$case)),
    lapply(case_result$methods, function(method) {
      entry <- method[[stage]]
      card_div(
        title_text = build_stage_title(entry, stage_name = stage_name),
        widget = make_widget(
          coords = entry$display_coords,
          edges = case_result$case$display_edges,
          triangles = case_result$case$display_triangles,
          vertex_col = vertex_col,
          face_col = face_col
        ),
        note = entry$metrics$note[[1L]]
      )
    })
  )
}

case_sections <- lapply(bundle$case_results, function(case_result) {
  trace_sections <- lapply(case_result$methods, function(method) {
    stage_blocks <- list()

    if (length(method$anchor$trace_selected)) {
      anchor_cards <- lapply(method$anchor$trace_selected, function(tr) {
        card_div(
          title_text = build_trace_title(method$method_label, "Anchor-GMDS", tr),
          widget = make_widget(
            coords = tr$display_coords,
            edges = case_result$case$display_edges,
            triangles = case_result$case$display_triangles,
            vertex_col = "#3d5a80",
            face_col = "#98c1d9"
          )
        )
      })
      stage_blocks[[length(stage_blocks) + 1L]] <- htmltools::tagList(
        htmltools::tags$h4("Anchor-GMDS snapshots"),
        htmltools::tags$div(class = "panel-grid", anchor_cards)
      )
    }

    if (length(method$pure$trace_selected)) {
      pure_cards <- lapply(method$pure$trace_selected, function(tr) {
        card_div(
          title_text = build_trace_title(method$method_label, "Pure follow-up", tr),
          widget = make_widget(
            coords = tr$display_coords,
            edges = case_result$case$display_edges,
            triangles = case_result$case$display_triangles,
            vertex_col = "#4d908e",
            face_col = "#90be6d"
          )
        )
      })
      stage_blocks[[length(stage_blocks) + 1L]] <- htmltools::tagList(
        htmltools::tags$h4("Pure-GMDS follow-up snapshots"),
        htmltools::tags$div(class = "panel-grid", pure_cards)
      )
    }

    if (!length(stage_blocks)) {
      return(NULL)
    }

    htmltools::tags$div(
      class = "trace-block",
      htmltools::tags$h3(method$method_label),
      stage_blocks
    )
  })
  trace_sections <- Filter(Negate(is.null), trace_sections)

  htmltools::tags$section(
    class = "case-section",
    htmltools::tags$h2(case_result$case$label),
    htmltools::tags$p(
      class = "case-lead",
      "Every non-reference panel is Procrustes-aligned to the reference active-set sample before rendering. The widgets show the same active-set triangulation across all stages so the coarse surface basin, the anchor-preserving correction, and the pure-GMDS follow-up can be compared directly."
    ),
    htmltools::tags$h3("Initializer Basins"),
    htmltools::tags$div(
      class = "panel-grid",
      make_stage_cards(case_result, stage = "initializer", stage_name = "Initializer", vertex_col = "#355070", face_col = "#84a98c")
    ),
    htmltools::tags$h3("After Anchor-Regularized GMDS"),
    htmltools::tags$div(
      class = "panel-grid",
      make_stage_cards(case_result, stage = "anchor", stage_name = "Anchor-GMDS", vertex_col = "#3d5a80", face_col = "#98c1d9")
    ),
    htmltools::tags$h3("After Pure-GMDS Follow-Up"),
    htmltools::tags$div(
      class = "panel-grid",
      make_stage_cards(case_result, stage = "pure", stage_name = "Pure follow-up", vertex_col = "#4d908e", face_col = "#90be6d")
    ),
    if (length(trace_sections)) {
      htmltools::tagList(
        htmltools::tags$h3("Selected GMDS Correction Snapshots"),
        trace_sections
      )
    }
  )
})

page <- htmltools::browsable(
  htmltools::tagList(
    htmltools::tags$head(
      htmltools::tags$title("Phase B: MISF-GMDS Top-Level Corrections 3D Gallery"),
      htmltools::tags$style(htmltools::HTML(
        "
        body {
          margin: 0;
          background: #f4efe7;
          color: #1f2933;
          font-family: 'Iowan Old Style', 'Palatino Linotype', 'Book Antiqua', Palatino, serif;
        }
        .page {
          max-width: 1680px;
          margin: 0 auto;
          padding: 28px 28px 40px;
        }
        .hero {
          background: linear-gradient(135deg, #fffdf9 0%, #f2e8db 100%);
          border: 1px solid #d7c8b6;
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
          max-width: 1220px;
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
          margin: 18px 0 10px;
          font-size: 22px;
          color: #5a4636;
        }
        .case-section h4 {
          margin: 12px 0 10px;
          font-size: 18px;
          color: #6b5b4d;
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
          margin-top: 10px;
          padding-top: 10px;
          border-top: 1px solid #d9d3c7;
        }
        "
      ))
    ),
    htmltools::tags$div(
      class = "page",
      htmltools::tags$section(
        class = "hero",
        htmltools::tags$h1("Phase B: MISF-GMDS Top-Level Correction Gallery"),
        htmltools::tags$p(
          sprintf(
            "This interactive companion shows the full Group B experiment on the coarse MISF active sets. Every Phase A initializer is reused as a basin, corrected by a short anchor-only GMDS stage with lambda_A %.2f -> %.2f, and then polished by a short pure-GMDS follow-up.",
            bundle$cfg$anchor_weight,
            bundle$cfg$anchor_weight_end
          )
        ),
        htmltools::tags$p(
          "Each case is organized into three final-stage galleries plus selected correction snapshots. The underlying scoring graph is the complete weighted top-level graph, while the displayed surface is the fixed active-set triangulation used throughout the report so the coarse 3D geometry is easier to inspect."
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

message("Wrote Phase B interactive HTML: ", output_html)
