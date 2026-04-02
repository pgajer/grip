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
tmp_dir <- file.path(manual_root, "tmp", "gmds-paraboloid-identifiability-2026-04-01")
interactive_dir <- file.path(manual_root, "interactive-prototypes")
output_html <- file.path(
  interactive_dir,
  "gmds_paraboloid_identifiability_all_embeddings_2026-04-01.html"
)
rds_path <- file.path(tmp_dir, "paraboloid_identifiability_embeddings.rds")

dir.create(interactive_dir, recursive = TRUE, showWarnings = FALSE)

for (pkg in c("htmltools", "htmlwidgets", "rgl")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required to render the interactive HTML.", pkg))
  }
}

if (!file.exists(rds_path)) {
  stop("Expected embeddings bundle at: ", rds_path)
}

bundle <- readRDS(rds_path)

fmt_num <- function(x, digits = 4L) {
  ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "NA")
}

make_widget <- function(coords,
                        edges,
                        vertex_col = "#3a5a40",
                        edge_col = "#adb5bd",
                        width = 320,
                        height = 245,
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
    rgl::lines3d(seg, color = grDevices::adjustcolor(edge_col, alpha.f = 0.78), lwd = 1.4)
  }
  if (any(keep)) {
    rgl::points3d(
      coords[keep, 1L],
      coords[keep, 2L],
      coords[keep, 3L],
      color = vertex_col,
      size = 5.5,
      alpha = 0.96
    )
  }

  rgl::aspect3d(1, 1, 1)
  rgl::view3d(theta = theta, phi = phi, zoom = zoom, fov = fov)
  rgl::rglwidget(width = width, height = height)
}

card_div <- function(title_lines, widget, note) {
  htmltools::tags$div(
    class = "panel-card",
    htmltools::tags$div(
      class = "panel-card-title",
      htmltools::tags$div(class = "panel-card-line panel-card-line-main", title_lines[[1L]]),
      if (length(title_lines) > 1L) {
        lapply(title_lines[-1L], function(x) {
          htmltools::tags$div(class = "panel-card-line panel-card-line-metric", x)
        })
      }
    ),
    htmltools::tags$div(class = "panel-widget", widget),
    htmltools::tags$p(class = "panel-card-note", note)
  )
}

embedding_title_lines <- function(embedding) {
  lines <- c(embedding$label)
  if (!is.null(embedding$angle_deg)) {
    lines <- c(lines, sprintf("fold angle %d deg", embedding$angle_deg))
  }
  if (!is.null(embedding$gmds_stress)) {
    lines <- c(lines, sprintf("sigma %s, rho %s", fmt_num(embedding$gmds_stress, 4L), fmt_num(embedding$rmse_to_truth, 4L)))
  }
  if (!is.null(embedding$elapsed_sec)) {
    lines <- c(lines, sprintf("elapsed %s s", fmt_num(embedding$elapsed_sec, 3L)))
  }
  lines
}

embedding_note <- function(embedding) {
  switch(
    embedding$kind,
    reference = "Exact scaled reference paraboloid under the graph generator's edge-weight normalization.",
    init = "Classical-MDS initialization aligned to the reference surface for visual comparison.",
    fold = "Controlled hinge-fold perturbation of the exact reference surface; no optimization applied.",
    recovery = "Final pure-GMDS embedding after 25 iterations from the stated start.",
    "Embedding from the identifiability experiment."
  )
}

case_sections <- lapply(bundle$cases, function(case) {
  cards <- lapply(case$embeddings, function(embedding) {
    widget <- make_widget(embedding$display_coords, case$edges)
    card_div(embedding_title_lines(embedding), widget, embedding_note(embedding))
  })
  htmltools::tags$section(
    class = "case-section",
    htmltools::tags$h2(case$label),
    htmltools::tags$p(
      class = "case-lead",
      "All embeddings from the identifiability diagnostic are shown below after alignment to the exact reference surface, so the local fold perturbations, the cMDS initialization, and the pure-GMDS recovery basins can be inspected directly."
    ),
    htmltools::tags$div(class = "panel-grid", cards)
  )
})

page <- htmltools::browsable(
  htmltools::tagList(
    htmltools::tags$head(
      htmltools::tags$title("Interactive 3D Panels: Paraboloid Identifiability Diagnostic"),
      htmltools::tags$style(htmltools::HTML(
        "
        body {
          margin: 0;
          background: #f5f1ea;
          color: #1f2933;
          font-family: 'Iowan Old Style', 'Palatino Linotype', 'Book Antiqua', Palatino, serif;
        }
        .page {
          max-width: 1680px;
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
          max-width: 1200px;
        }
        .case-section {
          margin-top: 28px;
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
          grid-template-columns: repeat(auto-fit, minmax(320px, 1fr));
          gap: 16px;
        }
        .panel-card {
          background: #fffdfb;
          border: 1px solid #d9d3c7;
          border-radius: 16px;
          box-shadow: 0 8px 22px rgba(51, 41, 28, 0.08);
          overflow: hidden;
        }
        .panel-card-title {
          padding: 14px 16px 8px;
          border-bottom: 1px solid #ede6da;
          background: linear-gradient(180deg, #fffdfa 0%, #faf5ed 100%);
          min-height: 104px;
        }
        .panel-card-line-main {
          font-size: 18px;
          font-weight: 700;
          color: #6b3e26;
        }
        .panel-card-line-metric {
          margin-top: 4px;
          font-size: 13px;
          color: #4b5563;
        }
        .panel-widget {
          padding: 12px 12px 8px;
          display: flex;
          justify-content: center;
          background: #fff;
        }
        .panel-card-note {
          margin: 0;
          padding: 0 16px 16px;
          color: #5f5445;
          font-size: 13px;
          line-height: 1.45;
        }
        "
      ))
    ),
    htmltools::tags$div(
      class = "page",
      htmltools::tags$section(
        class = "hero",
        htmltools::tags$h1("Interactive 3D Panels: Paraboloid Identifiability Diagnostic"),
        htmltools::tags$p(
          "This page collects every generated embedding from the identifiability experiment: the exact scaled reference surfaces, the classical-MDS starts, the full fold-angle sweep, and the pure-GMDS recovery endpoints from both cMDS and a folded initialization."
        )
      ),
      case_sections
    )
  )
)

htmltools::save_html(page, file = output_html, background = "white")
message("Wrote interactive HTML: ", output_html)
