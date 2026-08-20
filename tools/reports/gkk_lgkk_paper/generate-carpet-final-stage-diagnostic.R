#!/usr/bin/env Rscript

output_root <- file.path("output", "gkk_lgkk_paper", "tmp", "carpet-final-stage-diagnostic")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE)
} else if (requireNamespace("grip", quietly = TRUE)) {
  library(grip)
} else {
  stop("Install 'devtools' or the 'grip' package to run this diagnostic.")
}

helper_env <- new.env(parent = baseenv())
sys.source(file.path("tools", "benchmarks", "gkk_lgkk_paper", "benchmark-sierpinski-baseline.R"), envir = helper_env)

level <- 4L
seed <- 6L
built <- helper_env$build_sierpinski_carpet(level)
edges <- built$edges
canonical <- built$coords
n <- nrow(canonical)

trace_fr <- grip.layout.trace(
  edges = edges,
  n = n,
  dim = 2,
  seed = seed,
  final_mode = "fr",
  trace = "round",
  trace.every = 1
)

pre_final_idx <- which(trace_fr$meta$phase == "level_start" &
                         trace_fr$meta$misf_level == 0L)
if (length(pre_final_idx) == 0L) {
  stop("Could not locate the pre-final full-graph snapshot in the trace.")
}
pre_final_coords <- trace_fr$frames[[pre_final_idx[[1L]]]]

coords_fr <- trace_fr$final
coords_kk_repulse <- grip.layout(
  edges = edges,
  n = n,
  dim = 2,
  seed = seed,
  final_mode = "kk_repulse"
)

diag_pre <- grip.geometry.diagnostics(
  coords = pre_final_coords,
  target.coords = canonical,
  edges = edges,
  family = "sierpinski.carpet",
  rng.seed = 1L
)
diag_fr <- grip.geometry.diagnostics(
  coords = coords_fr,
  target.coords = canonical,
  edges = edges,
  family = "sierpinski.carpet",
  rng.seed = 1L
)
diag_kk <- grip.geometry.diagnostics(
  coords = coords_kk_repulse,
  target.coords = canonical,
  edges = edges,
  family = "sierpinski.carpet",
  rng.seed = 1L
)

fmt <- function(x) {
  ifelse(is.finite(x), format(round(x, 4L), nsmall = 4L, trim = TRUE), "NA")
}

panel_subtitle <- function(diag) {
  paste0(
    "RMSE ", fmt(diag$procrustes.rmse[[1L]]),
    " | axis dev ", fmt(diag$edge.axis.deviation[[1L]]),
    " | boundary ", fmt(diag$boundary.waviness[[1L]])
  )
}

panels <- list(
  list(
    coords = canonical,
    title = "Canonical",
    subtitle = "target embedding"
  ),
  list(
    coords = pre_final_coords,
    title = "No Final FR",
    subtitle = paste0("pre-final snapshot | ", panel_subtitle(diag_pre))
  ),
  list(
    coords = coords_fr,
    title = "Final FR",
    subtitle = panel_subtitle(diag_fr)
  ),
  list(
    coords = coords_kk_repulse,
    title = "KK + Explicit Repulsion",
    subtitle = panel_subtitle(diag_kk)
  )
)

png_path <- file.path(
  output_root,
  "sierpinski-carpet-level-4-canonical-no-final-fr-final-fr-kk-repulse.png"
)
pdf_path <- file.path(
  output_root,
  "sierpinski-carpet-level-4-canonical-no-final-fr-final-fr-kk-repulse.pdf"
)
csv_path <- file.path(
  output_root,
  "sierpinski-carpet-level-4-canonical-no-final-fr-final-fr-kk-repulse-metrics.csv"
)
md_path <- file.path(
  output_root,
  "sierpinski-carpet-level-4-canonical-no-final-fr-final-fr-kk-repulse-summary.md"
)

draw_sheet <- function(device_fun) {
  device_fun()
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(2, 2), mar = c(0, 0, 3.3, 0), xaxs = "i", yaxs = "i")
  for (panel in panels) {
    helper_env$plot_layout_panel(
      panel$coords,
      edges,
      title_text = panel$title,
      subtitle_text = panel$subtitle
    )
  }
}

draw_sheet(function() {
  grDevices::png(png_path, width = 2600, height = 2200, res = 180, bg = "#f7f3ea")
})
draw_sheet(function() {
  grDevices::pdf(pdf_path, width = 14.5, height = 12.0,
                 paper = "special", bg = "#f7f3ea", useDingbats = FALSE)
})

metrics <- rbind(
  data.frame(mode = "no_final_fr", diag_pre, stringsAsFactors = FALSE),
  data.frame(mode = "final_fr", diag_fr, stringsAsFactors = FALSE),
  data.frame(mode = "kk_repulse", diag_kk, stringsAsFactors = FALSE)
)
utils::write.csv(metrics, csv_path, row.names = FALSE)

writeLines(
  c(
    "# Carpet Final-Stage Diagnostic",
    "",
    sprintf("- graph: `Sierpinski carpet level %d`", level),
    sprintf("- seed: `%d`", seed),
    "- panel 2 (`No Final FR`) uses the trace snapshot taken at the start of the full-graph level, before any final-stage FR rounds are applied.",
    "",
    "## Outputs",
    "",
    sprintf("- PNG: `%s`", png_path),
    sprintf("- PDF: `%s`", pdf_path),
    sprintf("- Metrics CSV: `%s`", csv_path),
    "",
    "## Geometry Diagnostics",
    "",
    paste(capture.output(print(metrics)), collapse = "\n")
  ),
  con = md_path
)

message("Wrote:")
message("  ", png_path)
message("  ", pdf_path)
message("  ", csv_path)
message("  ", md_path)
