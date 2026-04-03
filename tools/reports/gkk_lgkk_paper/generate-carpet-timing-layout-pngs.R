#!/usr/bin/env Rscript

manual_root <- file.path("dev", "manual", "tmp", "carpet-layout-timing-level4-5")
png_dir <- file.path(manual_root, "png")
dir.create(manual_root, recursive = TRUE, showWarnings = FALSE)
dir.create(png_dir, recursive = TRUE, showWarnings = FALSE)

repo_root <- normalizePath(".", winslash = "/", mustWork = TRUE)

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Package 'devtools' is required to generate the carpet timing PNGs.")
}
if (!requireNamespace("igraph", quietly = TRUE)) {
  stop("Package 'igraph' is required to generate the carpet timing PNGs.")
}
if (!requireNamespace("callr", quietly = TRUE)) {
  stop("Package 'callr' is required to generate the carpet timing PNGs.")
}

devtools::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE)
sys.source(file.path("tools", "benchmarks", "gkk_lgkk_paper", "benchmark-sierpinski-baseline.R"), envir = environment())

layout_cases <- data.frame(
  level = c(4L, 4L, 4L, 5L, 5L, 5L),
  method_id = c("grip.layout", "igraph_fr", "igraph_kk",
                "grip.layout", "igraph_fr", "igraph_kk"),
  method_label = c("grip.layout()", "igraph::layout_with_fr()", "igraph::layout_with_kk()",
                   "grip.layout()", "igraph::layout_with_fr()", "igraph::layout_with_kk()"),
  timeout_sec = c(180L, 180L, 180L, 300L, 300L, 300L),
  seed = c(1L, 1L, 1L, 1L, 1L, 1L),
  stringsAsFactors = FALSE
)

run_layout_case <- function(method_id, edges, n, seed, timeout_sec, repo_root) {
  tryCatch(
    {
      result <- callr::r(
        function(method_id, edges, n, seed, repo_root) {
          suppressPackageStartupMessages({
            if (identical(method_id, "grip.layout")) {
              devtools::load_all(repo_root, quiet = TRUE, export_all = FALSE, helpers = FALSE)
            } else {
              library(igraph)
            }
          })

          gc()
          if (identical(method_id, "grip.layout")) {
            grip_fun <- get("grip.layout", envir = asNamespace("grip"))
            started <- proc.time()[["elapsed"]]
            coords <- grip_fun(edges = edges, n = n, dim = 2, seed = seed)
            elapsed <- proc.time()[["elapsed"]] - started
          } else {
            graph <- igraph::graph_from_edgelist(as.matrix(edges), directed = FALSE)
            set.seed(seed)
            started <- proc.time()[["elapsed"]]
            coords <- switch(
              method_id,
              igraph_fr = igraph::layout_with_fr(graph),
              igraph_kk = igraph::layout_with_kk(graph),
              stop(sprintf("Unknown method_id: %s", method_id))
            )
            elapsed <- proc.time()[["elapsed"]] - started
          }

          list(
            status = "ok",
            elapsed_sec = elapsed,
            coords = coords
          )
        },
        args = list(
          method_id = method_id,
          edges = edges,
          n = n,
          seed = seed,
          repo_root = repo_root
        ),
        timeout = timeout_sec,
        spinner = FALSE
      )
      result
    },
    error = function(e) {
      msg <- conditionMessage(e)
      status <- if (inherits(e, "callr_timeout_error") ||
                    grepl("timed out|timeout", msg, ignore.case = TRUE)) {
        "timeout"
      } else {
        "error"
      }
      list(
        status = status,
        elapsed_sec = NA_real_,
        coords = NULL,
        message = msg
      )
    }
  )
}

format_elapsed <- function(x) {
  if (!is.finite(x)) return("elapsed=NA")
  sprintf("elapsed=%.3fs", x)
}

plot_timeout_panel <- function(title_text, subtitle_text, color = "#8f1d21") {
  plot.new()
  graphics::par(usr = c(0, 1, 0, 1))
  graphics::rect(0, 0, 1, 1, col = "#f7f3ea", border = NA)
  graphics::text(0.5, 0.62, title_text, cex = 1.25, font = 2, col = "#16324f")
  graphics::text(0.5, 0.46, subtitle_text, cex = 0.95, col = "#466074")
  graphics::text(0.5, 0.24, "No layout produced within timeout", cex = 1.0, col = color)
}

draw_layout_png <- function(path, spec, result, aligned_coords = NULL) {
  grDevices::png(path, width = 1200, height = 1200, res = 180, bg = "#f7f3ea")
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mar = c(0, 0, 2.8, 0), xaxs = "i", yaxs = "i")

  title_text <- sprintf("Sierpinski carpet level %d", spec$level)
  subtitle_text <- sprintf("%s | seed=%d | %s",
                           spec$method_label,
                           spec$seed,
                           format_elapsed(result$elapsed_sec))

  if (identical(result$status, "ok") && !is.null(aligned_coords)) {
    plot_layout_panel(aligned_coords, spec$edges, title_text, subtitle_text)
  } else {
    plot_timeout_panel(title_text, sprintf("%s | seed=%d", spec$method_label, spec$seed))
  }
}

draw_contact_sheet <- function(path, case_results) {
  grDevices::png(path, width = 2400, height = 1600, res = 180, bg = "#f7f3ea")
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(2, 3), mar = c(0, 0, 2.8, 0), xaxs = "i", yaxs = "i")

  for (item in case_results) {
    spec <- item$spec
    result <- item$result
    title_text <- sprintf("Sierpinski carpet level %d", spec$level)
    subtitle_text <- sprintf("%s | seed=%d | %s",
                             spec$method_label,
                             spec$seed,
                             format_elapsed(result$elapsed_sec))
    if (identical(result$status, "ok") && !is.null(item$aligned_coords)) {
      plot_layout_panel(item$aligned_coords, spec$edges, title_text, subtitle_text)
    } else {
      plot_timeout_panel(title_text, sprintf("%s | seed=%d", spec$method_label, spec$seed))
    }
  }
}

case_results <- list()
summary_rows <- list()

for (i in seq_len(nrow(layout_cases))) {
  spec_row <- layout_cases[i, , drop = FALSE]
  level <- spec_row$level[[1L]]
  seed <- spec_row$seed[[1L]]
  timeout_sec <- spec_row$timeout_sec[[1L]]
  built <- build_sierpinski_carpet(level)
  edges <- built$edges
  n <- nrow(built$coords)

  message(sprintf(
    "Rendering %s on carpet level %d with seed %d (timeout %ds)...",
    spec_row$method_label[[1L]],
    level,
    seed,
    timeout_sec
  ))

  result <- run_layout_case(
    method_id = spec_row$method_id[[1L]],
    edges = edges,
    n = n,
    seed = seed,
    timeout_sec = timeout_sec,
    repo_root = repo_root
  )

  aligned_coords <- NULL
  if (identical(result$status, "ok") && !is.null(result$coords)) {
    aligned_coords <- align_to_target(result$coords, built$coords)$aligned
  }

  png_path <- file.path(
    png_dir,
    sprintf("sierpinski-carpet-level-%d-%s-seed-%d.png",
            level,
            gsub("[^a-z0-9]+", "-", tolower(spec_row$method_id[[1L]])),
            seed)
  )
  draw_layout_png(
    path = png_path,
    spec = list(
      level = level,
      method_label = spec_row$method_label[[1L]],
      seed = seed,
      edges = edges
    ),
    result = result,
    aligned_coords = aligned_coords
  )

  case_results[[length(case_results) + 1L]] <- list(
    spec = list(
      level = level,
      method_label = spec_row$method_label[[1L]],
      seed = seed,
      edges = edges
    ),
    result = result,
    aligned_coords = aligned_coords,
    png_path = png_path
  )

  summary_rows[[length(summary_rows) + 1L]] <- data.frame(
    level = level,
    method_id = spec_row$method_id[[1L]],
    method_label = spec_row$method_label[[1L]],
    seed = seed,
    timeout_sec = timeout_sec,
    status = result$status,
    elapsed_sec = result$elapsed_sec,
    png_path = png_path,
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summary_rows)
sheet_path <- file.path(manual_root, "sierpinski-carpet-level4-5-layouts-seed1-grid.png")
draw_contact_sheet(sheet_path, case_results)

summary_csv_path <- file.path(manual_root, "sierpinski-carpet-level4-5-layout-png-summary.csv")
summary_md_path <- file.path(manual_root, "sierpinski-carpet-level4-5-layout-png-summary.md")
utils::write.csv(summary_df, summary_csv_path, row.names = FALSE)

lines <- c(
  "# Carpet Layout PNGs",
  "",
  "This report renders the seed-1 layouts corresponding to the timing comparison for",
  "`grip.layout()`, `igraph::layout_with_fr()`, and `igraph::layout_with_kk()` on",
  "Sierpinski carpet levels 4 and 5.",
  "",
  sprintf("- contact sheet: `%s`", sheet_path),
  "",
  "| Level | Method | Status | Elapsed sec | PNG |",
  "| ---: | --- | --- | ---: | --- |"
)
for (i in seq_len(nrow(summary_df))) {
  row <- summary_df[i, , drop = FALSE]
  lines <- c(lines, sprintf(
    "| %d | %s | %s | %s | `%s` |",
    row$level,
    row$method_label,
    row$status,
    if (is.finite(row$elapsed_sec)) format(round(row$elapsed_sec, 3), nsmall = 3, trim = TRUE) else "NA",
    row$png_path
  ))
}
writeLines(lines, con = summary_md_path)

message(sprintf("Layout PNG summary written to %s", summary_csv_path))
message(sprintf("Markdown report written to %s", summary_md_path))
message(sprintf("Contact sheet written to %s", sheet_path))
