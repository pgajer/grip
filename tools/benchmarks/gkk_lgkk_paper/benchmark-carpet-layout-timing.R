#!/usr/bin/env Rscript

manual_root <- file.path("output", "gkk_lgkk_paper", "tmp", "carpet-layout-timing-level4-5")
dir.create(manual_root, recursive = TRUE, showWarnings = FALSE)

if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE)
} else if (requireNamespace("grip", quietly = TRUE)) {
  library(grip)
} else {
  stop("Install 'devtools' or the 'grip' package to run the carpet timing benchmark.")
}

if (!requireNamespace("igraph", quietly = TRUE)) {
  stop("Package 'igraph' is required to run the carpet timing benchmark.")
}
if (!requireNamespace("callr", quietly = TRUE)) {
  stop("Package 'callr' is required to run the carpet timing benchmark.")
}

repo_root <- normalizePath(".", winslash = "/", mustWork = TRUE)

graph_configs <- data.frame(
  level = c(4L, 5L),
  reps = c(3L, 1L),
  timeout_sec = c(180L, 300L),
  stringsAsFactors = FALSE
)

method_configs <- data.frame(
  method_id = c("grip.layout", "igraph_fr", "igraph_kk"),
  method_label = c("grip.layout()", "igraph::layout_with_fr()", "igraph::layout_with_kk()"),
  stringsAsFactors = FALSE
)

format_num <- function(x, digits = 3L) {
  ifelse(is.finite(x), format(round(x, digits), nsmall = digits, trim = TRUE), "NA")
}

run_timed_layout <- function(method_id, edges, n, seed, timeout_sec, repo_root) {
  tryCatch(
    {
      result <- callr::r(
        function(method_id, edges, n, seed, repo_root) {
          suppressPackageStartupMessages({
            if (identical(method_id, "grip.layout")) {
              if (requireNamespace("devtools", quietly = TRUE)) {
                devtools::load_all(repo_root, quiet = TRUE, export_all = FALSE, helpers = FALSE)
              } else {
                library(grip)
              }
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
            elapsed_sec = elapsed,
            n_vertices = nrow(coords),
            dim = ncol(coords)
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

      list(
        status = "ok",
        elapsed_sec = result$elapsed_sec,
        n_vertices = result$n_vertices,
        dim = result$dim,
        message = NA_character_
      )
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
        n_vertices = n,
        dim = 2L,
        message = msg
      )
    }
  )
}

write_summary_markdown <- function(raw_metrics, summary_metrics, graph_info, path) {
  lines <- c(
    "# Carpet Layout Timing",
    "",
    "Methods:",
    "- `grip.layout()` with the current primary defaults",
    "- `igraph::layout_with_fr()` with igraph defaults (`niter = 500`)",
    "- `igraph::layout_with_kk()` with igraph defaults (`maxiter = 50 * vcount(graph)`)",
    "",
    "Graph set:",
    sprintf("- level 4 carpet: `%d` vertices, `%d` edges, `%d` timing runs per method, timeout `%ds` per run",
            graph_info$vertices[[1L]], graph_info$edges[[1L]], graph_info$reps[[1L]], graph_info$timeout_sec[[1L]]),
    sprintf("- level 5 carpet: `%d` vertices, `%d` edges, `%d` timing run per method, timeout `%ds` per run",
            graph_info$vertices[[2L]], graph_info$edges[[2L]], graph_info$reps[[2L]], graph_info$timeout_sec[[2L]]),
    "",
    "Timing summary:",
    "",
    "| Level | Method | Requested runs | Completed | Mean sec | Median sec | Min sec | Max sec | SD sec | Timeout count | Error count |",
    "| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"
  )

  for (i in seq_len(nrow(summary_metrics))) {
    row <- summary_metrics[i, , drop = FALSE]
    lines <- c(lines, sprintf(
      "| %d | %s | %d | %d | %s | %s | %s | %s | %s | %d | %d |",
      row$level,
      row$method_label,
      row$runs_requested,
      row$runs_completed,
      format_num(row$elapsed_mean, 3L),
      format_num(row$elapsed_median, 3L),
      format_num(row$elapsed_min, 3L),
      format_num(row$elapsed_max, 3L),
      format_num(row$elapsed_sd, 3L),
      row$timeout_count,
      row$error_count
    ))
  }

  lines <- c(
    lines,
    "",
    "Per-run results:",
    "",
    "| Level | Method | Seed | Status | Elapsed sec | Message |",
    "| ---: | --- | ---: | --- | ---: | --- |"
  )

  for (i in seq_len(nrow(raw_metrics))) {
    row <- raw_metrics[i, , drop = FALSE]
    lines <- c(lines, sprintf(
      "| %d | %s | %d | %s | %s | %s |",
      row$level,
      row$method_label,
      row$seed,
      row$status,
      format_num(row$elapsed_sec, 3L),
      if (is.na(row$message) || !nzchar(row$message)) "" else gsub("\\|", "/", row$message)
    ))
  }

  writeLines(lines, con = path)
}

graph_info <- do.call(
  rbind,
  lapply(seq_len(nrow(graph_configs)), function(i) {
    cfg <- graph_configs[i, , drop = FALSE]
    edges <- edges.sierpinski.carpet(cfg$level[[1L]])
    data.frame(
      level = cfg$level[[1L]],
      vertices = max(edges),
      edges = nrow(edges),
      reps = cfg$reps[[1L]],
      timeout_sec = cfg$timeout_sec[[1L]],
      stringsAsFactors = FALSE
    )
  })
)

raw_metrics <- list()
for (i in seq_len(nrow(graph_configs))) {
  cfg <- graph_configs[i, , drop = FALSE]
  level <- cfg$level[[1L]]
  reps <- cfg$reps[[1L]]
  timeout_sec <- cfg$timeout_sec[[1L]]
  edges <- edges.sierpinski.carpet(level)
  n <- max(edges)
  m <- nrow(edges)

  for (j in seq_len(nrow(method_configs))) {
    method <- method_configs[j, , drop = FALSE]
    for (seed in seq_len(reps)) {
      message(sprintf(
        "Timing %s on carpet level %d (n=%d, m=%d), seed %d/%d...",
        method$method_label[[1L]], level, n, m, seed, reps
      ))
      timed <- run_timed_layout(
        method_id = method$method_id[[1L]],
        edges = edges,
        n = n,
        seed = seed,
        timeout_sec = timeout_sec,
        repo_root = repo_root
      )
      raw_metrics[[length(raw_metrics) + 1L]] <- data.frame(
        level = level,
        vertices = n,
        edges = m,
        method_id = method$method_id[[1L]],
        method_label = method$method_label[[1L]],
        seed = seed,
        timeout_sec = timeout_sec,
        status = timed$status,
        elapsed_sec = timed$elapsed_sec,
        message = timed$message,
        stringsAsFactors = FALSE
      )
    }
  }
}

raw_metrics <- do.call(rbind, raw_metrics)
raw_metrics <- raw_metrics[order(raw_metrics$level, raw_metrics$method_label, raw_metrics$seed), , drop = FALSE]

summary_metrics <- do.call(
  rbind,
  lapply(split(raw_metrics, paste(raw_metrics$level, raw_metrics$method_label, sep = "::")), function(df) {
    ok <- df[df$status == "ok" & is.finite(df$elapsed_sec), , drop = FALSE]
    data.frame(
      level = df$level[[1L]],
      vertices = df$vertices[[1L]],
      edges = df$edges[[1L]],
      method_id = df$method_id[[1L]],
      method_label = df$method_label[[1L]],
      runs_requested = nrow(df),
      runs_completed = nrow(ok),
      timeout_count = sum(df$status == "timeout"),
      error_count = sum(df$status == "error"),
      elapsed_mean = if (nrow(ok) > 0L) mean(ok$elapsed_sec) else NA_real_,
      elapsed_median = if (nrow(ok) > 0L) stats::median(ok$elapsed_sec) else NA_real_,
      elapsed_min = if (nrow(ok) > 0L) min(ok$elapsed_sec) else NA_real_,
      elapsed_max = if (nrow(ok) > 0L) max(ok$elapsed_sec) else NA_real_,
      elapsed_sd = if (nrow(ok) > 1L) stats::sd(ok$elapsed_sec) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
)
summary_metrics <- summary_metrics[order(summary_metrics$level, summary_metrics$method_label), , drop = FALSE]

raw_csv_path <- file.path(manual_root, "carpet-layout-timing-raw.csv")
summary_csv_path <- file.path(manual_root, "carpet-layout-timing-summary.csv")
summary_md_path <- file.path(manual_root, "carpet-layout-timing-summary.md")

utils::write.csv(raw_metrics, raw_csv_path, row.names = FALSE)
utils::write.csv(summary_metrics, summary_csv_path, row.names = FALSE)
write_summary_markdown(raw_metrics, summary_metrics, graph_info, summary_md_path)

message(sprintf("Raw timing results written to %s", raw_csv_path))
message(sprintf("Timing summary written to %s", summary_csv_path))
message(sprintf("Markdown report written to %s", summary_md_path))
