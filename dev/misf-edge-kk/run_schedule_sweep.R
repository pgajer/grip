#!/usr/bin/env Rscript

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Package 'pkgload' is required for schedule sweep", call. = FALSE)
}
pkgload::load_all(".", quiet = TRUE, export_all = TRUE, helpers = FALSE)

source("dev/misf-edge-kk/misf_edge_kk_objective_reference.R")
source("dev/misf-edge-kk/misf_edge_kk_constraint_reference.R")
source("dev/misf-edge-kk/misf_edge_kk_armijo_reference.R")
source("dev/misf-edge-kk/misf_edge_kk_level_reference.R")
source("dev/misf-edge-kk/misf_edge_kk_multilevel_reference.R")
source("dev/misf-edge-kk/misf_edge_kk_polish_reference.R")
source("dev/misf-edge-kk/misf_edge_kk_layout_reference.R")
source("dev/misf-edge-kk/misf_edge_kk_baseline_reference.R")

out_dir <- "output/misf-edge-kk-validation"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

parameter_grid <- expand.grid(
  rho_start = c(0, 0.2, 0.6, 1.0),
  lambda_start = c(0, 0.03, 0.1, 0.2),
  anchor_start = c(0, 0.1, 0.25),
  level_max_iter = c(8L, 20L),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

empty_sweep_row <- function(fixture, params, status, message,
                            elapsed = NA_real_) {
  data.frame(
    fixture = fixture$name,
    n = fixture$n,
    edge_count = nrow(fixture$edges),
    dim = fixture$dim,
    rho_start = params$rho_start,
    lambda_start = params$lambda_start,
    anchor_start = params$anchor_start,
    level_max_iter = params$level_max_iter,
    status = status,
    message = message,
    finite = FALSE,
    edge_rel_rmse = NA_real_,
    edge_rmse = NA_real_,
    gmds_stress = NA_real_,
    edge_energy = NA_real_,
    shape_rmse_vs_mds_edge_kk = NA_real_,
    elapsed_sec = elapsed,
    final_edge_rmse = NA_real_,
    final_scale = NA_real_,
    stringsAsFactors = FALSE
  )
}

run_schedule_on_fixture <- function(fixture, prepared, mds_coords, params) {
  tryCatch({
    level_count <- length(fixture$level_sizes)
    timed <- mek_timed(mek_layout_misf_edge_kk_reference(
      n = fixture$n,
      edges = fixture$edges,
      edge_lengths = fixture$edge_lengths,
      coords_post_insertion = fixture$coords_post_insertion,
      misf_order = fixture$misf_order,
      level_sizes = fixture$level_sizes,
      candidate_metric_pairs = fixture$candidate_metric_pairs,
      rho_schedule = seq(params$rho_start, 0, length.out = level_count),
      lambda_schedule = seq(params$lambda_start, 0, length.out = level_count),
      anchor_weight_schedule = seq(params$anchor_start, 0,
                                   length.out = level_count),
      exact_repulsion_below = 16L,
      repulsion_sample_count = 32L,
      repulsion_seed = 11L,
      level_max_iter = as.integer(params$level_max_iter),
      level_initial_step = 0.1,
      polish_max_iter = 30L,
      polish_initial_step = 0.1,
      polish_scale_mode = "profiled"
    ))
    fit <- timed$value
    score <- mek_score_baseline_coords(
      coords = fit$coords,
      prepared = prepared,
      reference_coords = mds_coords
    )
    data.frame(
      fixture = fixture$name,
      n = fixture$n,
      edge_count = nrow(fixture$edges),
      dim = fixture$dim,
      rho_start = params$rho_start,
      lambda_start = params$lambda_start,
      anchor_start = params$anchor_start,
      level_max_iter = params$level_max_iter,
      status = "ok",
      message = "",
      finite = score$finite,
      edge_rel_rmse = score$edge_rel_rmse,
      edge_rmse = score$edge_rmse,
      gmds_stress = score$gmds_stress,
      edge_energy = score$edge_energy,
      shape_rmse_vs_mds_edge_kk = score$shape_rmse_vs_mds_edge_kk,
      elapsed_sec = timed$elapsed,
      final_edge_rmse = fit$polish_trace$edge_rmse[[1L]],
      final_scale = fit$polish_trace$scale_value[[1L]],
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    empty_sweep_row(
      fixture = fixture,
      params = params,
      status = "error",
      message = conditionMessage(e)
    )
  })
}

reference_coords_for_fixture <- function(fixture, prepared) {
  tryCatch({
    mek_run_mds_edge_kk_method(fixture, prepared)$coords
  }, error = function(e) {
    warning(sprintf(
      "MDS+edge-KK reference failed for fixture %s: %s",
      fixture$name,
      conditionMessage(e)
    ), call. = FALSE)
    NULL
  })
}

fixtures <- mek_baseline_fixtures(dim = 2L)
fixture_results <- vector("list", length(fixtures))

for (fixture_index in seq_along(fixtures)) {
  fixture <- fixtures[[fixture_index]]
  prepared <- grip.prepare.graph.geodesic.mds(
    edges = fixture$edges,
    n = fixture$n,
    edge_weights = fixture$edge_lengths
  )
  mds_coords <- reference_coords_for_fixture(fixture, prepared)
  if (is.null(mds_coords)) {
    fixture_results[[fixture_index]] <- do.call(rbind, lapply(
      seq_len(nrow(parameter_grid)),
      function(i) empty_sweep_row(
        fixture,
        parameter_grid[i, ],
        status = "error",
        message = "MDS+edge-KK shape reference failed"
      )
    ))
    next
  }
  fixture_results[[fixture_index]] <- do.call(rbind, lapply(
    seq_len(nrow(parameter_grid)),
    function(i) run_schedule_on_fixture(
      fixture = fixture,
      prepared = prepared,
      mds_coords = mds_coords,
      params = parameter_grid[i, ]
    )
  ))
}

sweep_results <- do.call(rbind, fixture_results)
rownames(sweep_results) <- NULL

write.csv(
  sweep_results,
  file = file.path(out_dir, "schedule_sweep_results.csv"),
  row.names = FALSE
)

ok <- sweep_results[sweep_results$status == "ok", , drop = FALSE]
schedule_keys <- c("rho_start", "lambda_start", "anchor_start", "level_max_iter")
schedule_summary <- aggregate(
  cbind(edge_rel_rmse, gmds_stress, edge_energy,
        shape_rmse_vs_mds_edge_kk, elapsed_sec) ~
    rho_start + lambda_start + anchor_start + level_max_iter,
  data = ok,
  FUN = median
)
names(schedule_summary)[names(schedule_summary) == "edge_rel_rmse"] <- "median_edge_rel_rmse"
names(schedule_summary)[names(schedule_summary) == "gmds_stress"] <- "median_gmds_stress"
names(schedule_summary)[names(schedule_summary) == "edge_energy"] <- "median_edge_energy"
names(schedule_summary)[names(schedule_summary) == "shape_rmse_vs_mds_edge_kk"] <- "median_shape_rmse"
names(schedule_summary)[names(schedule_summary) == "elapsed_sec"] <- "median_elapsed_sec"

max_by_schedule <- aggregate(
  cbind(edge_rel_rmse, shape_rmse_vs_mds_edge_kk) ~
    rho_start + lambda_start + anchor_start + level_max_iter,
  data = ok,
  FUN = max
)
names(max_by_schedule)[names(max_by_schedule) == "edge_rel_rmse"] <- "max_edge_rel_rmse"
names(max_by_schedule)[names(max_by_schedule) == "shape_rmse_vs_mds_edge_kk"] <- "max_shape_rmse"
schedule_summary <- merge(
  schedule_summary,
  max_by_schedule,
  by = schedule_keys,
  sort = FALSE
)
schedule_summary$completed_fixtures <- as.integer(aggregate(
  fixture ~ rho_start + lambda_start + anchor_start + level_max_iter,
  data = ok,
  FUN = length
)$fixture)
schedule_summary$rank_edge <- rank(
  schedule_summary$median_edge_rel_rmse,
  ties.method = "min"
)
schedule_summary$rank_shape <- rank(
  schedule_summary$median_shape_rmse,
  ties.method = "min"
)
schedule_summary$rank_balanced <- rank(
  schedule_summary$rank_edge + schedule_summary$rank_shape,
  ties.method = "min"
)
schedule_summary <- schedule_summary[order(
  schedule_summary$rank_balanced,
  schedule_summary$rank_edge,
  schedule_summary$rank_shape
), ]
rownames(schedule_summary) <- NULL

write.csv(
  schedule_summary,
  file = file.path(out_dir, "schedule_sweep_summary.csv"),
  row.names = FALSE
)

best_edge <- schedule_summary[order(schedule_summary$rank_edge), ][seq_len(5L), ]
best_shape <- schedule_summary[order(schedule_summary$rank_shape), ][seq_len(5L), ]
best_balanced <- schedule_summary[order(schedule_summary$rank_balanced), ][seq_len(5L), ]

format_number <- function(x) {
  ifelse(is.na(x), "NA", formatC(x, digits = 4, format = "fg"))
}

schedule_line <- function(row) {
  sprintf(
    "| %s | %s | %s | %s | %s | %s | %s | %s |",
    format_number(as.numeric(row[["rho_start"]])),
    format_number(as.numeric(row[["lambda_start"]])),
    format_number(as.numeric(row[["anchor_start"]])),
    as.integer(row[["level_max_iter"]]),
    format_number(as.numeric(row[["median_edge_rel_rmse"]])),
    format_number(as.numeric(row[["median_shape_rmse"]])),
    format_number(as.numeric(row[["median_gmds_stress"]])),
    format_number(as.numeric(row[["median_elapsed_sec"]]))
  )
}

top_table <- function(df) {
  apply(df, 1L, schedule_line)
}

fixture_best <- do.call(rbind, lapply(split(ok, ok$fixture), function(df) {
  df[order(df$edge_rel_rmse, df$shape_rmse_vs_mds_edge_kk), ][1L, ]
}))
fixture_lines <- apply(fixture_best, 1L, function(row) {
  sprintf(
    "| %s | %s | %s | %s | %s | %s | %s |",
    row[["fixture"]],
    format_number(as.numeric(row[["rho_start"]])),
    format_number(as.numeric(row[["lambda_start"]])),
    format_number(as.numeric(row[["anchor_start"]])),
    as.integer(row[["level_max_iter"]]),
    format_number(as.numeric(row[["edge_rel_rmse"]])),
    format_number(as.numeric(row[["shape_rmse_vs_mds_edge_kk"]]))
  )
})

gate_results <- data.frame(
  gate = c(
    "MEK-SWEEP-001-complete",
    "MEK-SWEEP-002-finite",
    "MEK-SWEEP-003-shortlist"
  ),
  passed = c(
    all(sweep_results$status == "ok"),
    all(is.finite(ok$edge_rel_rmse)) &&
      all(is.finite(ok$shape_rmse_vs_mds_edge_kk)),
    nrow(schedule_summary) > 0L &&
      min(schedule_summary$median_edge_rel_rmse) < 0.001
  ),
  detail = c(
    "all schedule/fixture runs completed",
    "all successful runs have finite edge and shape metrics",
    "at least one schedule has median edge relative RMSE below 0.001"
  ),
  stringsAsFactors = FALSE
)
write.csv(
  gate_results,
  file = file.path(out_dir, "schedule_sweep_gates.csv"),
  row.names = FALSE
)

report_lines <- c(
  "# MISF Edge-KK Schedule Sweep",
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "This development-only sweep explores how the auxiliary metric, repulsion, and anchor schedules affect the current MISF edge-KK reference prototype on the tiny baseline fixtures.",
  "",
  "## Sweep Grid",
  "",
  sprintf("- Fixtures: %s", paste(vapply(fixtures, `[[`, character(1L), "name"), collapse = ", ")),
  "- `rho_start`: 0, 0.2, 0.6, 1.0",
  "- `lambda_start`: 0, 0.03, 0.1, 0.2",
  "- `anchor_start`: 0, 0.1, 0.25",
  "- `level_max_iter`: 8, 20",
  "- All schedules decay linearly to zero at the final MISF level.",
  "",
  "## Gates",
  "",
  "| Gate | Passed | Detail |",
  "|---|:---:|---|",
  apply(gate_results, 1L, function(row) {
    sprintf("| %s | %s | %s |",
            row[["gate"]],
            if (identical(row[["passed"]], "TRUE")) "yes" else "no",
            row[["detail"]])
  }),
  "",
  "## Top Balanced Schedules",
  "",
  "| rho | lambda | anchor | iter | Median Edge Rel RMSE | Median Shape RMSE | Median GMDS Stress | Median Sec |",
  "|---:|---:|---:|---:|---:|---:|---:|---:|",
  top_table(best_balanced),
  "",
  "## Best By Edge RMSE",
  "",
  "| rho | lambda | anchor | iter | Median Edge Rel RMSE | Median Shape RMSE | Median GMDS Stress | Median Sec |",
  "|---:|---:|---:|---:|---:|---:|---:|---:|",
  top_table(best_edge),
  "",
  "## Best By Shape RMSE",
  "",
  "| rho | lambda | anchor | iter | Median Edge Rel RMSE | Median Shape RMSE | Median GMDS Stress | Median Sec |",
  "|---:|---:|---:|---:|---:|---:|---:|---:|",
  top_table(best_shape),
  "",
  "## Best Per Fixture By Edge RMSE",
  "",
  "| Fixture | rho | lambda | anchor | iter | Edge Rel RMSE | Shape RMSE |",
  "|---|---:|---:|---:|---:|---:|---:|",
  fixture_lines,
  "",
  "## Interpretation",
  "",
  "This sweep is a diagnostic screen, not a default-selection study. If the best schedules cluster near low auxiliary terms, the current prototype is mostly behaving like final edge-KK polish from the insertion coordinates. If moderate auxiliary metric or repulsion improves shape without hurting edge RMSE, those regions are candidates for a second, narrower sweep on richer fixtures.",
  "",
  "## Output Files",
  "",
  "- `schedule_sweep_results.csv`",
  "- `schedule_sweep_summary.csv`",
  "- `schedule_sweep_gates.csv`",
  "- `schedule_sweep_report.md`"
)
writeLines(report_lines, con = file.path(out_dir, "schedule_sweep_report.md"))

if (!all(gate_results$passed)) {
  warning(
    sprintf("Schedule sweep gates not all passed: %s",
            paste(gate_results$gate[!gate_results$passed], collapse = ", ")),
    call. = FALSE
  )
}

message("MISF edge-KK schedule sweep complete.")
