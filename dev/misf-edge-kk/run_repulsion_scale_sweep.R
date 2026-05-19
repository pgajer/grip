#!/usr/bin/env Rscript

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Package 'pkgload' is required for repulsion scale sweep", call. = FALSE)
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
  rho_start = c(0.1, 0.25),
  lambda_start = c(0, 0.002, 0.005, 0.01, 0.02, 0.05),
  repulsion_weight_mode = c("active_count", "sqrt_active", "unit"),
  anchor_start = 0.25,
  level_max_iter = 20L,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

empty_repulsion_row <- function(fixture, params, status, message,
                                elapsed = NA_real_) {
  data.frame(
    fixture = fixture$name,
    n = fixture$n,
    edge_count = nrow(fixture$edges),
    dim = fixture$dim,
    rho_start = params$rho_start,
    lambda_start = params$lambda_start,
    repulsion_weight_mode = params$repulsion_weight_mode,
    anchor_start = params$anchor_start,
    level_max_iter = params$level_max_iter,
    status = status,
    message = message,
    finite = FALSE,
    edge_rel_rmse = NA_real_,
    gmds_stress = NA_real_,
    edge_energy = NA_real_,
    shape_rmse_vs_mds_edge_kk = NA_real_,
    min_pair_distance = NA_real_,
    elapsed_sec = elapsed,
    final_scale = NA_real_,
    stringsAsFactors = FALSE
  )
}

run_repulsion_schedule <- function(fixture, prepared, mds_coords, params) {
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
      polish_scale_mode = "profiled",
      repulsion_weight_mode = params$repulsion_weight_mode
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
      repulsion_weight_mode = params$repulsion_weight_mode,
      anchor_start = params$anchor_start,
      level_max_iter = params$level_max_iter,
      status = "ok",
      message = "",
      finite = score$finite,
      edge_rel_rmse = score$edge_rel_rmse,
      gmds_stress = score$gmds_stress,
      edge_energy = score$edge_energy,
      shape_rmse_vs_mds_edge_kk = score$shape_rmse_vs_mds_edge_kk,
      min_pair_distance = mek_min_pair_distance(fit$coords),
      elapsed_sec = timed$elapsed,
      final_scale = fit$polish_trace$scale_value[[1L]],
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    empty_repulsion_row(
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
      function(i) empty_repulsion_row(
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
    function(i) run_repulsion_schedule(
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
  file = file.path(out_dir, "repulsion_scale_sweep_results.csv"),
  row.names = FALSE
)

ok <- sweep_results[sweep_results$status == "ok", , drop = FALSE]
summary <- aggregate(
  cbind(edge_rel_rmse, gmds_stress, edge_energy,
        shape_rmse_vs_mds_edge_kk, min_pair_distance, elapsed_sec) ~
    rho_start + lambda_start + repulsion_weight_mode,
  data = ok,
  FUN = median
)
names(summary)[names(summary) == "edge_rel_rmse"] <- "median_edge_rel_rmse"
names(summary)[names(summary) == "gmds_stress"] <- "median_gmds_stress"
names(summary)[names(summary) == "edge_energy"] <- "median_edge_energy"
names(summary)[names(summary) == "shape_rmse_vs_mds_edge_kk"] <- "median_shape_rmse"
names(summary)[names(summary) == "min_pair_distance"] <- "median_min_pair_distance"
names(summary)[names(summary) == "elapsed_sec"] <- "median_elapsed_sec"
summary$rank_edge <- rank(summary$median_edge_rel_rmse, ties.method = "min")
summary$rank_shape <- rank(summary$median_shape_rmse, ties.method = "min")
summary$rank_spacing <- rank(-summary$median_min_pair_distance,
                             ties.method = "min")
summary$rank_balanced <- rank(
  summary$rank_edge + summary$rank_shape,
  ties.method = "min"
)
summary <- summary[order(summary$rank_balanced, summary$rank_edge), ]
rownames(summary) <- NULL

write.csv(
  summary,
  file = file.path(out_dir, "repulsion_scale_sweep_summary.csv"),
  row.names = FALSE
)

lambda0 <- ok[ok$lambda_start == 0, c(
  "fixture", "rho_start", "repulsion_weight_mode", "edge_rel_rmse",
  "shape_rmse_vs_mds_edge_kk", "min_pair_distance"
)]
names(lambda0)[4:6] <- c("base_edge_rel_rmse", "base_shape_rmse",
                         "base_min_pair_distance")
with_base <- merge(
  ok,
  lambda0,
  by = c("fixture", "rho_start", "repulsion_weight_mode"),
  all.x = TRUE,
  sort = FALSE
)
with_base$edge_delta_vs_lambda0 <- with_base$edge_rel_rmse -
  with_base$base_edge_rel_rmse
with_base$shape_delta_vs_lambda0 <- with_base$shape_rmse_vs_mds_edge_kk -
  with_base$base_shape_rmse
with_base$min_pair_delta_vs_lambda0 <- with_base$min_pair_distance -
  with_base$base_min_pair_distance

delta_summary <- aggregate(
  cbind(edge_delta_vs_lambda0, shape_delta_vs_lambda0,
        min_pair_delta_vs_lambda0) ~ lambda_start + repulsion_weight_mode,
  data = with_base,
  FUN = median
)
write.csv(
  delta_summary,
  file = file.path(out_dir, "repulsion_scale_sweep_deltas.csv"),
  row.names = FALSE
)

format_number <- function(x) {
  ifelse(is.na(x), "NA", formatC(x, digits = 4, format = "fg"))
}

summary_line <- function(row) {
  sprintf(
    "| %s | %s | %s | %s | %s | %s | %s |",
    format_number(as.numeric(row[["rho_start"]])),
    format_number(as.numeric(row[["lambda_start"]])),
    row[["repulsion_weight_mode"]],
    format_number(as.numeric(row[["median_edge_rel_rmse"]])),
    format_number(as.numeric(row[["median_shape_rmse"]])),
    format_number(as.numeric(row[["median_min_pair_distance"]])),
    format_number(as.numeric(row[["median_elapsed_sec"]]))
  )
}

delta_line <- function(row) {
  sprintf(
    "| %s | %s | %s | %s | %s |",
    format_number(as.numeric(row[["lambda_start"]])),
    row[["repulsion_weight_mode"]],
    format_number(as.numeric(row[["edge_delta_vs_lambda0"]])),
    format_number(as.numeric(row[["shape_delta_vs_lambda0"]])),
    format_number(as.numeric(row[["min_pair_delta_vs_lambda0"]]))
  )
}

best_balanced <- summary[order(summary$rank_balanced), ][seq_len(5L), ]
best_edge <- summary[order(summary$rank_edge), ][seq_len(5L), ]
best_spacing <- summary[order(summary$rank_spacing), ][seq_len(5L), ]

gate_results <- data.frame(
  gate = c(
    "MEK-RSCALE-001-complete",
    "MEK-RSCALE-002-finite",
    "MEK-RSCALE-003-mode-coverage"
  ),
  passed = c(
    all(sweep_results$status == "ok"),
    all(is.finite(ok$edge_rel_rmse)) &&
      all(is.finite(ok$shape_rmse_vs_mds_edge_kk)) &&
      all(is.finite(ok$min_pair_distance)),
    setequal(unique(ok$repulsion_weight_mode),
             c("active_count", "sqrt_active", "unit"))
  ),
  detail = c(
    "all repulsion scale schedule/fixture runs completed",
    "all successful runs have finite edge, shape, and spacing metrics",
    "all requested repulsion weight modes were exercised"
  ),
  stringsAsFactors = FALSE
)
write.csv(
  gate_results,
  file = file.path(out_dir, "repulsion_scale_sweep_gates.csv"),
  row.names = FALSE
)

report_lines <- c(
  "# MISF Edge-KK Repulsion Scale Sweep",
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "This development-only sweep tests whether the earlier negative repulsion signal is caused by the current active-count weight normalization or by repulsion itself.",
  "",
  "## Sweep Grid",
  "",
  sprintf("- Fixtures: %s", paste(vapply(fixtures, `[[`, character(1L), "name"), collapse = ", ")),
  "- `rho_start`: 0.1, 0.25",
  "- `lambda_start`: 0, 0.002, 0.005, 0.01, 0.02, 0.05",
  "- `repulsion_weight_mode`: active_count, sqrt_active, unit",
  "- `anchor_start`: 0.25",
  "- `level_max_iter`: 20",
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
  "| rho | lambda | mode | Median Edge Rel RMSE | Median Shape RMSE | Median Min Pair Dist | Median Sec |",
  "|---:|---:|---|---:|---:|---:|---:|",
  apply(best_balanced, 1L, summary_line),
  "",
  "## Best By Edge RMSE",
  "",
  "| rho | lambda | mode | Median Edge Rel RMSE | Median Shape RMSE | Median Min Pair Dist | Median Sec |",
  "|---:|---:|---|---:|---:|---:|---:|",
  apply(best_edge, 1L, summary_line),
  "",
  "## Best By Minimum Pair Distance",
  "",
  "| rho | lambda | mode | Median Edge Rel RMSE | Median Shape RMSE | Median Min Pair Dist | Median Sec |",
  "|---:|---:|---|---:|---:|---:|---:|",
  apply(best_spacing, 1L, summary_line),
  "",
  "## Median Delta Versus Lambda Zero",
  "",
  "| lambda | mode | Edge Delta | Shape Delta | Min Pair Distance Delta |",
  "|---:|---|---:|---:|---:|",
  apply(delta_summary, 1L, delta_line),
  "",
  "## Interpretation",
  "",
  "A useful repulsion mode should increase minimum pair distance without materially worsening edge RMSE or shape RMSE. If all positive lambda settings have non-negative edge and shape deltas, then repulsion is not helping these fixtures even after weight rescaling.",
  "",
  "## Output Files",
  "",
  "- `repulsion_scale_sweep_results.csv`",
  "- `repulsion_scale_sweep_summary.csv`",
  "- `repulsion_scale_sweep_deltas.csv`",
  "- `repulsion_scale_sweep_gates.csv`",
  "- `repulsion_scale_sweep_report.md`"
)
writeLines(report_lines, con = file.path(out_dir, "repulsion_scale_sweep_report.md"))

if (!all(gate_results$passed)) {
  warning(
    sprintf("Repulsion scale sweep gates not all passed: %s",
            paste(gate_results$gate[!gate_results$passed], collapse = ", ")),
    call. = FALSE
  )
}

message("MISF edge-KK repulsion scale sweep complete.")
