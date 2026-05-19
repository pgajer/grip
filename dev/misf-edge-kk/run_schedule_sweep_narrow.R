#!/usr/bin/env Rscript

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Package 'pkgload' is required for narrow schedule sweep", call. = FALSE)
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
  rho_start = c(0, 0.1, 0.25, 0.5, 0.75, 1.0, 1.5),
  lambda_start = c(0, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3),
  anchor_start = 0.25,
  level_max_iter = 20L,
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
  file = file.path(out_dir, "schedule_sweep_narrow_results.csv"),
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
  file = file.path(out_dir, "schedule_sweep_narrow_summary.csv"),
  row.names = FALSE
)

rho_marginal <- aggregate(
  cbind(edge_rel_rmse, shape_rmse_vs_mds_edge_kk, gmds_stress) ~ rho_start,
  data = ok,
  FUN = median
)
rho_marginal$factor <- "rho"
names(rho_marginal)[names(rho_marginal) == "rho_start"] <- "value"

lambda_marginal <- aggregate(
  cbind(edge_rel_rmse, shape_rmse_vs_mds_edge_kk, gmds_stress) ~ lambda_start,
  data = ok,
  FUN = median
)
lambda_marginal$factor <- "lambda"
names(lambda_marginal)[names(lambda_marginal) == "lambda_start"] <- "value"

marginals <- rbind(
  rho_marginal[, c("factor", "value", "edge_rel_rmse",
                   "shape_rmse_vs_mds_edge_kk", "gmds_stress")],
  lambda_marginal[, c("factor", "value", "edge_rel_rmse",
                      "shape_rmse_vs_mds_edge_kk", "gmds_stress")]
)
names(marginals) <- c(
  "factor", "value", "median_edge_rel_rmse",
  "median_shape_rmse", "median_gmds_stress"
)
write.csv(
  marginals,
  file = file.path(out_dir, "schedule_sweep_narrow_marginals.csv"),
  row.names = FALSE
)

best_under_cap <- function(cap) {
  candidates <- schedule_summary[schedule_summary$median_shape_rmse <= cap, ]
  if (nrow(candidates) == 0L) {
    return(NULL)
  }
  candidates[order(candidates$median_edge_rel_rmse), ][1L, ]
}
shape_caps <- c(0.012, 0.015, 0.02)
cap_table <- do.call(rbind, lapply(shape_caps, function(cap) {
  row <- best_under_cap(cap)
  if (is.null(row)) {
    return(data.frame(
      shape_cap = cap,
      rho_start = NA_real_,
      lambda_start = NA_real_,
      anchor_start = NA_real_,
      level_max_iter = NA_integer_,
      median_edge_rel_rmse = NA_real_,
      median_shape_rmse = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    shape_cap = cap,
    rho_start = row$rho_start,
    lambda_start = row$lambda_start,
    anchor_start = row$anchor_start,
    level_max_iter = row$level_max_iter,
    median_edge_rel_rmse = row$median_edge_rel_rmse,
    median_shape_rmse = row$median_shape_rmse,
    stringsAsFactors = FALSE
  )
}))
write.csv(
  cap_table,
  file = file.path(out_dir, "schedule_sweep_narrow_shape_caps.csv"),
  row.names = FALSE
)

format_number <- function(x) {
  ifelse(is.na(x), "NA", formatC(x, digits = 4, format = "fg"))
}

schedule_line <- function(row) {
  sprintf(
    "| %s | %s | %s | %s | %s | %s |",
    format_number(as.numeric(row[["rho_start"]])),
    format_number(as.numeric(row[["lambda_start"]])),
    format_number(as.numeric(row[["median_edge_rel_rmse"]])),
    format_number(as.numeric(row[["median_shape_rmse"]])),
    format_number(as.numeric(row[["median_gmds_stress"]])),
    format_number(as.numeric(row[["median_elapsed_sec"]]))
  )
}

marginal_line <- function(row) {
  sprintf(
    "| %s | %s | %s | %s | %s |",
    row[["factor"]],
    format_number(as.numeric(row[["value"]])),
    format_number(as.numeric(row[["median_edge_rel_rmse"]])),
    format_number(as.numeric(row[["median_shape_rmse"]])),
    format_number(as.numeric(row[["median_gmds_stress"]]))
  )
}

cap_line <- function(row) {
  sprintf(
    "| %s | %s | %s | %s | %s |",
    format_number(as.numeric(row[["shape_cap"]])),
    format_number(as.numeric(row[["rho_start"]])),
    format_number(as.numeric(row[["lambda_start"]])),
    format_number(as.numeric(row[["median_edge_rel_rmse"]])),
    format_number(as.numeric(row[["median_shape_rmse"]]))
  )
}

best_edge <- schedule_summary[order(schedule_summary$rank_edge), ][seq_len(5L), ]
best_shape <- schedule_summary[order(schedule_summary$rank_shape), ][seq_len(5L), ]
best_balanced <- schedule_summary[order(schedule_summary$rank_balanced), ][seq_len(5L), ]

gate_results <- data.frame(
  gate = c(
    "MEK-NSWEEP-001-complete",
    "MEK-NSWEEP-002-finite",
    "MEK-NSWEEP-003-shape-cap"
  ),
  passed = c(
    all(sweep_results$status == "ok"),
    all(is.finite(ok$edge_rel_rmse)) &&
      all(is.finite(ok$shape_rmse_vs_mds_edge_kk)),
    any(schedule_summary$median_shape_rmse <= 0.012)
  ),
  detail = c(
    "all narrow sweep schedule/fixture runs completed",
    "all successful narrow sweep runs have finite edge and shape metrics",
    "at least one schedule satisfies median shape RMSE <= 0.012"
  ),
  stringsAsFactors = FALSE
)
write.csv(
  gate_results,
  file = file.path(out_dir, "schedule_sweep_narrow_gates.csv"),
  row.names = FALSE
)

report_lines <- c(
  "# MISF Edge-KK Narrow Schedule Sweep",
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "This development-only sweep fixes anchors and iteration budget, then varies auxiliary metric stress (`rho`) and repulsion (`lambda`) more finely to separate their effects.",
  "",
  "## Sweep Grid",
  "",
  sprintf("- Fixtures: %s", paste(vapply(fixtures, `[[`, character(1L), "name"), collapse = ", ")),
  "- `rho_start`: 0, 0.1, 0.25, 0.5, 0.75, 1.0, 1.5",
  "- `lambda_start`: 0, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3",
  "- `anchor_start`: 0.25",
  "- `level_max_iter`: 20",
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
  "| rho | lambda | Median Edge Rel RMSE | Median Shape RMSE | Median GMDS Stress | Median Sec |",
  "|---:|---:|---:|---:|---:|---:|",
  apply(best_balanced, 1L, schedule_line),
  "",
  "## Best By Edge RMSE",
  "",
  "| rho | lambda | Median Edge Rel RMSE | Median Shape RMSE | Median GMDS Stress | Median Sec |",
  "|---:|---:|---:|---:|---:|---:|",
  apply(best_edge, 1L, schedule_line),
  "",
  "## Best By Shape RMSE",
  "",
  "| rho | lambda | Median Edge Rel RMSE | Median Shape RMSE | Median GMDS Stress | Median Sec |",
  "|---:|---:|---:|---:|---:|---:|",
  apply(best_shape, 1L, schedule_line),
  "",
  "## Marginal Effects",
  "",
  "| Factor | Value | Median Edge Rel RMSE | Median Shape RMSE | Median GMDS Stress |",
  "|---|---:|---:|---:|---:|",
  apply(marginals, 1L, marginal_line),
  "",
  "## Best Edge RMSE Under Shape Caps",
  "",
  "| Shape Cap | rho | lambda | Median Edge Rel RMSE | Median Shape RMSE |",
  "|---:|---:|---:|---:|---:|",
  apply(cap_table, 1L, cap_line),
  "",
  "## Interpretation",
  "",
  "This sweep is designed to distinguish whether auxiliary metric stress or repulsion is responsible for the edge-RMSE gains seen in the broad sweep. The shape-cap table is the practical default-selection view: it asks how much edge accuracy can be recovered while keeping the embedding close to the MDS+edge-KK shape reference.",
  "",
  "## Output Files",
  "",
  "- `schedule_sweep_narrow_results.csv`",
  "- `schedule_sweep_narrow_summary.csv`",
  "- `schedule_sweep_narrow_marginals.csv`",
  "- `schedule_sweep_narrow_shape_caps.csv`",
  "- `schedule_sweep_narrow_gates.csv`",
  "- `schedule_sweep_narrow_report.md`"
)
writeLines(report_lines, con = file.path(out_dir, "schedule_sweep_narrow_report.md"))

if (!all(gate_results$passed)) {
  warning(
    sprintf("Narrow schedule sweep gates not all passed: %s",
            paste(gate_results$gate[!gate_results$passed], collapse = ", ")),
    call. = FALSE
  )
}

message("MISF edge-KK narrow schedule sweep complete.")
