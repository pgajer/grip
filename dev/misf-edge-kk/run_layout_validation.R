#!/usr/bin/env Rscript

source("dev/misf-edge-kk/misf_edge_kk_objective_reference.R")
source("dev/misf-edge-kk/misf_edge_kk_constraint_reference.R")
source("dev/misf-edge-kk/misf_edge_kk_armijo_reference.R")
source("dev/misf-edge-kk/misf_edge_kk_level_reference.R")
source("dev/misf-edge-kk/misf_edge_kk_multilevel_reference.R")
source("dev/misf-edge-kk/misf_edge_kk_polish_reference.R")
source("dev/misf-edge-kk/misf_edge_kk_layout_reference.R")

out_dir <- "output/misf-edge-kk-validation"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

layout_gate <- function(gate, passed, detail = "", value = NA_real_,
                        tolerance = NA_real_) {
  data.frame(
    area = "MEK-LAYOUT",
    gate = gate,
    passed = isTRUE(passed),
    detail = detail,
    value = as.double(value),
    tolerance = as.double(tolerance),
    stringsAsFactors = FALSE
  )
}

finite_nonmissing <- function(x) {
  y <- x[!is.na(x)]
  length(y) == 0L || all(is.finite(y))
}

trace_without_timing <- function(trace) {
  trace[, setdiff(names(trace), "elapsed_time"), drop = FALSE]
}

run_layout_gates <- function() {
  fixture <- mek_level_fixture(dim = 3L)
  rho_schedule <- c(0.8, 0.3, 0)
  lambda_schedule <- c(0.12, 0.04, 0)
  anchor_schedule <- c(0.3, 0.1, 0)

  fit <- mek_layout_misf_edge_kk_reference(
    n = fixture$n,
    edges = fixture$edges,
    edge_lengths = fixture$edge_lengths,
    coords_post_insertion = fixture$coords_post_insertion,
    misf_order = fixture$misf_order,
    level_sizes = fixture$level_sizes,
    candidate_metric_pairs = fixture$candidate_metric_pairs,
    rho_schedule = rho_schedule,
    lambda_schedule = lambda_schedule,
    anchor_weight_schedule = anchor_schedule,
    exact_repulsion_below = 3L,
    repulsion_sample_count = 4L,
    repulsion_seed = 77L,
    level_max_iter = 3L,
    level_initial_step = 0.1,
    polish_max_iter = 5L,
    polish_initial_step = 0.1,
    polish_scale_mode = "profiled"
  )
  fit_repeat <- mek_layout_misf_edge_kk_reference(
    n = fixture$n,
    edges = fixture$edges,
    edge_lengths = fixture$edge_lengths,
    coords_post_insertion = fixture$coords_post_insertion,
    misf_order = fixture$misf_order,
    level_sizes = fixture$level_sizes,
    candidate_metric_pairs = fixture$candidate_metric_pairs,
    rho_schedule = rho_schedule,
    lambda_schedule = lambda_schedule,
    anchor_weight_schedule = anchor_schedule,
    exact_repulsion_below = 3L,
    repulsion_sample_count = 4L,
    repulsion_seed = 77L,
    level_max_iter = 3L,
    level_initial_step = 0.1,
    polish_max_iter = 5L,
    polish_initial_step = 0.1,
    polish_scale_mode = "profiled"
  )

  trace <- fit$trace
  final_row <- trace[nrow(trace), , drop = FALSE]
  expected_rows <- length(fixture$level_sizes) + 1L
  stage_ok <- identical(
    trace$stage,
    c(rep("misf_refinement", length(fixture$level_sizes)), "final_polish")
  )
  polish_start_error <- max(abs(fit$polish_input_coords - fit$multilevel$coords))
  numeric_trace <- trace[
    setdiff(names(trace), c("stage", "scale_mode", "stop_reason", "accepted_step",
                            "metric_rmse", "repulsion_weight_mode"))
  ]
  coords_repeat_error <- max(abs(fit$coords - fit_repeat$coords))
  trace_repeat_equal <- isTRUE(all.equal(
    trace_without_timing(fit$trace),
    trace_without_timing(fit_repeat$trace),
    tolerance = 1e-12,
    check.attributes = FALSE
  ))
  schedule_ends_off <- tail(fit$multilevel$rho_schedule, 1L) == 0 &&
    tail(fit$multilevel$lambda_schedule, 1L) == 0 &&
    tail(fit$multilevel$anchor_weight_schedule, 1L) == 0 &&
    final_row$rho == 0 &&
    final_row$lambda == 0 &&
    final_row$anchor_weight_max == 0

  rbind(
    layout_gate(
      "MEK-LAYOUT-001",
      all(is.finite(fit$coords)) &&
        finite_nonmissing(unlist(numeric_trace)),
      "complete pipeline returns finite coordinates and trace values"
    ),
    layout_gate(
      "MEK-LAYOUT-002",
      nrow(trace) == expected_rows && stage_ok,
      "combined trace has all MISF levels followed by one final polish row",
      value = nrow(trace),
      tolerance = expected_rows
    ),
    layout_gate(
      "MEK-LAYOUT-003",
      final_row$stage == "final_polish" &&
        final_row$metric_count == 0L &&
        final_row$repulsion_count == 0L &&
        final_row$rho == 0 &&
        final_row$lambda == 0 &&
        final_row$anchor_weight_max == 0,
      "final combined-trace row is auxiliary-free"
    ),
    layout_gate(
      "MEK-LAYOUT-004",
      polish_start_error <= 1e-12,
      "final polish starts from the multilevel output coordinates",
      value = polish_start_error,
      tolerance = 1e-12
    ),
    layout_gate(
      "MEK-LAYOUT-005",
      schedule_ends_off,
      "MISF schedules and final polish end with auxiliary terms off"
    ),
    layout_gate(
      "MEK-LAYOUT-006",
      is.finite(final_row$edge_rmse) &&
        identical(final_row$scale_mode, "profiled") &&
        is.finite(final_row$scale_value) &&
        final_row$scale_value > 0,
      "final edge RMSE and profiled scale are finite and reported",
      value = final_row$edge_rmse
    ),
    layout_gate(
      "MEK-LAYOUT-007",
      length(fit$armijo_traces) == expected_rows &&
        all(vapply(fit$armijo_traces, nrow, integer(1L)) ==
              trace$armijo_iterations),
      "one Armijo trace is retained for each stage row"
    ),
    layout_gate(
      "MEK-LAYOUT-008",
      coords_repeat_error <= 1e-12 && trace_repeat_equal,
      "running the full reference pipeline twice is deterministic",
      value = coords_repeat_error,
      tolerance = 1e-12
    )
  )
}

layout_results <- run_layout_gates()

write.csv(
  layout_results,
  file = file.path(out_dir, "layout_gate_results.csv"),
  row.names = FALSE
)

read_area <- function(path, area) {
  if (!file.exists(path)) {
    return(data.frame())
  }
  df <- read.csv(path)
  if (!("area" %in% names(df))) {
    df$area <- area
  }
  if (!("detail" %in% names(df))) {
    df$detail <- ""
  }
  if (!("value" %in% names(df))) {
    df$value <- if ("abs_error" %in% names(df)) df$abs_error else NA_real_
  }
  if (!("tolerance" %in% names(df))) {
    df$tolerance <- NA_real_
  }
  df <- df[, c("area", "gate", "passed", "detail", "value", "tolerance")]
  df$area <- area
  df
}

combined_results <- rbind(
  read_area(file.path(out_dir, "objective_gate_results.csv"), "MEK-OBJ"),
  read_area(file.path(out_dir, "scale_policy_gate_results.csv"), "MEK-SCL"),
  read_area(file.path(out_dir, "constraint_gate_results.csv"), "MEK-CON"),
  read_area(file.path(out_dir, "armijo_gate_results.csv"), "MEK-ARM"),
  read_area(file.path(out_dir, "level_gate_results.csv"), "MEK-LVL"),
  read_area(file.path(out_dir, "multilevel_gate_results.csv"), "MEK-MLVL"),
  read_area(file.path(out_dir, "polish_gate_results.csv"), "MEK-POL"),
  layout_results
)
combined_results$area <- sub("^((MEK-[A-Z]+).*)$", "\\2", combined_results$gate)
combined_results <- combined_results[order(combined_results$area, combined_results$gate), ]
rownames(combined_results) <- NULL

write.csv(
  combined_results,
  file = file.path(out_dir, "initial_gate_results.csv"),
  row.names = FALSE
)
saveRDS(
  list(
    all = combined_results,
    layout = layout_results
  ),
  file = file.path(out_dir, "initial_gate_results.rds")
)

summary_by_area <- aggregate(
  passed ~ area,
  data = combined_results,
  FUN = function(x) c(passed = sum(x), total = length(x))
)
summary_df <- data.frame(
  area = summary_by_area$area,
  passed = summary_by_area$passed[, "passed"],
  total = summary_by_area$passed[, "total"],
  stringsAsFactors = FALSE
)
summary_df$all_passed <- summary_df$passed == summary_df$total
write.csv(
  summary_df,
  file = file.path(out_dir, "initial_gate_summary.csv"),
  row.names = FALSE
)

report_lines <- c(
  "# MISF Edge-KK Initial Validation",
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "This development-only run validates the initial reference pieces for the proposed `grip.layout.misf.edge.kk()` algorithm.",
  "",
  "## Gate Summary",
  "",
  "| Area | Passed | Total | All Passed |",
  "|---|---:|---:|:---:|",
  apply(summary_df, 1L, function(row) {
    sprintf("| %s | %s | %s | %s |",
            row[["area"]], row[["passed"]], row[["total"]],
            if (identical(row[["all_passed"]], "TRUE")) "yes" else "no")
  }),
  "",
  "## Output Files",
  "",
  "- `objective_gate_results.csv`",
  "- `scale_policy_gate_results.csv`",
  "- `constraint_gate_results.csv`",
  "- `armijo_gate_results.csv`",
  "- `level_gate_results.csv`",
  "- `multilevel_gate_results.csv`",
  "- `polish_gate_results.csv`",
  "- `layout_gate_results.csv`",
  "- `initial_gate_results.csv`",
  "- `initial_gate_summary.csv`",
  "- `initial_gate_results.rds`",
  "",
  "## Scope",
  "",
  "This run covers objective values, objective gradients, total-objective composition, scale-policy guards, MISF active-set contracts, constraint construction, Armijo active-level refinement, one-level integrated refinement, multilevel orchestration, final auxiliary-free edge-stress polish, and the end-to-end reference wrapper. It does not yet validate equivalence to the exported edge-KK optimizer."
)
writeLines(report_lines, con = file.path(out_dir, "initial_validation_report.md"))
writeLines(report_lines, con = file.path(out_dir, "initial_objective_validation_report.md"))

if (!all(combined_results$passed)) {
  failed <- combined_results$gate[!combined_results$passed]
  stop(sprintf("Initial validation failed: %s", paste(failed, collapse = ", ")),
       call. = FALSE)
}

message("Initial MISF edge-KK layout wrapper validation passed.")
