#!/usr/bin/env Rscript

source("dev/misf-edge-kk/misf_edge_kk_objective_reference.R")
source("dev/misf-edge-kk/misf_edge_kk_constraint_reference.R")
source("dev/misf-edge-kk/misf_edge_kk_armijo_reference.R")
source("dev/misf-edge-kk/misf_edge_kk_level_reference.R")
source("dev/misf-edge-kk/misf_edge_kk_multilevel_reference.R")
source("dev/misf-edge-kk/misf_edge_kk_polish_reference.R")

out_dir <- "output/misf-edge-kk-validation"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

polish_gate <- function(gate, passed, detail = "", value = NA_real_,
                        tolerance = NA_real_) {
  data.frame(
    area = "MEK-POL",
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

run_polish_gates <- function() {
  fixture <- mek_level_fixture(dim = 3L)
  multilevel <- mek_refine_multilevel_reference(
    n = fixture$n,
    edges = fixture$edges,
    edge_lengths = fixture$edge_lengths,
    coords_post_insertion = fixture$coords_post_insertion,
    misf_order = fixture$misf_order,
    level_sizes = fixture$level_sizes,
    candidate_metric_pairs = fixture$candidate_metric_pairs,
    rho_schedule = c(0.8, 0.3, 0),
    lambda_schedule = c(0.12, 0.04, 0),
    anchor_weight_schedule = c(0.3, 0.1, 0),
    exact_repulsion_below = 3L,
    repulsion_sample_count = 4L,
    repulsion_seed = 77L,
    max_iter = 3L,
    initial_step = 0.1
  )
  polish <- mek_final_edge_kk_polish_reference(
    coords = multilevel$coords,
    edges = fixture$edges,
    edge_lengths = fixture$edge_lengths,
    max_iter = 5L,
    initial_step = 0.1,
    scale_mode = "profiled"
  )
  identity_polish <- mek_final_edge_kk_polish_reference(
    coords = multilevel$coords,
    edges = fixture$edges,
    edge_lengths = fixture$edge_lengths,
    max_iter = 1L,
    initial_step = 0.1,
    scale_mode = "identity"
  )

  trace <- polish$polish_trace
  energy_delta <- trace$total_energy - trace$initial_total_energy
  center_error <- max(abs(colMeans(polish$coords)))
  numeric_trace <- trace[
    setdiff(names(trace), c("stage", "scale_mode", "stop_reason", "accepted_step",
                            "metric_rmse"))
  ]

  rbind(
    polish_gate(
      "MEK-POL-001",
      trace$rho == 0 &&
        trace$lambda == 0 &&
        trace$metric_count == 0L &&
        trace$repulsion_count == 0L &&
        trace$anchor_weight_max == 0 &&
        trace$metric_energy == 0 &&
        trace$repulsion_energy == 0 &&
        trace$anchor_energy == 0,
      "final polish objective contains only original edge stress"
    ),
    polish_gate(
      "MEK-POL-002",
      identical(polish$active, seq_len(fixture$n)) &&
        trace$active_count == fixture$n &&
        trace$edge_count == nrow(fixture$edges),
      "final polish optimizes all vertices over original graph edges",
      value = trace$edge_count,
      tolerance = nrow(fixture$edges)
    ),
    polish_gate(
      "MEK-POL-003",
      energy_delta <= 1e-10,
      "final polish does not increase the profiled edge objective",
      value = energy_delta,
      tolerance = 1e-10
    ),
    polish_gate(
      "MEK-POL-004",
      all(is.finite(polish$coords)) &&
        finite_nonmissing(unlist(numeric_trace)) &&
        finite_nonmissing(polish$armijo_trace$energy_before) &&
        finite_nonmissing(polish$armijo_trace$energy_trial),
      "polish coordinates, trace values, and Armijo energies remain finite"
    ),
    polish_gate(
      "MEK-POL-005",
      identical(trace$scale_mode, "profiled") &&
        is.finite(trace$scale_value) &&
        trace$scale_value > 0 &&
        identical(identity_polish$polish_trace$scale_mode, "identity") &&
        identity_polish$polish_trace$scale_value == 1,
      "final polish marks scale mode and reports a valid scale value",
      value = trace$scale_value
    ),
    polish_gate(
      "MEK-POL-006",
      is.finite(trace$edge_rmse) &&
        is.na(trace$metric_rmse),
      "final edge RMSE is finite and auxiliary metric RMSE is absent"
    ),
    polish_gate(
      "MEK-POL-007",
      nrow(polish$metric_pairs) == 0L &&
        nrow(polish$repulsion_pairs) == 0L &&
        all(polish$anchor_weights == 0),
      "returned polish state contains no auxiliary metric, repulsion, or anchors"
    ),
    polish_gate(
      "MEK-POL-008",
      nrow(polish$armijo_trace) == trace$armijo_iterations &&
        trace$armijo_iterations > 0L &&
        trace$armijo_accepted >= 0L,
      "Armijo trace is retained and matches polish trace accounting",
      value = trace$armijo_iterations
    ),
    polish_gate(
      "MEK-POL-009",
      center_error <= 1e-12,
      "polish coordinates are recentered after accepted steps",
      value = center_error,
      tolerance = 1e-12
    )
  )
}

polish_results <- run_polish_gates()

write.csv(
  polish_results,
  file = file.path(out_dir, "polish_gate_results.csv"),
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
  polish_results
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
    polish = polish_results
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
  "- `initial_gate_results.csv`",
  "- `initial_gate_summary.csv`",
  "- `initial_gate_results.rds`",
  "",
  "## Scope",
  "",
  "This run covers objective values, objective gradients, total-objective composition, scale-policy guards, MISF active-set contracts, constraint construction, Armijo active-level refinement, one-level integrated refinement, multilevel orchestration, and final auxiliary-free edge-stress polish. It does not yet validate equivalence to the exported edge-KK optimizer."
)
writeLines(report_lines, con = file.path(out_dir, "initial_validation_report.md"))
writeLines(report_lines, con = file.path(out_dir, "initial_objective_validation_report.md"))

if (!all(combined_results$passed)) {
  failed <- combined_results$gate[!combined_results$passed]
  stop(sprintf("Initial validation failed: %s", paste(failed, collapse = ", ")),
       call. = FALSE)
}

message("Initial MISF edge-KK polish validation passed.")
