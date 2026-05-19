#!/usr/bin/env Rscript

source("dev/misf-edge-kk/misf_edge_kk_objective_reference.R")
source("dev/misf-edge-kk/misf_edge_kk_constraint_reference.R")
source("dev/misf-edge-kk/misf_edge_kk_armijo_reference.R")
source("dev/misf-edge-kk/misf_edge_kk_level_reference.R")
source("dev/misf-edge-kk/misf_edge_kk_multilevel_reference.R")

out_dir <- "output/misf-edge-kk-validation"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

multilevel_gate <- function(gate, passed, detail = "", value = NA_real_,
                            tolerance = NA_real_) {
  data.frame(
    area = "MEK-MLVL",
    gate = gate,
    passed = isTRUE(passed),
    detail = detail,
    value = as.double(value),
    tolerance = as.double(tolerance),
    stringsAsFactors = FALSE
  )
}

all_level_coordinates_finite <- function(coords_by_level) {
  all(vapply(coords_by_level, function(coords) all(is.finite(coords)), logical(1L)))
}

inactive_unchanged_error <- function(fit, misf_order, level_sizes) {
  errors <- vapply(seq_along(level_sizes), function(level_index) {
    active <- mek_misf_active_set(misf_order, level_sizes, level_index)
    inactive <- setdiff(seq_len(nrow(fit$insertion_coords)), active)
    if (!length(inactive)) {
      return(0)
    }
    max(abs(
      fit$coords_by_level[[level_index]][inactive, , drop = FALSE] -
        fit$insertion_coords[inactive, , drop = FALSE]
    ))
  }, numeric(1L))
  max(errors)
}

carry_forward_error <- function(fit, misf_order, level_sizes) {
  if (length(level_sizes) < 2L) {
    return(0)
  }
  errors <- vapply(seq_len(length(level_sizes) - 1L) + 1L, function(level_index) {
    previous_active <- mek_misf_active_set(
      misf_order, level_sizes, level_index - 1L
    )
    anchors <- fit$level_results[[level_index]]$anchors
    idx <- match(previous_active, anchors$vertices)
    max(abs(
      anchors$coords[idx, , drop = FALSE] -
        fit$coords_by_level[[level_index - 1L]][previous_active, , drop = FALSE]
    ))
  }, numeric(1L))
  max(errors)
}

term_recomposition_error <- function(level_trace) {
  recomposed <- level_trace$edge_energy +
    level_trace$rho * level_trace$metric_energy +
    level_trace$lambda * level_trace$repulsion_energy +
    level_trace$anchor_energy
  max(abs(level_trace$total_energy - recomposed))
}

run_multilevel_gates <- function() {
  fixture <- mek_level_fixture(dim = 3L)
  rho_schedule <- c(0.8, 0.3, 0)
  lambda_schedule <- c(0.12, 0.04, 0)
  anchor_schedule <- c(0.3, 0.1, 0)
  fit <- mek_refine_multilevel_reference(
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
    max_iter = 3L,
    initial_step = 0.1
  )

  trace <- fit$level_trace
  inactive_error <- inactive_unchanged_error(
    fit, fixture$misf_order, fixture$level_sizes
  )
  carry_error <- carry_forward_error(fit, fixture$misf_order, fixture$level_sizes)
  recomposition_error <- term_recomposition_error(trace)
  numeric_trace <- trace[
    setdiff(names(trace), c("scale_mode", "repulsion_weight_mode",
                            "stop_reason", "accepted_step"))
  ]
  finite_trace <- all(is.finite(unlist(numeric_trace)))
  final_active <- mek_misf_active_set(
    fixture$misf_order,
    fixture$level_sizes,
    length(fixture$level_sizes)
  )

  rbind(
    multilevel_gate(
      "MEK-MLVL-001",
      identical(trace$active_count, as.integer(fixture$level_sizes)),
      "active counts follow misf.level_size at every level"
    ),
    multilevel_gate(
      "MEK-MLVL-002",
      identical(trace$level, seq_along(fixture$level_sizes)),
      "level traces are ordered and complete"
    ),
    multilevel_gate(
      "MEK-MLVL-003",
      identical(sort(final_active), seq_len(fixture$n)) &&
        trace$active_count[[nrow(trace)]] == fixture$n,
      "final level activates the full graph",
      value = trace$active_count[[nrow(trace)]],
      tolerance = fixture$n
    ),
    multilevel_gate(
      "MEK-MLVL-004",
      all_level_coordinates_finite(fit$coords_by_level) &&
        all(is.finite(fit$coords)) &&
        finite_trace,
      "coordinates and multilevel trace remain finite at every level"
    ),
    multilevel_gate(
      "MEK-MLVL-005",
      inactive_error <= 1e-12,
      "vertices remain at insertion coordinates until their active level",
      value = inactive_error,
      tolerance = 1e-12
    ),
    multilevel_gate(
      "MEK-MLVL-006",
      all.equal(trace$rho, rho_schedule, tolerance = 1e-12) == TRUE &&
        all.equal(trace$lambda, lambda_schedule, tolerance = 1e-12) == TRUE &&
        all.equal(trace$anchor_weight_mean, anchor_schedule,
                  tolerance = 1e-12) == TRUE,
      "rho, lambda, and anchor schedules are applied exactly"
    ),
    multilevel_gate(
      "MEK-MLVL-007",
      carry_error <= 1e-12,
      "previously active coordinates carry forward into the next level",
      value = carry_error,
      tolerance = 1e-12
    ),
    multilevel_gate(
      "MEK-MLVL-008",
      recomposition_error <= 1e-10,
      "each level total energy recomposes from objective terms",
      value = recomposition_error,
      tolerance = 1e-10
    ),
    multilevel_gate(
      "MEK-MLVL-009",
      length(fit$armijo_traces) == length(fixture$level_sizes) &&
        all(vapply(fit$armijo_traces, nrow, integer(1L)) ==
              trace$armijo_iterations),
      "one Armijo trace is retained per level and matches trace accounting"
    )
  )
}

multilevel_results <- run_multilevel_gates()

write.csv(
  multilevel_results,
  file = file.path(out_dir, "multilevel_gate_results.csv"),
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
  multilevel_results
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
    multilevel = multilevel_results
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
  "- `initial_gate_results.csv`",
  "- `initial_gate_summary.csv`",
  "- `initial_gate_results.rds`",
  "",
  "## Scope",
  "",
  "This run covers objective values, objective gradients, total-objective composition, scale-policy guards, MISF active-set contracts, constraint construction, Armijo active-level refinement, one-level integrated refinement, and multilevel orchestration. It does not yet validate edge-KK equivalence."
)
writeLines(report_lines, con = file.path(out_dir, "initial_validation_report.md"))
writeLines(report_lines, con = file.path(out_dir, "initial_objective_validation_report.md"))

if (!all(combined_results$passed)) {
  failed <- combined_results$gate[!combined_results$passed]
  stop(sprintf("Initial validation failed: %s", paste(failed, collapse = ", ")),
       call. = FALSE)
}

message("Initial MISF edge-KK multilevel validation passed.")
