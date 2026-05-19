#!/usr/bin/env Rscript

source("dev/misf-edge-kk/misf_edge_kk_objective_reference.R")

out_dir <- "output/misf-edge-kk-validation"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

gate_result <- function(gate, passed, detail = "", value = NA_real_,
                        tolerance = NA_real_) {
  data.frame(
    gate = gate,
    passed = isTRUE(passed),
    detail = detail,
    value = as.double(value),
    tolerance = as.double(tolerance),
    stringsAsFactors = FALSE
  )
}

run_value_gate <- function() {
  coords <- matrix(c(0, 0, 3, 4), nrow = 2L, byrow = TRUE)
  edges <- matrix(c(1L, 2L), ncol = 2L)
  state <- mek_edge_stress_state(
    coords = coords,
    edges = edges,
    edge_lengths = 4,
    stiffness = 2,
    scale = 1,
    eps = 0
  )
  expected <- 0.5 * 2 * (5 - 4)^2
  gate_result(
    gate = "MEK-OBJ-001",
    passed = abs(state$energy - expected) <= 1e-12,
    detail = "edge-stress value matches hand-computed two-vertex fixture",
    value = abs(state$energy - expected),
    tolerance = 1e-12
  )
}

run_gradient_gates <- function() {
  dims <- c(2L, 3L, 5L)
  rows <- list()
  idx <- 1L
  for (dim in dims) {
    fixture <- mek_objective_fixture(dim = dim)
    rows[[idx]] <- mek_gradient_check(
      name = sprintf("MEK-OBJ-002-dim-%d", dim),
      coords = fixture$coords,
      state_fn = function(x) mek_edge_stress_state(
        coords = x,
        edges = fixture$edges,
        edge_lengths = fixture$edge_lengths,
        stiffness = fixture$edge_stiffness,
        scale = 1,
        eps = 1e-8
      )
    )
    idx <- idx + 1L
    rows[[idx]] <- mek_gradient_check(
      name = sprintf("MEK-OBJ-003-dim-%d", dim),
      coords = fixture$coords,
      state_fn = function(x) mek_metric_stress_state(
        coords = x,
        metric_pairs = fixture$metric_pairs,
        metric_distances = fixture$metric_distances,
        stiffness = fixture$metric_stiffness,
        scale = 1,
        eps = 1e-8
      )
    )
    idx <- idx + 1L
    repulsion_fixture <- mek_objective_fixture(dim = dim, repeated = TRUE)
    rows[[idx]] <- mek_gradient_check(
      name = sprintf("MEK-OBJ-004-dim-%d", dim),
      coords = repulsion_fixture$coords,
      state_fn = function(x) mek_log_repulsion_state(
        coords = x,
        pairs = repulsion_fixture$repulsion_pairs,
        pair_weights = repulsion_fixture$repulsion_weights,
        eps = 1e-4
      ),
      tolerance = 1e-4
    )
    idx <- idx + 1L
    rows[[idx]] <- mek_gradient_check(
      name = sprintf("MEK-OBJ-005-dim-%d", dim),
      coords = fixture$coords,
      state_fn = function(x) mek_anchor_state(
        coords = x,
        anchor_targets = fixture$anchor_targets,
        anchor_weights = fixture$anchor_weights
      )
    )
    idx <- idx + 1L
    rows[[idx]] <- mek_gradient_check(
      name = sprintf("MEK-OBJ-006-dim-%d", dim),
      coords = fixture$coords,
      state_fn = function(x) mek_total_state(
        coords = x,
        edges = fixture$edges,
        edge_lengths = fixture$edge_lengths,
        edge_stiffness = fixture$edge_stiffness,
        metric_pairs = fixture$metric_pairs,
        metric_distances = fixture$metric_distances,
        metric_stiffness = fixture$metric_stiffness,
        repulsion_pairs = fixture$repulsion_pairs,
        repulsion_weights = fixture$repulsion_weights,
        anchor_targets = fixture$anchor_targets,
        anchor_weights = fixture$anchor_weights,
        rho = 0.25,
        lambda = 0.1,
        edge_scale = 1,
        metric_scale = 1,
        edge_eps = 1e-8,
        repulsion_eps = 1e-4
      ),
      tolerance = 1e-4
    )
    idx <- idx + 1L
  }
  do.call(rbind, rows)
}

run_composition_gate <- function() {
  fixture <- mek_objective_fixture(dim = 3L)
  total <- mek_total_state(
    coords = fixture$coords,
    edges = fixture$edges,
    edge_lengths = fixture$edge_lengths,
    edge_stiffness = fixture$edge_stiffness,
    metric_pairs = fixture$metric_pairs,
    metric_distances = fixture$metric_distances,
    metric_stiffness = fixture$metric_stiffness,
    repulsion_pairs = fixture$repulsion_pairs,
    repulsion_weights = fixture$repulsion_weights,
    anchor_targets = fixture$anchor_targets,
    anchor_weights = fixture$anchor_weights,
    rho = 0.25,
    lambda = 0.1,
    repulsion_eps = 1e-4
  )
  expected_energy <- total$edge$energy +
    0.25 * total$metric$energy +
    0.1 * total$repulsion$energy +
    total$anchor$energy
  expected_gradient <- total$edge$gradient +
    0.25 * total$metric$gradient +
    0.1 * total$repulsion$gradient +
    total$anchor$gradient
  data.frame(
    gate = c("MEK-OBJ-006-energy-composition",
             "MEK-OBJ-006-gradient-composition"),
    passed = c(
      abs(total$energy - expected_energy) <= 1e-12,
      max(abs(total$gradient - expected_gradient)) <= 1e-12
    ),
    abs_error = c(
      abs(total$energy - expected_energy),
      max(abs(total$gradient - expected_gradient))
    ),
    rel_error = c(0, 0),
    tolerance = c(1e-12, 1e-12),
    stringsAsFactors = FALSE
  )
}

run_scale_policy_gates <- function() {
  ok_identity <- tryCatch({
    mek_validate_misf_scale_policy("identity", rho = 0.1, lambda = 0.2)
  }, error = function(e) FALSE)
  ok_fixed <- tryCatch({
    mek_validate_misf_scale_policy("global_fixed", rho = 0.1, lambda = 0.2)
  }, error = function(e) FALSE)
  bad_profiled <- tryCatch({
    mek_validate_misf_scale_policy("profiled", rho = 0.1, lambda = 0)
    FALSE
  }, error = function(e) TRUE)
  data.frame(
    gate = c("MEK-SCL-001-identity", "MEK-SCL-001-global-fixed",
             "MEK-SCL-001-profiled-rejected"),
    passed = c(ok_identity, ok_fixed, bad_profiled),
    detail = c(
      "identity scale allowed with active regularization",
      "global fixed scale allowed with active regularization",
      "profiled scale rejected when rho > 0 or lambda > 0"
    ),
    stringsAsFactors = FALSE
  )
}

value_results <- run_value_gate()
gradient_results <- run_gradient_gates()
composition_results <- run_composition_gate()
scale_results <- run_scale_policy_gates()

objective_results <- rbind(
  data.frame(
    gate = value_results$gate,
    passed = value_results$passed,
    abs_error = value_results$value,
    rel_error = value_results$value,
    tolerance = value_results$tolerance,
    stringsAsFactors = FALSE
  ),
  gradient_results,
  composition_results
)

scale_results_out <- data.frame(
  gate = scale_results$gate,
  passed = scale_results$passed,
  abs_error = NA_real_,
  rel_error = NA_real_,
  tolerance = NA_real_,
  stringsAsFactors = FALSE
)

all_results <- rbind(objective_results, scale_results_out)
all_results$area <- sub("^((MEK-[A-Z]+).*)$", "\\2", all_results$gate)
all_results <- all_results[, c("area", "gate", "passed", "abs_error",
                               "rel_error", "tolerance")]

write.csv(
  objective_results,
  file = file.path(out_dir, "objective_gate_results.csv"),
  row.names = FALSE
)
write.csv(
  scale_results,
  file = file.path(out_dir, "scale_policy_gate_results.csv"),
  row.names = FALSE
)
write.csv(
  all_results,
  file = file.path(out_dir, "initial_gate_results.csv"),
  row.names = FALSE
)
saveRDS(
  list(
    objective = objective_results,
    scale = scale_results,
    all = all_results
  ),
  file = file.path(out_dir, "initial_gate_results.rds")
)

summary_by_area <- aggregate(
  passed ~ area,
  data = all_results,
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
  "# MISF Edge-KK Initial Objective Validation",
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "This development-only run validates the reference R objective pieces for the proposed `grip.layout.misf.edge.kk()` algorithm.",
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
  "- `initial_gate_results.csv`",
  "- `initial_gate_summary.csv`",
  "- `initial_gate_results.rds`",
  "",
  "## Scope",
  "",
  "This run covers the first prototype milestone only: objective values, objective gradients, total-objective composition, and the cheap scale-policy guard. It does not yet validate MISF active sets, constraint construction, Armijo refinement, or edge-KK equivalence."
)
writeLines(report_lines, con = file.path(out_dir, "initial_objective_validation_report.md"))

if (!all(all_results$passed)) {
  failed <- all_results$gate[!all_results$passed]
  stop(sprintf("Initial validation failed: %s", paste(failed, collapse = ", ")),
       call. = FALSE)
}

message("Initial MISF edge-KK objective validation passed.")
