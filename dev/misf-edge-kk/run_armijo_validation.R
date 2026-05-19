#!/usr/bin/env Rscript

source("dev/misf-edge-kk/misf_edge_kk_objective_reference.R")
source("dev/misf-edge-kk/misf_edge_kk_constraint_reference.R")
source("dev/misf-edge-kk/misf_edge_kk_armijo_reference.R")

out_dir <- "output/misf-edge-kk-validation"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

armijo_gate <- function(gate, passed, detail = "", value = NA_real_,
                        tolerance = NA_real_) {
  data.frame(
    area = "MEK-ARM",
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

run_fixed_objective_gate <- function() {
  fixture <- mek_armijo_fixture(dim = 3L)
  pair_keys <- mek_pair_key(fixture$repulsion_pairs)
  repulsion_weights <- fixture$repulsion_weights
  state_calls <- 0L
  checked_state_fn <- function(coords) {
    state_calls <<- state_calls + 1L
    if (!identical(pair_keys, mek_pair_key(fixture$repulsion_pairs))) {
      stop("repulsion pairs changed during line search", call. = FALSE)
    }
    if (!identical(repulsion_weights, fixture$repulsion_weights)) {
      stop("repulsion weights changed during line search", call. = FALSE)
    }
    fixture$state_fn(coords)
  }
  fit <- mek_armijo_refine(
    fixture$coords,
    fixture$active,
    checked_state_fn,
    max_iter = 2L,
    initial_step = 0.5
  )
  armijo_gate(
    "MEK-ARM-001",
    state_calls >= 2L &&
      identical(pair_keys, mek_pair_key(fixture$repulsion_pairs)) &&
      identical(repulsion_weights, fixture$repulsion_weights) &&
      nrow(fit$trace) > 0L,
    "objective ingredients remain fixed across Armijo trial evaluations",
    value = state_calls
  )
}

run_decrease_gate <- function() {
  fixture <- mek_armijo_fixture(dim = 3L)
  fit <- mek_armijo_refine(
    fixture$coords,
    fixture$active,
    fixture$state_fn,
    max_iter = 5L,
    initial_step = 0.5
  )
  accepted <- fit$trace[fit$trace$accepted, , drop = FALSE]
  errors <- if (nrow(accepted)) {
    accepted$energy_trial - accepted$armijo_target
  } else {
    numeric(0L)
  }
  armijo_gate(
    "MEK-ARM-002",
    nrow(accepted) > 0L &&
      all(accepted$armijo_satisfied) &&
      max(errors) <= 1e-12,
    "accepted trial steps satisfy the Armijo decrease inequality",
    value = if (length(errors)) max(errors) else NA_real_,
    tolerance = 1e-12
  )
}

run_shrink_gate <- function() {
  fixture <- mek_armijo_fixture(dim = 3L)
  fit <- mek_armijo_refine(
    fixture$coords,
    fixture$active,
    fixture$state_fn,
    max_iter = 1L,
    initial_step = 100,
    step_shrink = 0.25
  )
  trace <- fit$trace
  armijo_gate(
    "MEK-ARM-003",
    nrow(trace) == 1L && trace$accepted[[1L]] &&
      trace$shrink_attempts[[1L]] > 0L &&
      abs(trace$step[[1L]] - 100 * 0.25^trace$shrink_attempts[[1L]]) <= 1e-14,
    "rejected trials shrink the step size geometrically before acceptance",
    value = if (nrow(trace)) trace$shrink_attempts[[1L]] else NA_real_,
    tolerance = 0
  )
}

run_recenter_gate <- function() {
  fixture <- mek_armijo_fixture(dim = 3L)
  fit <- mek_armijo_refine(
    fixture$coords,
    fixture$active,
    fixture$state_fn,
    max_iter = 1L,
    initial_step = 0.5,
    recenter = TRUE
  )
  center_error <- max(abs(colMeans(fit$coords[fixture$active, , drop = FALSE])))
  armijo_gate(
    "MEK-ARM-004",
    nrow(fit$trace) > 0L && any(fit$trace$accepted) &&
      center_error <= 1e-12,
    "active coordinates are recentered after accepted steps",
    value = center_error,
    tolerance = 1e-12
  )
}

run_stopping_gate <- function() {
  fixture <- mek_armijo_fixture(dim = 2L)
  zero_state_fn <- function(coords) {
    list(
      energy = 0,
      gradient = matrix(0, nrow(coords), ncol(coords)),
      gradient_norm = 0
    )
  }
  budget_fit <- mek_armijo_refine(
    fixture$coords,
    fixture$active,
    fixture$state_fn,
    max_iter = 0L
  )
  gradient_fit <- mek_armijo_refine(
    fixture$coords,
    fixture$active,
    zero_state_fn,
    max_iter = 5L
  )
  rejected_fit <- mek_armijo_refine(
    fixture$coords,
    fixture$active,
    fixture$state_fn,
    max_iter = 1L,
    initial_step = 1e-12,
    min_step = 1e-6
  )
  observed <- c(
    budget_fit$stop_reason,
    gradient_fit$stop_reason,
    rejected_fit$stop_reason
  )
  expected <- c("iteration_budget", "gradient_tolerance", "no_accepted_step")
  armijo_gate(
    "MEK-ARM-005",
    identical(observed, expected),
    sprintf("stopping reasons recorded as %s", paste(observed, collapse = ", "))
  )
}

run_finite_stress_gate <- function() {
  fixture <- mek_armijo_fixture(dim = 3L, repeated = TRUE)
  fit <- mek_armijo_refine(
    fixture$coords,
    fixture$active,
    fixture$state_fn,
    max_iter = 3L,
    initial_step = 0.05,
    min_step = 1e-12
  )
  energies <- c(
    fit$trace$energy_before,
    fit$trace$energy_trial,
    fit$trace$energy_after_recenter,
    fit$final_state$energy
  )
  armijo_gate(
    "MEK-ARM-006",
    all(is.finite(fit$coords)) &&
      finite_nonmissing(energies) &&
      is.finite(fit$final_state$gradient_norm),
    "degenerate repeated-coordinate stress case remains finite"
  )
}

armijo_results <- do.call(rbind, list(
  run_fixed_objective_gate(),
  run_decrease_gate(),
  run_shrink_gate(),
  run_recenter_gate(),
  run_stopping_gate(),
  run_finite_stress_gate()
))

write.csv(
  armijo_results,
  file = file.path(out_dir, "armijo_gate_results.csv"),
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
  armijo_results
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
    armijo = armijo_results
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
  "- `initial_gate_results.csv`",
  "- `initial_gate_summary.csv`",
  "- `initial_gate_results.rds`",
  "",
  "## Scope",
  "",
  "This run covers objective values, objective gradients, total-objective composition, scale-policy guards, MISF active-set contracts, constraint construction, and Armijo active-level refinement. It does not yet validate edge-KK equivalence."
)
writeLines(report_lines, con = file.path(out_dir, "initial_validation_report.md"))
writeLines(report_lines, con = file.path(out_dir, "initial_objective_validation_report.md"))

if (!all(combined_results$passed)) {
  failed <- combined_results$gate[!combined_results$passed]
  stop(sprintf("Initial validation failed: %s", paste(failed, collapse = ", ")),
       call. = FALSE)
}

message("Initial MISF edge-KK Armijo validation passed.")
