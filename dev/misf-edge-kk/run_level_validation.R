#!/usr/bin/env Rscript

source("dev/misf-edge-kk/misf_edge_kk_objective_reference.R")
source("dev/misf-edge-kk/misf_edge_kk_constraint_reference.R")
source("dev/misf-edge-kk/misf_edge_kk_armijo_reference.R")
source("dev/misf-edge-kk/misf_edge_kk_level_reference.R")

out_dir <- "output/misf-edge-kk-validation"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

level_gate <- function(gate, passed, detail = "", value = NA_real_,
                       tolerance = NA_real_) {
  data.frame(
    area = "MEK-LVL",
    gate = gate,
    passed = isTRUE(passed),
    detail = detail,
    value = as.double(value),
    tolerance = as.double(tolerance),
    stringsAsFactors = FALSE
  )
}

same_pair_set <- function(x, y) {
  identical(sort(mek_pair_key(x)), sort(mek_pair_key(y)))
}

run_level_gates <- function() {
  fixture <- mek_level_fixture(dim = 3L)
  level2 <- mek_refine_level_reference(
    n = fixture$n,
    edges = fixture$edges,
    edge_lengths = fixture$edge_lengths,
    coords_post_insertion = fixture$coords_post_insertion,
    misf_order = fixture$misf_order,
    level_sizes = fixture$level_sizes,
    level_index = 2L,
    candidate_metric_pairs = fixture$candidate_metric_pairs,
    rho = 0.4,
    lambda = 0.08,
    anchor_weight = 0.2,
    exact_repulsion_below = 3L,
    repulsion_sample_count = 4L,
    repulsion_seed = 42L,
    max_iter = 4L,
    initial_step = 0.1
  )
  level2_repeat <- mek_refine_level_reference(
    n = fixture$n,
    edges = fixture$edges,
    edge_lengths = fixture$edge_lengths,
    coords_post_insertion = fixture$coords_post_insertion,
    misf_order = fixture$misf_order,
    level_sizes = fixture$level_sizes,
    level_index = 2L,
    candidate_metric_pairs = fixture$candidate_metric_pairs,
    rho = 0.4,
    lambda = 0.08,
    anchor_weight = 0.2,
    exact_repulsion_below = 3L,
    repulsion_sample_count = 4L,
    repulsion_seed = 42L,
    max_iter = 2L,
    initial_step = 0.1
  )
  level1 <- mek_refine_level_reference(
    n = fixture$n,
    edges = fixture$edges,
    edge_lengths = fixture$edge_lengths,
    coords_post_insertion = fixture$coords_post_insertion,
    misf_order = fixture$misf_order,
    level_sizes = fixture$level_sizes,
    level_index = 1L,
    candidate_metric_pairs = fixture$candidate_metric_pairs,
    rho = 0.5,
    lambda = 0.1,
    anchor_weight = 0.25,
    max_iter = 2L,
    initial_step = 0.1
  )

  required_fields <- c(
    "level", "active_count", "inserted_count", "edge_count",
    "metric_count", "repulsion_count", "rho", "lambda",
    "anchor_weight_mean", "anchor_weight_max", "scale_mode",
    "edge_scale", "metric_scale", "total_energy", "edge_energy",
    "metric_energy", "repulsion_energy", "anchor_energy",
    "gradient_norm", "accepted_step", "edge_rmse", "metric_rmse",
    "armijo_iterations", "armijo_accepted", "state_calls",
    "stop_reason", "elapsed_time"
  )
  trace <- level2$level_trace
  missing_fields <- setdiff(required_fields, names(trace))
  recomposed <- trace$edge_energy +
    trace$rho * trace$metric_energy +
    trace$lambda * trace$repulsion_energy +
    trace$anchor_energy
  recomposition_error <- abs(trace$total_energy - recomposed)

  inactive <- setdiff(seq_len(fixture$n), level2$active)
  inactive_absent <- !any(level2$active_edges[, 1L] %in% inactive |
                            level2$active_edges[, 2L] %in% inactive) &&
    !any(level2$metric_constraints$i %in% inactive |
           level2$metric_constraints$j %in% inactive) &&
    !any(level2$repulsion_pairs[, 1L] %in% inactive |
           level2$repulsion_pairs[, 2L] %in% inactive) &&
    !any(level2$anchors$vertices %in% inactive)
  inactive_coord_error <- max(abs(
    level2$coords[inactive, , drop = FALSE] -
      fixture$coords_post_insertion[inactive, , drop = FALSE]
  ))
  active_center_error <- max(abs(colMeans(
    level2$coords[level2$active, , drop = FALSE]
  )))
  finite_trace <- all(is.finite(unlist(trace[setdiff(
    names(trace),
    c("scale_mode", "repulsion_weight_mode", "stop_reason")
  )])))
  expected_inserted <- mek_misf_new_vertices(
    fixture$misf_order,
    fixture$level_sizes,
    2L
  )

  rbind(
    level_gate(
      "MEK-LVL-001",
      length(missing_fields) == 0L,
      sprintf("required level trace fields present; missing: %s",
              if (length(missing_fields)) paste(missing_fields, collapse = ", ") else "none")
    ),
    level_gate(
      "MEK-LVL-002",
      recomposition_error <= 1e-10,
      "level total energy equals edge + rho metric + lambda repulsion + anchors",
      value = recomposition_error,
      tolerance = 1e-10
    ),
    level_gate(
      "MEK-LVL-003",
      inactive_absent,
      "inactive vertices are absent from all assembled level constraints"
    ),
    level_gate(
      "MEK-LVL-004",
      identical(sort(level2$inserted), sort(expected_inserted)) &&
        trace$inserted_count == length(expected_inserted),
      "inserted vertices equal A_q minus A_{q+1} and are counted in trace",
      value = trace$inserted_count,
      tolerance = length(expected_inserted)
    ),
    level_gate(
      "MEK-LVL-005",
      nrow(level2$repulsion_pairs) == 4L &&
        same_pair_set(level2$repulsion_pairs, level2_repeat$repulsion_pairs),
      "sampled active repulsion pairs are deterministic and fixed for the level",
      value = nrow(level2$repulsion_pairs),
      tolerance = 4
    ),
    level_gate(
      "MEK-LVL-006",
      inactive_coord_error <= 1e-12,
      "inactive vertex coordinates remain unchanged by active-level refinement",
      value = inactive_coord_error,
      tolerance = 1e-12
    ),
    level_gate(
      "MEK-LVL-007",
      active_center_error <= 1e-12,
      "active coordinates are centered after level refinement",
      value = active_center_error,
      tolerance = 1e-12
    ),
    level_gate(
      "MEK-LVL-008",
      finite_trace && all(is.finite(level2$coords)) &&
        trace$armijo_iterations > 0L &&
        trace$state_calls >= trace$armijo_iterations,
      "level trace, coordinates, and Armijo accounting are finite and populated",
      value = trace$state_calls,
      tolerance = trace$armijo_iterations
    ),
    level_gate(
      "MEK-LVL-009",
      level1$level_trace$active_count < level2$level_trace$active_count,
      "deterministic fixture exercises increasing active-set levels",
      value = level2$level_trace$active_count - level1$level_trace$active_count,
      tolerance = 1
    )
  )
}

level_results <- run_level_gates()

write.csv(
  level_results,
  file = file.path(out_dir, "level_gate_results.csv"),
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
  level_results
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
    level = level_results
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
  "- `initial_gate_results.csv`",
  "- `initial_gate_summary.csv`",
  "- `initial_gate_results.rds`",
  "",
  "## Scope",
  "",
  "This run covers objective values, objective gradients, total-objective composition, scale-policy guards, MISF active-set contracts, constraint construction, Armijo active-level refinement, and one-level integrated refinement. It does not yet validate edge-KK equivalence."
)
writeLines(report_lines, con = file.path(out_dir, "initial_validation_report.md"))
writeLines(report_lines, con = file.path(out_dir, "initial_objective_validation_report.md"))

if (!all(combined_results$passed)) {
  failed <- combined_results$gate[!combined_results$passed]
  stop(sprintf("Initial validation failed: %s", paste(failed, collapse = ", ")),
       call. = FALSE)
}

message("Initial MISF edge-KK level validation passed.")
