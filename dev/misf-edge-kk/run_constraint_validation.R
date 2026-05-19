#!/usr/bin/env Rscript

source("dev/misf-edge-kk/misf_edge_kk_objective_reference.R")
source("dev/misf-edge-kk/misf_edge_kk_constraint_reference.R")

out_dir <- "output/misf-edge-kk-validation"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

contract_gate <- function(gate, passed, detail = "", value = NA_real_,
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

same_set <- function(x, y) {
  identical(sort(as.integer(x)), sort(as.integer(y)))
}

same_pairs <- function(x, y) {
  x <- mek_unique_pairs(x)
  y <- mek_unique_pairs(y)
  identical(sort(mek_pair_key(x)), sort(mek_pair_key(y)))
}

run_misf_gates <- function() {
  fixture <- mek_contract_fixture()
  active1 <- mek_misf_active_set(
    fixture$misf_order, fixture$level_sizes, 1L
  )
  active2 <- mek_misf_active_set(
    fixture$misf_order, fixture$level_sizes, 2L
  )
  active3 <- mek_misf_active_set(
    fixture$misf_order, fixture$level_sizes, 3L
  )
  new2 <- mek_misf_new_vertices(
    fixture$misf_order, fixture$level_sizes, 2L
  )
  new3 <- mek_misf_new_vertices(
    fixture$misf_order, fixture$level_sizes, 3L
  )

  active_edges2 <- mek_active_edges(fixture$edges, active2)
  metric2 <- mek_metric_constraints(
    n = fixture$n,
    edges = fixture$edges,
    edge_lengths = fixture$edge_lengths,
    active = active2,
    candidate_pairs = fixture$candidate_metric_pairs
  )
  repulsion2 <- mek_repulsion_pairs(active2, exact_below = 10L)
  anchors2 <- mek_capture_anchors(
    fixture$coords_post_insertion,
    active2,
    weights = seq_along(active2) / 10
  )
  inactive2 <- setdiff(seq_len(fixture$n), active2)

  coords <- fixture$coords_post_insertion
  state <- mek_edge_stress_state(
    coords = coords,
    edges = active_edges2,
    edge_lengths = rep(1, nrow(active_edges2)),
    stiffness = rep(1, nrow(active_edges2)),
    scale = 1,
    eps = 1e-8
  )
  previous_active_moved <- any(rowSums(abs(state$gradient[active1, , drop = FALSE])) > 0)

  rbind(
    contract_gate(
      "MEK-MISF-001-level-1",
      identical(active1, fixture$misf_order[seq_len(fixture$level_sizes[[1L]])]),
      "level 1 active set is exact MISF-order prefix"
    ),
    contract_gate(
      "MEK-MISF-001-level-2",
      identical(active2, fixture$misf_order[seq_len(fixture$level_sizes[[2L]])]),
      "level 2 active set is exact MISF-order prefix"
    ),
    contract_gate(
      "MEK-MISF-001-level-3",
      identical(active3, fixture$misf_order[seq_len(fixture$level_sizes[[3L]])]),
      "level 3 active set is exact MISF-order prefix"
    ),
    contract_gate(
      "MEK-MISF-002-level-2",
      same_set(new2, c(1L, 3L)),
      "new vertices equal A_2 minus A_1"
    ),
    contract_gate(
      "MEK-MISF-002-level-3",
      same_set(new3, c(5L, 7L)),
      "new vertices equal A_3 minus A_2"
    ),
    contract_gate(
      "MEK-MISF-003-active-edges",
      !any(active_edges2[, 1L] %in% inactive2 | active_edges2[, 2L] %in% inactive2),
      "inactive vertices excluded from original active edges"
    ),
    contract_gate(
      "MEK-MISF-003-metric",
      !any(metric2$i %in% inactive2 | metric2$j %in% inactive2),
      "inactive vertices excluded from metric constraints"
    ),
    contract_gate(
      "MEK-MISF-003-repulsion",
      !any(repulsion2[, 1L] %in% inactive2 | repulsion2[, 2L] %in% inactive2),
      "inactive vertices excluded from repulsion pairs"
    ),
    contract_gate(
      "MEK-MISF-003-anchors",
      !any(anchors2$vertices %in% inactive2),
      "inactive vertices excluded from anchors"
    ),
    contract_gate(
      "MEK-MISF-004",
      previous_active_moved,
      "previously active vertices can have nonzero active-level gradients"
    )
  )
}

run_constraint_gates <- function() {
  fixture <- mek_contract_fixture()
  active2 <- mek_misf_active_set(
    fixture$misf_order, fixture$level_sizes, 2L
  )
  active3 <- mek_misf_active_set(
    fixture$misf_order, fixture$level_sizes, 3L
  )
  active_edges2 <- mek_active_edges(fixture$edges, active2)
  expected_edges2 <- matrix(c(
    1L, 2L,
    1L, 3L,
    2L, 4L,
    4L, 6L
  ), ncol = 2L, byrow = TRUE)
  metric2 <- mek_metric_constraints(
    n = fixture$n,
    edges = fixture$edges,
    edge_lengths = fixture$edge_lengths,
    active = active2,
    candidate_pairs = fixture$candidate_metric_pairs
  )
  metric_keys <- paste(metric2$i, metric2$j, sep = "-")
  active_edge_keys <- mek_pair_key(active_edges2)
  dist <- mek_graph_distance_matrix(
    fixture$n, fixture$edges, fixture$edge_lengths
  )
  metric_expected_distances <- vapply(seq_len(nrow(metric2)), function(idx) {
    dist[metric2$i[[idx]], metric2$j[[idx]]]
  }, numeric(1L))

  exact_repulsion <- mek_repulsion_pairs(active2, exact_below = 10L)
  expected_exact_n <- length(active2) * (length(active2) - 1L) / 2L
  sampled1 <- mek_repulsion_pairs(
    active3, exact_below = 3L, sample_count = 4L, seed = 99L, level_index = 3L
  )
  sampled2 <- mek_repulsion_pairs(
    active3, exact_below = 3L, sample_count = 4L, seed = 99L, level_index = 3L
  )
  repulsion_weights <- mek_normalize_repulsion_weights(sampled1, length(active3))
  anchors <- mek_capture_anchors(fixture$coords_post_insertion, active2)
  mutated_coords <- fixture$coords_post_insertion
  mutated_coords[active2, ] <- mutated_coords[active2, ] + 100

  rbind(
    contract_gate(
      "MEK-CON-001",
      same_pairs(active_edges2, expected_edges2),
      "original active edge set contains all and only edges with endpoints in A_q"
    ),
    contract_gate(
      "MEK-CON-002-deduplicated",
      !any(duplicated(metric_keys)),
      "metric constraints are deduplicated"
    ),
    contract_gate(
      "MEK-CON-002-disjoint",
      !any(metric_keys %in% active_edge_keys),
      "metric constraints are disjoint from original active edges"
    ),
    contract_gate(
      "MEK-CON-003",
      max(abs(metric2$distance - metric_expected_distances)) <= 1e-12,
      "metric targets equal weighted graph-geodesic distances",
      value = max(abs(metric2$distance - metric_expected_distances)),
      tolerance = 1e-12
    ),
    contract_gate(
      "MEK-CON-004-exact",
      nrow(exact_repulsion) == expected_exact_n,
      "exact repulsion contains all unordered active pairs",
      value = nrow(exact_repulsion),
      tolerance = expected_exact_n
    ),
    contract_gate(
      "MEK-CON-004-sampled-deterministic",
      identical(sampled1, sampled2) && nrow(sampled1) == 4L,
      "sampled repulsion is deterministic for fixed seed and level"
    ),
    contract_gate(
      "MEK-CON-005",
      abs(sum(repulsion_weights) - length(active3)) <= 1e-12,
      "repulsion weights normalize total pair mass to active count",
      value = abs(sum(repulsion_weights) - length(active3)),
      tolerance = 1e-12
    ),
    contract_gate(
      "MEK-CON-006",
      identical(anchors$coords, fixture$coords_post_insertion[active2, , drop = FALSE]) &&
        !identical(anchors$coords, mutated_coords[active2, , drop = FALSE]),
      "anchor targets are copied from post-insertion coordinates"
    )
  )
}

normalize_existing <- function(df, area) {
  if (is.null(df) || nrow(df) == 0L) {
    return(data.frame())
  }
  out <- data.frame(
    area = area,
    gate = df$gate,
    passed = df$passed,
    detail = if ("detail" %in% names(df)) df$detail else "",
    value = if ("abs_error" %in% names(df)) df$abs_error else NA_real_,
    tolerance = if ("tolerance" %in% names(df)) df$tolerance else NA_real_,
    stringsAsFactors = FALSE
  )
  out
}

constraint_results <- rbind(run_misf_gates(), run_constraint_gates())
constraint_results$area <- sub("^((MEK-[A-Z]+).*)$", "\\2", constraint_results$gate)
constraint_results <- constraint_results[, c("area", "gate", "passed",
                                             "detail", "value", "tolerance")]

write.csv(
  constraint_results,
  file = file.path(out_dir, "constraint_gate_results.csv"),
  row.names = FALSE
)

objective_results <- if (file.exists(file.path(out_dir, "objective_gate_results.csv"))) {
  read.csv(file.path(out_dir, "objective_gate_results.csv"))
} else {
  NULL
}
scale_results <- if (file.exists(file.path(out_dir, "scale_policy_gate_results.csv"))) {
  read.csv(file.path(out_dir, "scale_policy_gate_results.csv"))
} else {
  NULL
}

combined_results <- rbind(
  normalize_existing(objective_results, "MEK-OBJ"),
  normalize_existing(scale_results, "MEK-SCL"),
  constraint_results
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
    objective = objective_results,
    scale = scale_results,
    constraint = constraint_results,
    all = combined_results
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
  "- `initial_gate_results.csv`",
  "- `initial_gate_summary.csv`",
  "- `initial_gate_results.rds`",
  "",
  "## Scope",
  "",
  "This run covers objective values, objective gradients, total-objective composition, scale-policy guards, MISF active-set contracts, and constraint construction. It does not yet validate Armijo refinement or edge-KK equivalence."
)
writeLines(report_lines, con = file.path(out_dir, "initial_validation_report.md"))
writeLines(report_lines, con = file.path(out_dir, "initial_objective_validation_report.md"))

if (!all(combined_results$passed)) {
  failed <- combined_results$gate[!combined_results$passed]
  stop(sprintf("Initial validation failed: %s", paste(failed, collapse = ", ")),
       call. = FALSE)
}

message("Initial MISF edge-KK constraint validation passed.")
