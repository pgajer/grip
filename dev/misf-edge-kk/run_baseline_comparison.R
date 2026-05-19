#!/usr/bin/env Rscript

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Package 'pkgload' is required for baseline comparison", call. = FALSE)
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

empty_result_row <- function(fixture_name, method, status, message,
                             elapsed = NA_real_) {
  data.frame(
    fixture = fixture_name,
    method = method,
    status = status,
    message = message,
    finite = FALSE,
    edge_rel_rmse = NA_real_,
    edge_rmse = NA_real_,
    gmds_stress = NA_real_,
    edge_energy = NA_real_,
    shape_rmse_vs_mds_edge_kk = NA_real_,
    elapsed_sec = elapsed,
    stringsAsFactors = FALSE
  )
}

result_row <- function(fixture_name, method, coords, prepared,
                       reference_coords, elapsed, message = "") {
  score <- mek_score_baseline_coords(
    coords = coords,
    prepared = prepared,
    reference_coords = reference_coords
  )
  data.frame(
    fixture = fixture_name,
    method = method,
    status = "ok",
    message = message,
    finite = score$finite,
    edge_rel_rmse = score$edge_rel_rmse,
    edge_rmse = score$edge_rmse,
    gmds_stress = score$gmds_stress,
    edge_energy = score$edge_energy,
    shape_rmse_vs_mds_edge_kk = score$shape_rmse_vs_mds_edge_kk,
    elapsed_sec = elapsed,
    stringsAsFactors = FALSE
  )
}

run_method_row <- function(fixture, prepared, method, method_fn,
                           reference_coords = NULL) {
  tryCatch({
    timed <- mek_timed(method_fn())
    coords <- timed$value$coords
    result_row(
      fixture_name = fixture$name,
      method = method,
      coords = coords,
      prepared = prepared,
      reference_coords = reference_coords,
      elapsed = timed$elapsed
    )
  }, error = function(e) {
    empty_result_row(
      fixture_name = fixture$name,
      method = method,
      status = "error",
      message = conditionMessage(e)
    )
  })
}

run_fixture_comparison <- function(fixture) {
  prepared <- grip.prepare.graph.geodesic.mds(
    edges = fixture$edges,
    n = fixture$n,
    edge_weights = fixture$edge_lengths
  )

  mds_coords <- NULL
  mds_timed <- tryCatch(
    mek_timed(mek_run_mds_edge_kk_method(fixture, prepared)),
    error = function(e) e
  )
  mds_row <- if (inherits(mds_timed, "error")) {
    empty_result_row(
      fixture_name = fixture$name,
      method = "metric_mds_edge_kk",
      status = "error",
      message = conditionMessage(mds_timed)
    )
  } else {
    mds_coords <- mds_timed$value$coords
    result_row(
      fixture_name = fixture$name,
      method = "metric_mds_edge_kk",
      coords = mds_coords,
      prepared = prepared,
      reference_coords = mds_coords,
      elapsed = mds_timed$elapsed,
      message = "shape reference"
    )
  }

  reference_row <- run_method_row(
    fixture = fixture,
    prepared = prepared,
    method = "misf_edge_kk_reference",
    method_fn = function() mek_run_reference_method(fixture),
    reference_coords = mds_coords
  )
  weighted_row <- run_method_row(
    fixture = fixture,
    prepared = prepared,
    method = "weighted_grip_edge_kk",
    method_fn = function() mek_run_weighted_edge_kk_method(fixture, prepared),
    reference_coords = mds_coords
  )

  rbind(mds_row, reference_row, weighted_row)
}

fixtures <- mek_baseline_fixtures(dim = 2L)
results <- do.call(rbind, lapply(fixtures, run_fixture_comparison))
rownames(results) <- NULL

write.csv(
  results,
  file = file.path(out_dir, "baseline_comparison_results.csv"),
  row.names = FALSE
)

gate_results <- data.frame(
  gate = c(
    "MEK-BL-001-finite",
    "MEK-BL-002-no-errors",
    "MEK-BL-003-edge-rmse-finite",
    "MEK-BL-004-reference-edge-rmse-reasonable"
  ),
  passed = c(
    all(results$finite[results$status == "ok"]),
    all(results$status == "ok"),
    all(is.finite(results$edge_rel_rmse[results$status == "ok"])),
    all(results$edge_rel_rmse[
      results$status == "ok" & results$method == "misf_edge_kk_reference"
    ] < 0.75)
  ),
  detail = c(
    "all successful methods returned finite coordinates",
    "all fixture/method combinations completed without errors",
    "all successful methods reported finite edge relative RMSE",
    "reference prototype edge relative RMSE is below the initial smoke threshold"
  ),
  stringsAsFactors = FALSE
)
write.csv(
  gate_results,
  file = file.path(out_dir, "baseline_comparison_gates.csv"),
  row.names = FALSE
)

summary_table <- aggregate(
  cbind(edge_rel_rmse, gmds_stress, edge_energy, elapsed_sec) ~ method,
  data = results[results$status == "ok", ],
  FUN = median
)

format_number <- function(x) {
  ifelse(is.na(x), "NA", formatC(x, digits = 4, format = "fg"))
}

result_lines <- apply(results, 1L, function(row) {
  sprintf(
    "| %s | %s | %s | %s | %s | %s | %s |",
    row[["fixture"]],
    row[["method"]],
    row[["status"]],
    format_number(as.numeric(row[["edge_rel_rmse"]])),
    format_number(as.numeric(row[["gmds_stress"]])),
    format_number(as.numeric(row[["shape_rmse_vs_mds_edge_kk"]])),
    format_number(as.numeric(row[["elapsed_sec"]]))
  )
})

summary_lines <- apply(summary_table, 1L, function(row) {
  sprintf(
    "| %s | %s | %s | %s | %s |",
    row[["method"]],
    format_number(as.numeric(row[["edge_rel_rmse"]])),
    format_number(as.numeric(row[["gmds_stress"]])),
    format_number(as.numeric(row[["edge_energy"]])),
    format_number(as.numeric(row[["elapsed_sec"]]))
  )
})

gate_lines <- apply(gate_results, 1L, function(row) {
  sprintf(
    "| %s | %s | %s |",
    row[["gate"]],
    if (identical(row[["passed"]], "TRUE")) "yes" else "no",
    row[["detail"]]
  )
})

report_lines <- c(
  "# MISF Edge-KK Baseline Comparison",
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "This development-only smoke comparison checks whether the end-to-end MISF edge-KK reference prototype behaves sensibly on tiny deterministic fixtures relative to the current package workflows.",
  "",
  "## Methods",
  "",
  "- `misf_edge_kk_reference`: dev reference wrapper with multilevel MISF refinement plus final profiled edge-stress polish.",
  "- `weighted_grip_edge_kk`: `grip.layout.weighted()` followed by `grip.optimize.edge.kk.layout()`.",
  "- `metric_mds_edge_kk`: `grip.metric.mds.layout()` followed by `grip.optimize.edge.kk.layout()`; also used as the shape-reference embedding.",
  "",
  "## Smoke Gates",
  "",
  "| Gate | Passed | Detail |",
  "|---|:---:|---|",
  gate_lines,
  "",
  "## Median Metrics By Method",
  "",
  "| Method | Edge Rel RMSE | GMDS Stress | Edge Energy | Elapsed Sec |",
  "|---|---:|---:|---:|---:|",
  summary_lines,
  "",
  "## Fixture Results",
  "",
  "| Fixture | Method | Status | Edge Rel RMSE | GMDS Stress | Shape RMSE vs MDS+edge-KK | Elapsed Sec |",
  "|---|---|---|---:|---:|---:|---:|",
  result_lines,
  "",
  "## Interpretation",
  "",
  "This is not yet a performance claim. The intended gate at this stage is simple: complete on small graphs, produce finite coordinates and finite edge-stress diagnostics, and avoid obvious pathological edge distortion. Broader quality comparisons should come after adding richer fixtures and matching production optimizer settings more closely.",
  "",
  "## Output Files",
  "",
  "- `baseline_comparison_results.csv`",
  "- `baseline_comparison_gates.csv`",
  "- `baseline_comparison_report.md`"
)
writeLines(report_lines, con = file.path(out_dir, "baseline_comparison_report.md"))

if (!all(gate_results$passed)) {
  warning(
    sprintf("Baseline smoke gates not all passed: %s",
            paste(gate_results$gate[!gate_results$passed], collapse = ", ")),
    call. = FALSE
  )
}

message("MISF edge-KK baseline comparison complete.")
