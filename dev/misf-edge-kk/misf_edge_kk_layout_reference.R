# End-to-end reference wrapper for grip.layout.misf.edge.kk().
#
# This development-only function runs the full proposed skeleton:
# multilevel MISF refinement followed by final auxiliary-free edge-stress polish.

mek_combined_layout_trace <- function(multilevel, polish) {
  level_trace <- multilevel$level_trace
  level_trace$stage <- "misf_refinement"
  level_trace$scale_value <- level_trace$edge_scale
  level_trace$initial_total_energy <- NA_real_

  polish_trace <- polish$polish_trace
  polish_trace$level <- NA_integer_
  polish_trace$inserted_count <- 0L
  polish_trace$edge_scale <- polish_trace$scale_value
  polish_trace$metric_scale <- NA_real_

  trace_cols <- c(
    "stage", "level", "active_count", "inserted_count", "edge_count",
    "metric_count", "repulsion_count", "repulsion_weight_mode",
    "repulsion_weight_mass", "rho", "lambda",
    "anchor_weight_mean", "anchor_weight_max", "scale_mode", "scale_value",
    "edge_scale", "metric_scale", "initial_total_energy", "total_energy",
    "edge_energy", "metric_energy", "repulsion_energy", "anchor_energy",
    "gradient_norm", "accepted_step", "edge_rmse", "metric_rmse",
    "armijo_iterations", "armijo_accepted", "stop_reason", "elapsed_time"
  )
  for (nm in setdiff(trace_cols, names(level_trace))) {
    level_trace[[nm]] <- NA
  }
  for (nm in setdiff(trace_cols, names(polish_trace))) {
    polish_trace[[nm]] <- NA
  }
  out <- rbind(level_trace[, trace_cols], polish_trace[, trace_cols])
  rownames(out) <- NULL
  out
}

mek_layout_misf_edge_kk_reference <- function(n,
                                              edges,
                                              edge_lengths,
                                              coords_post_insertion,
                                              misf_order,
                                              level_sizes,
                                              candidate_metric_pairs,
                                              rho_schedule = NULL,
                                              lambda_schedule = NULL,
                                              anchor_weight_schedule = NULL,
                                              exact_repulsion_below = 64L,
                                              repulsion_sample_count = 128L,
                                              repulsion_seed = 1L,
                                              level_max_iter = 10L,
                                              level_initial_step = 0.5,
                                              polish_max_iter = 25L,
                                              polish_initial_step = 0.25,
                                              level_scale_mode = "identity",
                                              polish_scale_mode = "profiled",
                                              edge_scale = 1,
                                              metric_scale = 1,
                                              polish_fixed_scale = 1,
                                              repulsion_weight_mode = c("active_count",
                                                                        "unit",
                                                                        "sqrt_active")) {
  repulsion_weight_mode <- match.arg(repulsion_weight_mode)
  multilevel <- mek_refine_multilevel_reference(
    n = n,
    edges = edges,
    edge_lengths = edge_lengths,
    coords_post_insertion = coords_post_insertion,
    misf_order = misf_order,
    level_sizes = level_sizes,
    candidate_metric_pairs = candidate_metric_pairs,
    rho_schedule = rho_schedule,
    lambda_schedule = lambda_schedule,
    anchor_weight_schedule = anchor_weight_schedule,
    exact_repulsion_below = exact_repulsion_below,
    repulsion_sample_count = repulsion_sample_count,
    repulsion_seed = repulsion_seed,
    max_iter = level_max_iter,
    initial_step = level_initial_step,
    edge_scale = edge_scale,
    metric_scale = metric_scale,
    scale_mode = level_scale_mode,
    repulsion_weight_mode = repulsion_weight_mode
  )
  polish_input_coords <- multilevel$coords
  polish <- mek_final_edge_kk_polish_reference(
    coords = polish_input_coords,
    edges = edges,
    edge_lengths = edge_lengths,
    max_iter = polish_max_iter,
    initial_step = polish_initial_step,
    scale_mode = polish_scale_mode,
    fixed_scale = polish_fixed_scale
  )

  list(
    coords = polish$coords,
    trace = mek_combined_layout_trace(multilevel, polish),
    level_trace = multilevel$level_trace,
    polish_trace = polish$polish_trace,
    multilevel = multilevel,
    polish = polish,
    polish_input_coords = polish_input_coords,
    armijo_traces = c(multilevel$armijo_traces, list(final_polish = polish$armijo_trace))
  )
}
