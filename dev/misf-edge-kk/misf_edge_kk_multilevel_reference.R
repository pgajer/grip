# Multilevel reference orchestration for grip.layout.misf.edge.kk().
#
# This development-only layer loops over MISF levels, carries coordinates
# forward, applies simple deterministic insertion coordinates, and collects
# level and Armijo traces.

mek_schedule_values <- function(x, n, name) {
  x <- as.double(x)
  if (length(x) == 1L) {
    x <- rep(x, n)
  }
  if (length(x) != n || any(!is.finite(x)) || any(x < 0)) {
    stop(sprintf("%s must be a finite non-negative scalar or length-n vector",
                 name),
         call. = FALSE)
  }
  x
}

mek_refine_multilevel_reference <- function(n,
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
                                            max_iter = 10L,
                                            initial_step = 0.5,
                                            edge_scale = 1,
                                            metric_scale = 1,
                                            scale_mode = "identity",
                                            repulsion_weight_mode = c("active_count",
                                                                      "unit",
                                                                      "sqrt_active")) {
  coords_post_insertion <- mek_validate_coords(coords_post_insertion)
  if (nrow(coords_post_insertion) != n) {
    stop("coords_post_insertion row count must equal n", call. = FALSE)
  }
  repulsion_weight_mode <- match.arg(repulsion_weight_mode)
  level_count <- length(level_sizes)
  if (level_count < 1L) {
    stop("level_sizes must contain at least one level", call. = FALSE)
  }
  if (is.null(rho_schedule)) {
    rho_schedule <- seq(0.5, 0, length.out = level_count)
  }
  if (is.null(lambda_schedule)) {
    lambda_schedule <- seq(0.1, 0, length.out = level_count)
  }
  if (is.null(anchor_weight_schedule)) {
    anchor_weight_schedule <- seq(0.25, 0, length.out = level_count)
  }
  rho_schedule <- mek_schedule_values(rho_schedule, level_count, "rho_schedule")
  lambda_schedule <- mek_schedule_values(
    lambda_schedule, level_count, "lambda_schedule"
  )
  anchor_weight_schedule <- mek_schedule_values(
    anchor_weight_schedule, level_count, "anchor_weight_schedule"
  )

  current <- coords_post_insertion
  level_results <- vector("list", level_count)
  coords_by_level <- vector("list", level_count)

  for (level_index in seq_len(level_count)) {
    inserted <- mek_misf_new_vertices(misf_order, level_sizes, level_index)
    current[inserted, ] <- coords_post_insertion[inserted, , drop = FALSE]

    level_results[[level_index]] <- mek_refine_level_reference(
      n = n,
      edges = edges,
      edge_lengths = edge_lengths,
      coords_post_insertion = current,
      misf_order = misf_order,
      level_sizes = level_sizes,
      level_index = level_index,
      candidate_metric_pairs = candidate_metric_pairs,
      rho = rho_schedule[[level_index]],
      lambda = lambda_schedule[[level_index]],
      anchor_weight = anchor_weight_schedule[[level_index]],
      exact_repulsion_below = exact_repulsion_below,
      repulsion_sample_count = repulsion_sample_count,
      repulsion_seed = repulsion_seed,
      max_iter = max_iter,
      initial_step = initial_step,
      edge_scale = edge_scale,
      metric_scale = metric_scale,
      scale_mode = scale_mode,
      repulsion_weight_mode = repulsion_weight_mode
    )
    current <- level_results[[level_index]]$coords
    coords_by_level[[level_index]] <- current
  }

  level_trace <- do.call(rbind, lapply(level_results, `[[`, "level_trace"))
  rownames(level_trace) <- NULL

  list(
    coords = current,
    level_trace = level_trace,
    level_results = level_results,
    armijo_traces = lapply(level_results, `[[`, "armijo_trace"),
    coords_by_level = coords_by_level,
    insertion_coords = coords_post_insertion,
    rho_schedule = rho_schedule,
    lambda_schedule = lambda_schedule,
    anchor_weight_schedule = anchor_weight_schedule,
    repulsion_weight_mode = repulsion_weight_mode
  )
}
