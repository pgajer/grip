#!/usr/bin/env Rscript

sys.source(file.path("tools", "benchmarks", "gkk_lgkk_paper", "benchmark-sierpinski-baseline.R"), envir = environment())
sys.source(file.path("tools", "reports", "gkk_lgkk_paper", "generate-sierpinski-diagnostics.R"), envir = environment())

search_space <- list(
  placement = c("barycenter"),
  rounds = c(64L, 96L, 128L, 160L, 192L, 224L),
  final_rounds = c(128L, 160L, 192L, 224L, 256L, 288L, 320L, 384L),
  num_init = c(12L, 20L, 28L, 36L, 49L, 64L),
  num_nbrs = c(8L, 12L, 16L, 20L, 24L, 28L, 32L),
  r = c(0.01, 0.03, 0.05, 0.07, 0.10, 0.15),
  s = c(1.5, 3.0, 4.5, 6.0, 7.5, 9.0),
  repulsion_factor = c(0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0)
)

score_weights <- c(
  procrustes_rmse = 0.20,
  edge_length_cv = 0.10,
  sampled_stress = 0.10,
  sampled_nonedge_sep_ratio = 0.10,
  torus_cross_section_circle_rmse = 0.20,
  torus_major_cycle_angle_rmse = 0.15,
  torus_minor_cycle_azimuth_sd = 0.15
)

candidate_fields <- c(
  "placement", "rounds", "final_rounds", "num_init",
  "num_nbrs", "r", "s", "repulsion_factor"
)

parse_named_args <- function(args) {
  out <- list()
  for (arg in args) {
    parts <- strsplit(arg, "=", fixed = TRUE)[[1L]]
    if (length(parts) != 2L) {
      stop("Arguments must use key=value format")
    }
    out[[parts[[1L]]]] <- parts[[2L]]
  }
  out
}

parse_int_scalar <- function(x, name) {
  val <- suppressWarnings(as.integer(x))
  if (length(val) != 1L || is.na(val)) {
    stop(sprintf("%s must be a single integer", name))
  }
  val
}

parse_int_vector <- function(x, name) {
  parts <- trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
  val <- suppressWarnings(as.integer(parts))
  if (length(val) == 0L || any(is.na(val))) {
    stop(sprintf("%s must be a comma-separated integer list", name))
  }
  val
}

parse_num_vector <- function(x, name) {
  parts <- trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
  val <- suppressWarnings(as.double(parts))
  if (length(val) == 0L || any(is.na(val) | !is.finite(val))) {
    stop(sprintf("%s must be a comma-separated numeric list", name))
  }
  val
}

parse_char_vector <- function(x, name) {
  parts <- trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0L) {
    stop(sprintf("%s must be a comma-separated list", name))
  }
  parts
}

parse_size_specs <- function(x) {
  parts <- trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0L) {
    stop("sizes must be a comma-separated list such as 8x8,12x12")
  }

  out <- vector("list", length(parts))
  for (i in seq_along(parts)) {
    part <- parts[[i]]
    if (grepl("^[0-9]+$", part)) {
      h <- as.integer(part)
      w <- h
    } else if (grepl("^[0-9]+x[0-9]+$", part, ignore.case = TRUE)) {
      dims <- strsplit(tolower(part), "x", fixed = TRUE)[[1L]]
      h <- as.integer(dims[[1L]])
      w <- as.integer(dims[[2L]])
    } else {
      stop("sizes must use n or hxw format, for example 8 or 8x12")
    }
    if (is.na(h) || is.na(w) || h <= 2L || w <= 2L) {
      stop("torus sizes must be integers >= 3")
    }
    out[[i]] <- list(h = h, w = w, label = sprintf("%dx%d", h, w))
  }
  out
}

validate_run_tag <- function(x) {
  if (!grepl("^[A-Za-z0-9._-]+$", x)) {
    stop("tag must match ^[A-Za-z0-9._-]+$")
  }
  x
}

apply_search_overrides <- function(space, args) {
  for (field in names(space)) {
    if (is.null(args[[field]])) {
      next
    }
    value <- switch(
      field,
      placement = parse_char_vector(args[[field]], field),
      rounds = parse_int_vector(args[[field]], field),
      final_rounds = parse_int_vector(args[[field]], field),
      num_init = parse_int_vector(args[[field]], field),
      num_nbrs = parse_int_vector(args[[field]], field),
      r = parse_num_vector(args[[field]], field),
      s = parse_num_vector(args[[field]], field),
      repulsion_factor = parse_num_vector(args[[field]], field)
    )
    if (field == "placement") {
      bad <- setdiff(value, c("barycenter"))
      if (length(bad) > 0L) {
        stop("placement values for torus tuning must be 'barycenter'")
      }
    }
    if (field %in% c("rounds", "final_rounds", "num_init", "num_nbrs") && any(value <= 0L)) {
      stop(sprintf("%s values must be positive integers", field))
    }
    if (field == "r" && any(value < 0 | value > 1)) {
      stop("r values must be in [0, 1]")
    }
    if (field == "s" && any(value < 0)) {
      stop("s values must be >= 0")
    }
    if (field == "repulsion_factor" && any(value < 0)) {
      stop("repulsion_factor values must be >= 0")
    }
    space[[field]] <- value
  }
  space
}

default_torus_sizes <- function() {
  parse_size_specs("8x8,12x12,16x16")
}

tetrahedron_baseline_profile <- function() {
  list(
    placement = "barycenter",
    rounds = 128L,
    final_rounds = 192L,
    num_init = 12L,
    num_nbrs = 16L,
    r = 0.07,
    s = 9.0,
    repulsion_factor = 1.5
  )
}

resolve_baseline_profile <- function(args) {
  profile <- tetrahedron_baseline_profile()
  overrides <- list(
    placement = if (!is.null(args$baseline_placement)) parse_char_vector(args$baseline_placement, "baseline_placement")[[1L]] else NULL,
    rounds = if (!is.null(args$baseline_rounds)) parse_int_scalar(args$baseline_rounds, "baseline_rounds") else NULL,
    final_rounds = if (!is.null(args$baseline_final_rounds)) parse_int_scalar(args$baseline_final_rounds, "baseline_final_rounds") else NULL,
    num_init = if (!is.null(args$baseline_num_init)) parse_int_scalar(args$baseline_num_init, "baseline_num_init") else NULL,
    num_nbrs = if (!is.null(args$baseline_num_nbrs)) parse_int_scalar(args$baseline_num_nbrs, "baseline_num_nbrs") else NULL,
    r = if (!is.null(args$baseline_r)) as.double(parse_num_vector(args$baseline_r, "baseline_r")[[1L]]) else NULL,
    s = if (!is.null(args$baseline_s)) as.double(parse_num_vector(args$baseline_s, "baseline_s")[[1L]]) else NULL,
    repulsion_factor = if (!is.null(args$baseline_repulsion_factor)) as.double(parse_num_vector(args$baseline_repulsion_factor, "baseline_repulsion_factor")[[1L]]) else NULL
  )
  for (field in names(overrides)) {
    value <- overrides[[field]]
    if (!is.null(value)) {
      profile[[field]] <- value
    }
  }

  if (!identical(profile$placement, "barycenter")) {
    stop("baseline_placement must be 'barycenter'")
  }
  if (profile$rounds <= 0L || profile$final_rounds <= 0L || profile$num_init <= 0L || profile$num_nbrs <= 0L) {
    stop("baseline integer parameters must be positive")
  }
  if (profile$r < 0 || profile$r > 1) {
    stop("baseline_r must be in [0, 1]")
  }
  if (profile$s < 0) {
    stop("baseline_s must be >= 0")
  }
  if (profile$repulsion_factor < 0) {
    stop("baseline_repulsion_factor must be >= 0")
  }
  profile
}

candidate_key <- function(candidate) {
  paste(
    candidate$placement,
    candidate$rounds,
    candidate$final_rounds,
    candidate$num_init,
    candidate$num_nbrs,
    sprintf("%.4f", candidate$r),
    sprintf("%.4f", candidate$s),
    sprintf("%.4f", candidate$repulsion_factor),
    sep = "|"
  )
}

candidate_row <- function(candidate) {
  out <- as.data.frame(candidate[c("candidate_id", "candidate_source", candidate_fields)],
                       stringsAsFactors = FALSE)
  out$rounds <- as.integer(out$rounds)
  out$final_rounds <- as.integer(out$final_rounds)
  out$num_init <- as.integer(out$num_init)
  out$num_nbrs <- as.integer(out$num_nbrs)
  out$r <- as.double(out$r)
  out$s <- as.double(out$s)
  out$repulsion_factor <- as.double(out$repulsion_factor)
  out
}

format_candidate <- function(row) {
  sprintf(
    "placement=%s, rounds=%d, final_rounds=%d, num_init=%d, num_nbrs=%d, r=%.2f, s=%.1f, repulsion_factor=%.2f",
    row$placement,
    row$rounds,
    row$final_rounds,
    row$num_init,
    row$num_nbrs,
    row$r,
    row$s,
    row$repulsion_factor
  )
}

sample_candidate <- function(space) {
  rounds <- sample(space$rounds, 1L)
  final_choices <- space$final_rounds[space$final_rounds >= rounds]
  if (length(final_choices) == 0L) {
    final_choices <- space$final_rounds
  }
  list(
    placement = "barycenter",
    rounds = as.integer(rounds),
    final_rounds = as.integer(sample(final_choices, 1L)),
    num_init = as.integer(sample(space$num_init, 1L)),
    num_nbrs = as.integer(sample(space$num_nbrs, 1L)),
    r = as.double(sample(space$r, 1L)),
    s = as.double(sample(space$s, 1L)),
    repulsion_factor = as.double(sample(space$repulsion_factor, 1L))
  )
}

generate_candidates <- function(space, n_random, search_seed, include_default_reference = TRUE,
                                include_carpet_reference = TRUE,
                                baseline_profile = tetrahedron_baseline_profile(),
                                baseline_id = "torus_tetrahedron_baseline") {
  baseline_candidate <- c(
    list(candidate_id = baseline_id, candidate_source = "baseline"),
    baseline_profile[candidate_fields]
  )
  candidates <- list(baseline_candidate)
  seen <- candidate_key(baseline_candidate)

  if (isTRUE(include_carpet_reference)) {
    carpet_candidate <- c(
      list(candidate_id = "torus_carpet_reference", candidate_source = "reference"),
      carpet_preset_profile()[candidate_fields]
    )
    candidates[[length(candidates) + 1L]] <- carpet_candidate
    seen <- c(seen, candidate_key(carpet_candidate))
  }

  if (isTRUE(include_default_reference)) {
    default_candidate <- c(
      list(candidate_id = "torus_default_reference", candidate_source = "reference"),
      default_profile()[candidate_fields]
    )
    candidates[[length(candidates) + 1L]] <- default_candidate
    seen <- c(seen, candidate_key(default_candidate))
  }

  set.seed(search_seed)
  attempts <- 0L
  max_attempts <- max(1000L, n_random * 250L)
  while (sum(vapply(candidates, function(x) identical(x$candidate_source, "random"), logical(1L))) < n_random &&
         attempts < max_attempts) {
    attempts <- attempts + 1L
    sampled <- sample_candidate(space)
    key <- candidate_key(sampled)
    if (key %in% seen) {
      next
    }
    idx <- sum(vapply(candidates, function(x) identical(x$candidate_source, "random"), logical(1L))) + 1L
    candidates[[length(candidates) + 1L]] <- c(
      list(candidate_id = sprintf("torus_rand_%03d", idx), candidate_source = "random"),
      sampled
    )
    seen <- c(seen, key)
  }

  candidates
}

coords_to_grid_array <- function(coords, h, w) {
  d <- ncol(coords)
  arr <- array(0, dim = c(h, w, d))
  idx <- 1L
  for (i in seq_len(h)) {
    for (j in seq_len(w)) {
      arr[i, j, ] <- coords[idx, ]
      idx <- idx + 1L
    }
  }
  arr
}

grid_array_to_coords <- function(arr) {
  d <- dim(arr)
  coords <- matrix(0, nrow = d[[1L]] * d[[2L]], ncol = d[[3L]])
  idx <- 1L
  for (i in seq_len(d[[1L]])) {
    for (j in seq_len(d[[2L]])) {
      coords[idx, ] <- arr[i, j, ]
      idx <- idx + 1L
    }
  }
  coords
}

shift_grid_array <- function(arr, shift, axis) {
  if (!shift) {
    return(arr)
  }
  if (axis == 1L) {
    idx <- ((seq_len(dim(arr)[1L]) - 1L + shift) %% dim(arr)[1L]) + 1L
    return(arr[idx, , , drop = FALSE])
  }
  idx <- ((seq_len(dim(arr)[2L]) - 1L + shift) %% dim(arr)[2L]) + 1L
  arr[, idx, , drop = FALSE]
}

align_to_target_orthogonal <- function(source, target, allow_reflection = TRUE) {
  src <- normalize_coords(source)
  dst <- normalize_coords(target)
  cross <- t(src) %*% dst
  sv <- svd(cross)
  rot <- sv$u %*% t(sv$v)
  if (!allow_reflection && det(rot) < 0) {
    fix <- diag(ncol(rot))
    fix[ncol(fix), ncol(fix)] <- -1
    rot <- sv$u %*% fix %*% t(sv$v)
  }
  aligned <- src %*% rot
  list(
    aligned = aligned,
    target = dst,
    rmse = sqrt(mean(rowSums((aligned - dst)^2)))
  )
}

build_torus_graph <- function(h, w, major_radius = 2.6, minor_radius = 1.0) {
  h <- as.integer(h)
  w <- as.integer(w)
  edges <- edges.torus(h, w)
  coords <- matrix(0, nrow = h * w, ncol = 3L)
  idx <- 1L
  for (i in 0:(h - 1L)) {
    theta <- 2 * pi * i / h
    for (j in 0:(w - 1L)) {
      phi <- 2 * pi * j / w
      ring <- major_radius + minor_radius * cos(theta)
      coords[idx, ] <- c(
        ring * cos(phi),
        ring * sin(phi),
        minor_radius * sin(theta)
      )
      idx <- idx + 1L
    }
  }
  list(edges = edges, coords = coords)
}

build_swapped_torus_coords <- function(h, w, major_radius = 2.6, minor_radius = 1.0) {
  swapped <- build_torus_graph(w, h, major_radius = major_radius, minor_radius = minor_radius)$coords
  swapped_arr <- coords_to_grid_array(swapped, w, h)
  out_arr <- array(0, dim = c(h, w, 3L))
  for (i in seq_len(h)) {
    for (j in seq_len(w)) {
      out_arr[i, j, ] <- swapped_arr[j, i, ]
    }
  }
  grid_array_to_coords(out_arr)
}

best_torus_alignment <- function(source_coords, spec) {
  variants <- list(
    list(name = "canonical", arr = spec$canonical_arr, axis = 1L, count = spec$h),
    list(name = "swapped_cycles", arr = spec$swapped_arr, axis = 2L, count = spec$w)
  )

  best <- NULL
  for (variant in variants) {
    for (shift in 0:(variant$count - 1L)) {
      target <- grid_array_to_coords(shift_grid_array(variant$arr, shift, variant$axis))
      fit <- align_to_target_orthogonal(source_coords, target, allow_reflection = TRUE)
      if (is.null(best) || fit$rmse < best$rmse) {
        best <- c(
          fit,
          list(
            variant = variant$name,
            shift = shift
          )
        )
      }
    }
  }
  best
}

build_graph_specs_for_sizes <- function(sizes) {
  lapply(sizes, function(size) {
    built <- build_torus_graph(size$h, size$w)
    package_edges <- edges.torus(size$h, size$w)
    if (!identical(unname(built$edges), unname(package_edges))) {
      stop(sprintf("Canonical torus builder does not match edges.torus(%d, %d)", size$h, size$w))
    }
    canonical <- built$coords
    swapped <- build_swapped_torus_coords(size$h, size$w)
    list(
      family = "torus",
      graph_label = size$label,
      h = size$h,
      w = size$w,
      edges = package_edges,
      canonical = canonical,
      canonical_arr = coords_to_grid_array(canonical, size$h, size$w),
      swapped_arr = coords_to_grid_array(swapped, size$h, size$w)
    )
  })
}

circular_diff <- function(a, b) {
  atan2(sin(a - b), cos(a - b))
}

circular_sd <- function(theta) {
  z <- mean(exp(1i * theta))
  rbar <- Mod(z)
  if (!is.finite(rbar) || rbar <= 1e-12) {
    return(pi / sqrt(3))
  }
  sqrt(max(-2 * log(rbar), 0))
}

fit_circle_rmse <- function(x, y) {
  if (length(x) < 3L || length(y) < 3L) {
    return(NA_real_)
  }
  A <- cbind(2 * x, 2 * y, 1)
  b <- x^2 + y^2
  coef <- tryCatch(qr.solve(A, b), error = function(e) NULL)
  if (is.null(coef)) {
    return(NA_real_)
  }
  cx <- coef[[1L]]
  cy <- coef[[2L]]
  c0 <- coef[[3L]]
  radius_sq <- c0 + cx^2 + cy^2
  if (!is.finite(radius_sq) || radius_sq <= 0) {
    return(NA_real_)
  }
  radius <- sqrt(radius_sq)
  dist <- sqrt((x - cx)^2 + (y - cy)^2)
  sqrt(mean((dist - radius)^2)) / radius
}

torus_metric_array <- function(aligned_coords, spec, variant) {
  arr <- coords_to_grid_array(aligned_coords, spec$h, spec$w)
  if (identical(variant, "swapped_cycles")) {
    arr <- aperm(arr, c(2L, 1L, 3L))
  }
  arr
}

torus_shape_metrics <- function(aligned_coords, spec, variant) {
  arr <- torus_metric_array(aligned_coords, spec, variant)
  dims <- dim(arr)
  n_minor <- dims[[1L]]
  n_major <- dims[[2L]]

  rho <- sqrt(aligned_coords[, 1L]^2 + aligned_coords[, 2L]^2)
  z <- aligned_coords[, 3L]
  cross_section_rmse <- fit_circle_rmse(rho, z)

  ideal_step <- 2 * pi / n_major
  major_cycle_rmse <- mean(vapply(seq_len(n_minor), function(i) {
    phi <- atan2(arr[i, , 2L], arr[i, , 1L])
    steps <- c(
      abs(vapply(seq_len(n_major - 1L), function(j) circular_diff(phi[j + 1L], phi[j]), numeric(1L))),
      abs(circular_diff(phi[1L], phi[n_major]))
    )
    sqrt(mean((steps - ideal_step)^2)) / ideal_step
  }, numeric(1L)))

  minor_cycle_azimuth_sd <- mean(vapply(seq_len(n_major), function(j) {
    phi <- atan2(arr[, j, 2L], arr[, j, 1L])
    circular_sd(phi)
  }, numeric(1L)))

  list(
    torus_cross_section_circle_rmse = cross_section_rmse,
    torus_major_cycle_angle_rmse = major_cycle_rmse,
    torus_minor_cycle_azimuth_sd = minor_cycle_azimuth_sd
  )
}

run_one_layout_torus <- function(spec, cfg, seed, stress_sample = 4000L, sep_sample = 8000L) {
  n <- max(spec$edges)
  adj <- make_adj_list(spec$edges, n)
  started <- proc.time()[["elapsed"]]
  coords <- legacy.grip(
    edges = spec$edges,
    n = n,
    dim = 3,
    placement = cfg$placement,
    rounds = cfg$rounds,
    final_rounds = cfg$final_rounds,
    num_init = cfg$num_init,
    num_nbrs = cfg$num_nbrs,
    r = cfg$r,
    s = cfg$s,
    repulsion_factor = cfg$repulsion_factor,
    seed = seed
  )
  elapsed <- proc.time()[["elapsed"]] - started
  aligned <- best_torus_alignment(coords, spec)
  edge_stats <- edge_length_stats(coords, spec$edges)
  shape <- torus_shape_metrics(aligned$aligned, spec, aligned$variant)

  data.frame(
    family = spec$family,
    graph_label = spec$graph_label,
    h = spec$h,
    w = spec$w,
    seed = seed,
    vertices = n,
    edges = nrow(spec$edges),
    placement = cfg$placement,
    rounds = cfg$rounds,
    final_rounds = cfg$final_rounds,
    num_init = cfg$num_init,
    num_nbrs = cfg$num_nbrs,
    r = cfg$r,
    s = cfg$s,
    repulsion_factor = cfg$repulsion_factor,
    procrustes_rmse = aligned$rmse,
    edge_length_cv = edge_stats$cv,
    median_edge_length = edge_stats$median,
    sampled_stress = sampled_stress(coords, adj, sample_size = stress_sample, rng_seed = 1000L + seed),
    sampled_nonedge_sep_ratio = sampled_nonedge_separation_ratio(
      coords, spec$edges, sample_size = sep_sample, rng_seed = 2000L + seed
    ),
    torus_cross_section_circle_rmse = shape$torus_cross_section_circle_rmse,
    torus_major_cycle_angle_rmse = shape$torus_major_cycle_angle_rmse,
    torus_minor_cycle_azimuth_sd = shape$torus_minor_cycle_azimuth_sd,
    elapsed_sec = elapsed,
    align_variant = aligned$variant,
    align_shift = aligned$shift,
    stringsAsFactors = FALSE
  )
}

run_one_layout_safe <- function(spec, candidate, seed) {
  started <- proc.time()[["elapsed"]]
  metrics <- tryCatch(
    run_one_layout_torus(spec, candidate, seed),
    error = function(e) {
      n <- max(spec$edges)
      data.frame(
        family = spec$family,
        graph_label = spec$graph_label,
        h = spec$h,
        w = spec$w,
        seed = seed,
        vertices = n,
        edges = nrow(spec$edges),
        placement = candidate$placement,
        rounds = candidate$rounds,
        final_rounds = candidate$final_rounds,
        num_init = candidate$num_init,
        num_nbrs = candidate$num_nbrs,
        r = candidate$r,
        s = candidate$s,
        repulsion_factor = candidate$repulsion_factor,
        procrustes_rmse = NA_real_,
        edge_length_cv = NA_real_,
        median_edge_length = NA_real_,
        sampled_stress = NA_real_,
        sampled_nonedge_sep_ratio = NA_real_,
        torus_cross_section_circle_rmse = NA_real_,
        torus_major_cycle_angle_rmse = NA_real_,
        torus_minor_cycle_azimuth_sd = NA_real_,
        elapsed_sec = NA_real_,
        align_variant = "",
        align_shift = NA_integer_,
        stringsAsFactors = FALSE
      )
    }
  )
  if (!"elapsed_sec" %in% names(metrics) || !is.finite(metrics$elapsed_sec[[1L]])) {
    metrics$elapsed_sec <- proc.time()[["elapsed"]] - started
  }
  metrics$candidate_id <- candidate$candidate_id
  metrics$candidate_source <- candidate$candidate_source
  metrics$status <- if (is.finite(metrics$procrustes_rmse[[1L]])) "ok" else "error"
  metrics$error_message <- if (identical(metrics$status[[1L]], "ok")) "" else "layout failed"
  metrics[, c(
    "family", "graph_label", "h", "w", "candidate_id", "candidate_source", "seed",
    "vertices", "edges", "placement", "rounds", "final_rounds", "num_init",
    "num_nbrs", "r", "s", "repulsion_factor", "status", "error_message",
    "procrustes_rmse", "edge_length_cv", "median_edge_length", "sampled_stress",
    "sampled_nonedge_sep_ratio", "torus_cross_section_circle_rmse",
    "torus_major_cycle_angle_rmse", "torus_minor_cycle_azimuth_sd",
    "elapsed_sec", "align_variant", "align_shift"
  )]
}

rank01 <- function(x, higher_better = FALSE) {
  n <- length(x)
  out <- rep(1, n)
  ok <- is.finite(x)
  if (!any(ok)) {
    return(out)
  }
  vals <- if (higher_better) -x[ok] else x[ok]
  if (length(vals) == 1L || all(abs(vals - vals[[1L]]) <= sqrt(.Machine$double.eps))) {
    out[ok] <- 0
    return(out)
  }
  ranks <- rank(vals, ties.method = "average")
  out[ok] <- (ranks - 1) / (length(ranks) - 1)
  out
}

summarize_candidate_graphs <- function(raw_metrics) {
  groups <- split(raw_metrics, paste(raw_metrics$graph_label, raw_metrics$candidate_id, sep = "|"))
  do.call(rbind, lapply(groups, function(df) {
    ok <- df$status == "ok"
    good <- df[ok, , drop = FALSE]
    best_seed <- if (nrow(good) > 0L) good$seed[[which.min(good$procrustes_rmse)]] else NA_integer_
    best_variant <- if (nrow(good) > 0L) good$align_variant[[which.min(good$procrustes_rmse)]] else ""
    best_shift <- if (nrow(good) > 0L) good$align_shift[[which.min(good$procrustes_rmse)]] else NA_integer_
    data.frame(
      family = df$family[[1L]],
      graph_label = df$graph_label[[1L]],
      h = df$h[[1L]],
      w = df$w[[1L]],
      candidate_id = df$candidate_id[[1L]],
      candidate_source = df$candidate_source[[1L]],
      placement = df$placement[[1L]],
      rounds = df$rounds[[1L]],
      final_rounds = df$final_rounds[[1L]],
      num_init = df$num_init[[1L]],
      num_nbrs = df$num_nbrs[[1L]],
      r = df$r[[1L]],
      s = df$s[[1L]],
      repulsion_factor = df$repulsion_factor[[1L]],
      vertices = df$vertices[[1L]],
      edges = df$edges[[1L]],
      n_runs = nrow(df),
      n_ok = sum(ok),
      n_fail = sum(!ok),
      procrustes_rmse_mean = if (nrow(good) > 0L) mean(good$procrustes_rmse) else NA_real_,
      procrustes_rmse_sd = if (nrow(good) > 1L) stats::sd(good$procrustes_rmse) else 0,
      edge_length_cv_mean = if (nrow(good) > 0L) mean(good$edge_length_cv) else NA_real_,
      edge_length_cv_sd = if (nrow(good) > 1L) stats::sd(good$edge_length_cv) else 0,
      sampled_stress_mean = if (nrow(good) > 0L) mean(good$sampled_stress) else NA_real_,
      sampled_stress_sd = if (nrow(good) > 1L) stats::sd(good$sampled_stress) else 0,
      sampled_nonedge_sep_ratio_mean = if (nrow(good) > 0L) mean(good$sampled_nonedge_sep_ratio) else NA_real_,
      sampled_nonedge_sep_ratio_sd = if (nrow(good) > 1L) stats::sd(good$sampled_nonedge_sep_ratio) else 0,
      torus_cross_section_circle_rmse_mean = if (nrow(good) > 0L) mean(good$torus_cross_section_circle_rmse) else NA_real_,
      torus_cross_section_circle_rmse_sd = if (nrow(good) > 1L) stats::sd(good$torus_cross_section_circle_rmse) else 0,
      torus_major_cycle_angle_rmse_mean = if (nrow(good) > 0L) mean(good$torus_major_cycle_angle_rmse) else NA_real_,
      torus_major_cycle_angle_rmse_sd = if (nrow(good) > 1L) stats::sd(good$torus_major_cycle_angle_rmse) else 0,
      torus_minor_cycle_azimuth_sd_mean = if (nrow(good) > 0L) mean(good$torus_minor_cycle_azimuth_sd) else NA_real_,
      torus_minor_cycle_azimuth_sd_sd = if (nrow(good) > 1L) stats::sd(good$torus_minor_cycle_azimuth_sd) else 0,
      elapsed_sec_mean = if (nrow(good) > 0L) mean(good$elapsed_sec) else NA_real_,
      best_seed = best_seed,
      best_align_variant = best_variant,
      best_align_shift = best_shift,
      stringsAsFactors = FALSE
    )
  }))
}

score_candidate_graphs <- function(graph_summary) {
  split_groups <- split(graph_summary, graph_summary$graph_label)
  scored <- lapply(split_groups, function(df) {
    df$rank_rmse <- rank01(df$procrustes_rmse_mean, higher_better = FALSE)
    df$rank_edge_cv <- rank01(df$edge_length_cv_mean, higher_better = FALSE)
    df$rank_stress <- rank01(df$sampled_stress_mean, higher_better = FALSE)
    df$rank_sep <- rank01(df$sampled_nonedge_sep_ratio_mean, higher_better = TRUE)
    df$rank_cross_section <- rank01(df$torus_cross_section_circle_rmse_mean, higher_better = FALSE)
    df$rank_major_cycle <- rank01(df$torus_major_cycle_angle_rmse_mean, higher_better = FALSE)
    df$rank_minor_cycle <- rank01(df$torus_minor_cycle_azimuth_sd_mean, higher_better = FALSE)
    df$score_graph <-
      score_weights[["procrustes_rmse"]] * df$rank_rmse +
      score_weights[["edge_length_cv"]] * df$rank_edge_cv +
      score_weights[["sampled_stress"]] * df$rank_stress +
      score_weights[["sampled_nonedge_sep_ratio"]] * df$rank_sep +
      score_weights[["torus_cross_section_circle_rmse"]] * df$rank_cross_section +
      score_weights[["torus_major_cycle_angle_rmse"]] * df$rank_major_cycle +
      score_weights[["torus_minor_cycle_azimuth_sd"]] * df$rank_minor_cycle
    df$score_graph[!is.finite(df$procrustes_rmse_mean)] <- Inf
    df[order(df$score_graph, df$procrustes_rmse_mean), , drop = FALSE]
  })
  do.call(rbind, scored)
}

summarize_family_rankings <- function(scored_graphs) {
  groups <- split(scored_graphs, scored_graphs$candidate_id)
  out <- do.call(rbind, lapply(groups, function(df) {
    data.frame(
      family = df$family[[1L]],
      candidate_id = df$candidate_id[[1L]],
      candidate_source = df$candidate_source[[1L]],
      placement = df$placement[[1L]],
      rounds = df$rounds[[1L]],
      final_rounds = df$final_rounds[[1L]],
      num_init = df$num_init[[1L]],
      num_nbrs = df$num_nbrs[[1L]],
      r = df$r[[1L]],
      s = df$s[[1L]],
      repulsion_factor = df$repulsion_factor[[1L]],
      n_graphs = nrow(df),
      n_fail_graphs = sum(!is.finite(df$score_graph)),
      score_mean = if (all(!is.finite(df$score_graph))) Inf else mean(df$score_graph[is.finite(df$score_graph)]),
      score_sd = if (sum(is.finite(df$score_graph)) > 1L) stats::sd(df$score_graph[is.finite(df$score_graph)]) else 0,
      procrustes_rmse_mean = mean(df$procrustes_rmse_mean, na.rm = TRUE),
      edge_length_cv_mean = mean(df$edge_length_cv_mean, na.rm = TRUE),
      sampled_stress_mean = mean(df$sampled_stress_mean, na.rm = TRUE),
      sampled_nonedge_sep_ratio_mean = mean(df$sampled_nonedge_sep_ratio_mean, na.rm = TRUE),
      torus_cross_section_circle_rmse_mean = mean(df$torus_cross_section_circle_rmse_mean, na.rm = TRUE),
      torus_major_cycle_angle_rmse_mean = mean(df$torus_major_cycle_angle_rmse_mean, na.rm = TRUE),
      torus_minor_cycle_azimuth_sd_mean = mean(df$torus_minor_cycle_azimuth_sd_mean, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  out[order(out$score_mean, out$procrustes_rmse_mean), , drop = FALSE]
}

format_num <- function(x, digits = 4) {
  ifelse(is.finite(x), format(round(x, digits), nsmall = digits, trim = TRUE), "NA")
}

build_graph_comparison <- function(scored_graphs, comparison_candidate_id, baseline_candidate_id) {
  labels <- unique(scored_graphs$graph_label)
  labels <- labels[order(match(labels, scored_graphs$graph_label))]
  do.call(rbind, lapply(labels, function(label) {
    graph_df <- scored_graphs[scored_graphs$graph_label == label, , drop = FALSE]
    baseline_row <- graph_df[graph_df$candidate_id == baseline_candidate_id, , drop = FALSE]
    tuned_row <- graph_df[graph_df$candidate_id == comparison_candidate_id, , drop = FALSE]
    data.frame(
      family = "torus",
      graph_label = label,
      h = graph_df$h[[1L]],
      w = graph_df$w[[1L]],
      baseline_score = baseline_row$score_graph[[1L]],
      tuned_score = tuned_row$score_graph[[1L]],
      baseline_rmse = baseline_row$procrustes_rmse_mean[[1L]],
      tuned_rmse = tuned_row$procrustes_rmse_mean[[1L]],
      baseline_edge_cv = baseline_row$edge_length_cv_mean[[1L]],
      tuned_edge_cv = tuned_row$edge_length_cv_mean[[1L]],
      baseline_stress = baseline_row$sampled_stress_mean[[1L]],
      tuned_stress = tuned_row$sampled_stress_mean[[1L]],
      baseline_sep = baseline_row$sampled_nonedge_sep_ratio_mean[[1L]],
      tuned_sep = tuned_row$sampled_nonedge_sep_ratio_mean[[1L]],
      baseline_cross_section = baseline_row$torus_cross_section_circle_rmse_mean[[1L]],
      tuned_cross_section = tuned_row$torus_cross_section_circle_rmse_mean[[1L]],
      baseline_major_cycle = baseline_row$torus_major_cycle_angle_rmse_mean[[1L]],
      tuned_major_cycle = tuned_row$torus_major_cycle_angle_rmse_mean[[1L]],
      baseline_minor_cycle = baseline_row$torus_minor_cycle_azimuth_sd_mean[[1L]],
      tuned_minor_cycle = tuned_row$torus_minor_cycle_azimuth_sd_mean[[1L]],
      stringsAsFactors = FALSE
    )
  }))
}

write_tuning_summary <- function(path,
                                 config,
                                 family_rankings,
                                 scored_graphs,
                                 graph_comparison,
                                 pdf_paths) {
  top_rows <- head(family_rankings, config$top_n)
  best_overall <- top_rows[1L, , drop = FALSE]
  best_by_graph <- do.call(rbind, lapply(split(scored_graphs, scored_graphs$graph_label), function(df) {
    df[order(df$score_graph, df$procrustes_rmse_mean), , drop = FALSE][1L, , drop = FALSE]
  }))
  best_by_graph <- best_by_graph[order(best_by_graph$h, best_by_graph$w), , drop = FALSE]

  lines <- c(
    "# Torus Parameter Tuning",
    "",
    "Search setup:",
    sprintf("- graph sizes: `%s`", paste(vapply(config$sizes, `[[`, "", "label"), collapse = ", ")),
    sprintf("- random candidates: `%d`", config$n_random),
    sprintf("- seeds: `%s`", paste(config$seeds, collapse = ", ")),
    sprintf("- search seed: `%d`", config$search_seed),
    sprintf("- baseline candidate: `%s`", config$baseline_candidate_id),
    sprintf("- comparison candidate: `%s`", config$comparison_candidate_id),
    sprintf("- reference candidates: `%s`", paste(config$reference_candidate_ids, collapse = ", ")),
    "",
    "Score definition (lower is better):",
    sprintf(
      "- `score_graph = %.2f * rank(RMSE) + %.2f * rank(edge_length_cv) + %.2f * rank(sampled_stress) + %.2f * rank(-nonedge_sep) + %.2f * rank(cross_section_circle_rmse) + %.2f * rank(major_cycle_angle_rmse) + %.2f * rank(minor_cycle_azimuth_sd)`",
      score_weights[["procrustes_rmse"]],
      score_weights[["edge_length_cv"]],
      score_weights[["sampled_stress"]],
      score_weights[["sampled_nonedge_sep_ratio"]],
      score_weights[["torus_cross_section_circle_rmse"]],
      score_weights[["torus_major_cycle_angle_rmse"]],
      score_weights[["torus_minor_cycle_azimuth_sd"]]
    ),
    "",
    "Search space:",
    sprintf("- placement: `%s`", paste(config$search_space$placement, collapse = ", ")),
    sprintf("- rounds: `%s`", paste(config$search_space$rounds, collapse = ", ")),
    sprintf("- final_rounds: `%s`", paste(config$search_space$final_rounds, collapse = ", ")),
    sprintf("- num_init: `%s`", paste(config$search_space$num_init, collapse = ", ")),
    sprintf("- num_nbrs: `%s`", paste(config$search_space$num_nbrs, collapse = ", ")),
    sprintf("- r: `%s`", paste(config$search_space$r, collapse = ", ")),
    sprintf("- s: `%s`", paste(config$search_space$s, collapse = ", ")),
    sprintf("- repulsion_factor: `%s`", paste(config$search_space$repulsion_factor, collapse = ", ")),
    "",
    "Best overall candidate:",
    sprintf("- `%s`", best_overall$candidate_id[[1L]]),
    sprintf("- %s", format_candidate(best_overall[1L, ])),
    sprintf("- mean score across graphs: %s", format_num(best_overall$score_mean[[1L]], digits = 4)),
    "",
    "Top candidates:",
    "",
    "| Candidate | Source | Mean score | RMSE | Edge CV | Stress | Non-edge sep | Cross-sec | Major-cycle | Minor-cycle |",
    "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"
  )

  for (i in seq_len(nrow(top_rows))) {
    row <- top_rows[i, ]
    lines <- c(lines, sprintf(
      "| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
      row$candidate_id,
      row$candidate_source,
      format_num(row$score_mean, 4),
      format_num(row$procrustes_rmse_mean, 4),
      format_num(row$edge_length_cv_mean, 4),
      format_num(row$sampled_stress_mean, 4),
      format_num(row$sampled_nonedge_sep_ratio_mean, 4),
      format_num(row$torus_cross_section_circle_rmse_mean, 4),
      format_num(row$torus_major_cycle_angle_rmse_mean, 4),
      format_num(row$torus_minor_cycle_azimuth_sd_mean, 4)
    ))
  }

  lines <- c(lines, "", "Best candidate by graph size:", "", "| Size | Candidate | Score | Params |", "| --- | --- | ---: | --- |")
  for (i in seq_len(nrow(best_by_graph))) {
    row <- best_by_graph[i, ]
    lines <- c(lines, sprintf(
      "| %s | %s | %s | %s |",
      row$graph_label,
      row$candidate_id,
      format_num(row$score_graph, 4),
      format_candidate(row)
    ))
  }

  lines <- c(lines, "", "Baseline vs comparison candidate:", "", "| Size | Baseline score | Tuned score | Baseline RMSE | Tuned RMSE | Baseline edge CV | Tuned edge CV | Baseline stress | Tuned stress | Baseline sep | Tuned sep | Baseline cross-sec | Tuned cross-sec | Baseline major | Tuned major | Baseline minor | Tuned minor |", "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
  for (i in seq_len(nrow(graph_comparison))) {
    row <- graph_comparison[i, ]
    lines <- c(lines, sprintf(
      "| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
      row$graph_label,
      format_num(row$baseline_score, 4),
      format_num(row$tuned_score, 4),
      format_num(row$baseline_rmse, 4),
      format_num(row$tuned_rmse, 4),
      format_num(row$baseline_edge_cv, 4),
      format_num(row$tuned_edge_cv, 4),
      format_num(row$baseline_stress, 4),
      format_num(row$tuned_stress, 4),
      format_num(row$baseline_sep, 4),
      format_num(row$tuned_sep, 4),
      format_num(row$baseline_cross_section, 4),
      format_num(row$tuned_cross_section, 4),
      format_num(row$baseline_major_cycle, 4),
      format_num(row$tuned_major_cycle, 4),
      format_num(row$baseline_minor_cycle, 4),
      format_num(row$tuned_minor_cycle, 4)
    ))
  }

  lines <- c(lines, "", "Comparison PDFs:", "")
  lines <- c(lines, sprintf("- `%s`", pdf_paths))
  writeLines(lines, con = path)
}

if (sys.nframe() == 0L) {
  args <- parse_named_args(commandArgs(trailingOnly = TRUE))
  run_tag <- validate_run_tag(if (!is.null(args$tag)) args$tag else "torus-tuning")
  sizes <- if (!is.null(args$sizes)) parse_size_specs(args$sizes) else default_torus_sizes()
  n_random <- if (!is.null(args$n_random)) parse_int_scalar(args$n_random, "n_random") else 120L
  top_n <- if (!is.null(args$top_n)) parse_int_scalar(args$top_n, "top_n") else 8L
  search_seed <- if (!is.null(args$search_seed)) parse_int_scalar(args$search_seed, "search_seed") else 20260322L
  seeds <- if (!is.null(args$seeds)) parse_int_vector(args$seeds, "seeds") else 1:10
  include_default_reference <- if (!is.null(args$include_default_reference)) {
    tolower(args$include_default_reference) %in% c("1", "true", "yes")
  } else {
    TRUE
  }
  include_carpet_reference <- if (!is.null(args$include_carpet_reference)) {
    tolower(args$include_carpet_reference) %in% c("1", "true", "yes")
  } else {
    TRUE
  }

  if (n_random < 0L) stop("n_random must be >= 0")
  if (top_n <= 0L) stop("top_n must be >= 1")

  search_space <- apply_search_overrides(search_space, args)
  baseline_profile <- resolve_baseline_profile(args)
  baseline_candidate_id <- if (!is.null(args$baseline_id)) args$baseline_id else "torus_tetrahedron_baseline"

  tuning_root <- file.path("output", "gkk_lgkk_paper")
  tuning_pdf_dir <- file.path(tuning_root, "reports", run_tag, "torus")
  tuning_tmp_dir <- file.path(tuning_root, "tmp", run_tag)
  tuning_preview_dir <- file.path(tuning_tmp_dir, "pdf-previews")
  dir.create(tuning_pdf_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tuning_tmp_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tuning_preview_dir, recursive = TRUE, showWarnings = FALSE)

  candidates <- generate_candidates(
    search_space,
    n_random,
    search_seed,
    include_default_reference = include_default_reference,
    include_carpet_reference = include_carpet_reference,
    baseline_profile = baseline_profile,
    baseline_id = baseline_candidate_id
  )
  graphs <- build_graph_specs_for_sizes(sizes)

  raw_metrics <- do.call(
    rbind,
    lapply(graphs, function(spec) {
      do.call(
        rbind,
        lapply(candidates, function(candidate) {
          do.call(rbind, lapply(seeds, function(seed) run_one_layout_safe(spec, candidate, seed)))
        })
      )
    })
  )

  graph_summary <- summarize_candidate_graphs(raw_metrics)
  scored_graphs <- score_candidate_graphs(graph_summary)
  family_rankings <- summarize_family_rankings(scored_graphs)

  reference_candidate_ids <- c(
    if (include_carpet_reference) "torus_carpet_reference" else character(),
    if (include_default_reference) "torus_default_reference" else character()
  )
  best_candidate_id <- family_rankings$candidate_id[[1L]]
  comparison_candidate_id <- if (!is.null(args$compare_candidate_id)) args$compare_candidate_id else best_candidate_id
  if (!comparison_candidate_id %in% scored_graphs$candidate_id) {
    stop(sprintf("compare_candidate_id '%s' was not evaluated in this run", comparison_candidate_id))
  }
  graph_comparison <- build_graph_comparison(scored_graphs, comparison_candidate_id, baseline_candidate_id)

  comparison_pdfs <- character()
  for (spec in graphs) {
    graph_df <- scored_graphs[scored_graphs$graph_label == spec$graph_label, , drop = FALSE]
    baseline_row <- graph_df[graph_df$candidate_id == baseline_candidate_id, , drop = FALSE]
    tuned_row <- graph_df[graph_df$candidate_id == comparison_candidate_id, , drop = FALSE]

    baseline_coords <- legacy.grip(
      edges = spec$edges,
      n = max(spec$edges),
      dim = 3,
      placement = baseline_row$placement[[1L]],
      rounds = baseline_row$rounds[[1L]],
      final_rounds = baseline_row$final_rounds[[1L]],
      num_init = baseline_row$num_init[[1L]],
      num_nbrs = baseline_row$num_nbrs[[1L]],
      r = baseline_row$r[[1L]],
      s = baseline_row$s[[1L]],
      repulsion_factor = baseline_row$repulsion_factor[[1L]],
      seed = baseline_row$best_seed[[1L]]
    )
    tuned_coords <- legacy.grip(
      edges = spec$edges,
      n = max(spec$edges),
      dim = 3,
      placement = tuned_row$placement[[1L]],
      rounds = tuned_row$rounds[[1L]],
      final_rounds = tuned_row$final_rounds[[1L]],
      num_init = tuned_row$num_init[[1L]],
      num_nbrs = tuned_row$num_nbrs[[1L]],
      r = tuned_row$r[[1L]],
      s = tuned_row$s[[1L]],
      repulsion_factor = tuned_row$repulsion_factor[[1L]],
      seed = tuned_row$best_seed[[1L]]
    )

    baseline_aligned <- best_torus_alignment(baseline_coords, spec)$aligned
    tuned_aligned <- best_torus_alignment(tuned_coords, spec)$aligned
    pdf_path <- file.path(
      tuning_pdf_dir,
      sprintf("torus-%s-baseline-vs-best-overall-3d.pdf", spec$graph_label)
    )
    subtitle <- sprintf(
      "baseline=%s (seed=%d) | tuned=%s (seed=%d)",
      baseline_row$candidate_id[[1L]], baseline_row$best_seed[[1L]],
      tuned_row$candidate_id[[1L]], tuned_row$best_seed[[1L]]
    )
    write_3d_diagnostic_pdf(
      path = pdf_path,
      canonical_coords = spec$canonical,
      baseline_coords = baseline_aligned,
      tuned_coords = tuned_aligned,
      edges = spec$edges,
      title_text = sprintf("Torus %s", spec$graph_label),
      subtitle_text = subtitle,
      tuned_label = "tuned"
    )
    comparison_pdfs <- c(comparison_pdfs, pdf_path)
  }

  raw_csv_path <- file.path(tuning_tmp_dir, "torus-tuning-raw-metrics.csv")
  graph_csv_path <- file.path(tuning_tmp_dir, "torus-tuning-graph-summary.csv")
  family_csv_path <- file.path(tuning_tmp_dir, "torus-tuning-family-ranking.csv")
  comparison_csv_path <- file.path(tuning_tmp_dir, "torus-tuning-baseline-vs-best.csv")
  summary_md_path <- file.path(tuning_tmp_dir, "torus-tuning-summary.md")

  utils::write.csv(raw_metrics, raw_csv_path, row.names = FALSE)
  utils::write.csv(scored_graphs, graph_csv_path, row.names = FALSE)
  utils::write.csv(family_rankings, family_csv_path, row.names = FALSE)
  utils::write.csv(graph_comparison, comparison_csv_path, row.names = FALSE)
  write_tuning_summary(
    path = summary_md_path,
    config = list(
      sizes = sizes,
      n_random = n_random,
      top_n = top_n,
      seeds = seeds,
      search_seed = search_seed,
      baseline_candidate_id = baseline_candidate_id,
      comparison_candidate_id = comparison_candidate_id,
      reference_candidate_ids = reference_candidate_ids,
      search_space = search_space
    ),
    family_rankings = family_rankings,
    scored_graphs = scored_graphs,
    graph_comparison = graph_comparison,
    pdf_paths = comparison_pdfs
  )
  render_pdf_previews(comparison_pdfs, tuning_preview_dir)

  message(sprintf("Raw tuning metrics written to %s", raw_csv_path))
  message(sprintf("Graph summary written to %s", graph_csv_path))
  message(sprintf("Family ranking written to %s", family_csv_path))
  message(sprintf("Baseline-vs-best comparison written to %s", comparison_csv_path))
  message(sprintf("Markdown summary written to %s", summary_md_path))
  message(sprintf("Comparison PDFs written to %s", tuning_pdf_dir))
}
