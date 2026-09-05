grip.mds.has.smacof <- function() {
  requireNamespace("smacof", quietly = TRUE)
}

#' Metric stress MDS using SMACOF
#'
#' `metric.mds()` minimizes unweighted raw distance stress on a graph's
#' all-pairs shortest-path distances using `smacof::mds(type = "ratio")`.
#' It requires the optional \pkg{smacof} package. Before version 0.2.0.9000,
#' this name performed classical scaling; use [classical.mds()] to retain that
#' behavior. The `add` and `eig` arguments belong to `classical.mds()` only.
#'
#' @details The objective is
#' \deqn{S(Z) = \sum_{i<j}(\|z_i-z_j\|_2-\delta_{ij})^2.}
#' Edge weights define graph distances, not pair stiffnesses in this objective.
#' SMACOF normalizes targets internally. Returned coordinates are rescaled to
#' minimize raw stress against the original input distances. `scale_mode`
#' controls the optional diagnostic panel only, not the optimization objective.
#'
#' Target-normalized raw stress and scale-profiled Stress-1 select the same
#' shapes when global scale is free and pair weights agree. Their values at a
#' fixed coordinate scale need not agree. Both the literal and target-profiled
#' Stress-1 values are independently calculated in `metadata`; the backend's
#' own reported stress is retained separately in the start summaries.
#'
#' The first run uses `init`; subsequent runs use random configurations.
#' The smallest achieved raw stress selects the result. A backend result that
#' increases stress beyond numerical tolerance is rejected in favor of its
#' scaled start and marked accordingly. Failed starts are recorded, and all
#' failing starts cause an error. Iteration limits are not convergence or global
#' optimality certificates. Collinear or planar starts may remain in their
#' initial span; use multiple starts to investigate this sensitivity.
#'
#' Targets must be finite, symmetric, and nonnegative, with zero diagonal and
#' at least one positive distance. Missing and infinite distances are rejected.
#' Supplied starts may have coincident points but must not be wholly collapsed.
#' The implementation uses dense all-pairs matrices; edge-only refinement with
#' [edge.kk()] is preferable when that preparation is too large.
#'
#' @inheritParams classical.mds
#' @param prepared An all-pairs prepared graph object containing
#'   `distance_matrix`. Edge-only preparations are not supported.
#' @param diagnostics Attach the common GMDS diagnostic panel. With `FALSE`
#'   and raw graph inputs, prepare only the distance matrix, without path caches.
#' @param scale_mode Diagnostic scale policy: `"profiled"` fits a separate
#'   scalar for each diagnostic family, and `"identity"` uses scale one.
#'   To evaluate user-specified scales, call [score.gmds()] on the returned
#'   coordinates separately. This argument never changes the fitted coordinates.
#' @param init `"classical"` (default), `"random"`, or a finite numeric
#'   matrix with `n` rows and `dim` columns in input-distance units.
#' @param n_init Positive integer number of starts, including the first start.
#' @param max_iter Positive integer iteration limit for each SMACOF run.
#' @param eps Positive tolerance for the backend's change in normalized stress.
#' @param seed Integer random seed, or `NULL` to use the current RNG stream.
#'   With a non-NULL seed, random starts do not change the caller's RNG state.
#' @return A `"grip_gmds_layout"` object with method `"metric_mds"`.
#'   `metadata` records the objective, backend/version, achieved raw stress,
#'   target-normalized RMSE, both Stress-1 conventions, selected start,
#'   coordinate scale multiplier, and per-start losses and stopping information.
#' @seealso [classical.mds()], [edge.kk()], [smacof::mds()]
#' @md
#' @export
metric.mds <- function(prepared = NULL,
                       edges = NULL,
                       n = NULL,
                       adj_list = NULL,
                       weight_list = NULL,
                       edge_weights = NULL,
                       dim = 2L,
                       init = c("classical", "random"),
                       n_init = 1L,
                       max_iter = 1000L,
                       eps = 1e-8,
                       seed = 1L,
                       diagnostics = TRUE,
                       scale_mode = c("profiled", "identity"),
                       distance_floor = 1e-8,
                       edge_length_epsilon = 1e-8,
                       band_quantiles = c(1 / 3, 2 / 3)) {
  if (!grip.mds.has.smacof()) {
    stop("metric.mds() requires the optional 'smacof' package; install it, ",
         "or use classical.mds() for classical scaling", call. = FALSE)
  }
  if (utils::packageVersion("smacof") < "2.1-7") {
    stop("metric.mds() requires smacof >= 2.1-7", call. = FALSE)
  }
  scale_mode <- match.arg(scale_mode)
  for (name in c("dim", "n_init", "max_iter")) {
    value <- get(name)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 1 || value != floor(value) || value > .Machine$integer.max) {
      stop(name, " must be a positive integer", call. = FALSE)
    }
  }
  if (dim < 2L) stop("dim must be at least 2", call. = FALSE)
  grip.validate.scalar(eps, "eps", lower = 0, open.lower = TRUE)
  if (!is.logical(diagnostics) || length(diagnostics) != 1L || is.na(diagnostics)) {
    stop("diagnostics must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L ||
      !is.finite(seed) || seed != floor(seed) || abs(seed) > .Machine$integer.max)) {
    stop("seed must be an integer or NULL", call. = FALSE)
  }
  prepared <- if (is.null(prepared) && !diagnostics) {
    grip.metric.mds.distance.prepared(edges, n, adj_list, weight_list, edge_weights)
  } else {
    grip.gmds.require.prepared(prepared = prepared, edges = edges, n = n,
      adj_list = adj_list, weight_list = weight_list, edge_weights = edge_weights)
  }
  delta <- prepared$distance_matrix
  if (is.null(delta)) {
    stop("metric.mds() requires an all-pairs prepared object with distance_matrix; ",
         "prepare.edge.kk() objects are edge-only", call. = FALSE)
  }
  if (!is.matrix(delta) || !is.numeric(delta) ||
      !identical(dim(delta), c(as.integer(prepared$n), as.integer(prepared$n))) ||
      any(!is.finite(delta)) || any(delta < 0) || any(diag(delta) != 0) ||
      !isTRUE(all.equal(delta, t(delta), check.attributes = FALSE, tolerance = 1e-12))) {
    stop("distance_matrix must be finite, symmetric and nonnegative with zero diagonal",
         call. = FALSE)
  }
  if (dim >= prepared$n) stop("dim must be less than the number of vertices", call. = FALSE)
  target <- as.double(stats::as.dist(delta))
  # Scale before squaring to avoid overflow in the backend normalization.
  target.max <- max(target)
  if (target.max <= 0) stop("at least one target distance must be positive", call. = FALSE)
  target.rms <- target.max * sqrt(mean((target / target.max)^2))
  delta.normalized <- delta / target.rms
  target.normalized <- target / target.rms
  supplied <- is.matrix(init)
  if (supplied) {
    if (!is.numeric(init) || !identical(dim(init), c(as.integer(prepared$n), as.integer(dim))) ||
        any(!is.finite(init))) {
      stop("init must be a finite n by dim numeric matrix", call. = FALSE)
    }
    first <- sweep(init, 2L, colMeans(init), "-") / target.rms
    init.name <- "supplied"
  } else {
    init <- match.arg(init)
    init.name <- init
    first <- NULL
    if (identical(init, "classical")) {
      normalized.prepared <- prepared
      normalized.prepared$distance_matrix <- delta.normalized
      first <- classical.mds(prepared = normalized.prepared, dim = dim,
                             diagnostics = FALSE)$coords
    }
  }
  if (!is.null(seed) && (n_init > 1L || identical(init.name, "random"))) {
    had.seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old.seed <- if (had.seed) get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      if (had.seed) assign(".Random.seed", old.seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(seed)
  }
  rescale <- function(x) {
    x <- sweep(x, 2L, colMeans(x), "-")
    d <- as.double(stats::dist(x))
    denominator <- sum(d^2)
    if (!is.finite(denominator) || denominator <= 0) {
      stop("MDS configuration is collapsed or has nonfinite distances", call. = FALSE)
    }
    multiplier <- sum(d * target.normalized) / denominator
    list(coords = multiplier * x, multiplier = multiplier,
         loss = sum((multiplier * d - target.normalized)^2))
  }
  best <- NULL
  best.loss <- Inf
  records <- vector("list", n_init)
  for (run in seq_len(n_init)) {
    start <- if (run == 1L && !is.null(first)) first else
      matrix(stats::rnorm(prepared$n * dim), nrow = prepared$n, ncol = dim)
    start <- rescale(start)
    notices <- character()
    fit <- tryCatch(withCallingHandlers(
      smacof::mds(delta.normalized, ndim = dim, type = "ratio",
                  init = start$coords, itmax = max_iter, eps = eps,
                  principal = FALSE, verbose = FALSE),
      warning = function(w) {
        notices <<- c(notices, conditionMessage(w))
        invokeRestart("muffleWarning")
      }), error = function(e) e)
    failed <- inherits(fit, "error")
    result <- if (failed) fit else tryCatch(rescale(fit$conf), error = function(e) e)
    failed <- inherits(result, "error")
    rejected <- !failed && result$loss > start$loss + 1e-10 * max(1, start$loss)
    if (rejected) result <- start
    reason <- if (failed) "backend_error" else if (rejected) "rejected_increase" else
      if (fit$niter >= max_iter) "iteration_limit" else "stress_tolerance"
    records[[run]] <- data.frame(
      start = run, initialization = if (run == 1L) init.name else "random",
      initial_raw_stress = start$loss * target.rms^2,
      raw_stress = if (failed) NA_real_ else result$loss * target.rms^2,
      normalized_stress = if (failed) NA_real_ else result$loss / length(target),
      iterations = if (inherits(fit, "error")) NA_integer_ else fit$niter,
      converged = identical(reason, "stress_tolerance"), termination = reason,
      backend_stress = if (inherits(fit, "error")) NA_real_ else fit$stress,
      warnings = paste(unique(notices), collapse = " | "),
      error = if (failed) conditionMessage(result) else "",
      stringsAsFactors = FALSE)
    if (!failed && result$loss < best.loss) {
      best <- result
      best.loss <- result$loss
      selected <- run
    }
  }
  runs <- do.call(rbind, records)
  if (is.null(best)) {
    stop("All SMACOF starts failed: ", paste(unique(runs$error), collapse = " | "),
         call. = FALSE)
  }
  if (!runs$converged[selected]) {
    warning("Selected metric.mds() start terminated with ", runs$termination[selected],
            "; inspect metadata$starts", call. = FALSE)
  }
  coords <- best$coords * target.rms
  d <- as.double(stats::dist(best$coords))
  target.scale <- sum(d * target.normalized) / sum(target.normalized^2)
  diag <- if (diagnostics) score.gmds(coords = coords, prepared = prepared,
    scale_mode = scale_mode, distance_floor = distance_floor,
    edge_length_epsilon = edge_length_epsilon, band_quantiles = band_quantiles) else NULL
  gmds.result(coords = coords, method = "metric_mds", prepared = prepared,
    trace = NULL, diagnostics = diag, metadata = list(
      engine = "smacof", grip_version = as.character(getNamespaceVersion("grip")),
      backend_version = as.character(utils::packageVersion("smacof")),
      objective = "raw_distance_stress", pair_weights = "uniform", type = "ratio",
      input_rms_distance = target.rms, coordinate_scale = best$multiplier * target.rms,
      raw_stress = best.loss * target.rms^2,
      target_normalized_rmse = sqrt(best.loss / length(target)),
      stress1_identity = sqrt(best.loss / sum(d^2)),
      stress1_profiled = sqrt(sum((d - target.scale * target.normalized)^2) / sum(d^2)),
      selected_start = selected, converged = runs$converged[selected],
      termination = runs$termination[selected], starts = runs,
      settings = list(init = init.name, n_init = n_init, max_iter = max_iter, eps = eps,
                      seed = seed, dimension = dim)))
}
