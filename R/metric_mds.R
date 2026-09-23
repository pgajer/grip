grip.mds.has.smacof <- function() {
  requireNamespace("smacof", quietly = TRUE)
}

#' Metric stress MDS using stochastic gradient descent or SMACOF
#'
#' `metric.mds()` minimizes raw distance stress on a graph's
#' all-pairs shortest-path distances using native stochastic gradient descent
#' (SGD, the default) or `smacof::mds(type = "ratio")`.
#' Only the SMACOF backend requires the optional \pkg{smacof} package.
#' Before version 0.2.0.9000,
#' this name performed classical scaling; use [classical.mds()] to retain that
#' behavior. The `add` and `eig` arguments belong to `classical.mds()` only.
#'
#' @details The objective is
#' \deqn{S(Z) = \sum_{i<j}w_{ij}(\|z_i-z_j\|_2-\delta_{ij})^2.}
#' By default, all pair weights are one. With `pair_weights = "inverse_squared"`,
#' \eqn{w_{ij}=1/\delta_{ij}^2}, as in Zheng et al. (2018); stress is then the
#' sum of squared relative distance errors. Edge weights define graph distances;
#' `pair_weights` separately specifies their importance in the stress objective.
#' Both backends, coordinate rescaling, checkpoint/start selection, and stress
#' metadata use the selected weights. The common `score.gmds()` diagnostic panel
#' retains its own definitions and does not become a weighted objective report.
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
#' Newly prepared graph objects supply symmetric strict shortest distances,
#' separately from retained-route lengths used in path diagnostics. Rebuild
#' older cached preparations if their near-tie distance matrix is asymmetric.
#'
#' Targets must be finite, symmetric, and nonnegative, with zero diagonal and
#' at least one positive distance. Missing and infinite distances are rejected.
#' Inverse-squared weighting requires strictly positive off-diagonal targets;
#' zero distances are rejected, not floored. Uniform weighting still allows zero
#' distances. Extreme ratios that overflow or underflow the normalized weights
#' are rejected. Supplied starts may have coincident points but must not be wholly collapsed.
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
#' @param max_iter Positive integer iteration limit per start. For SGD, the
#'   number of complete passes over all unordered pairs. The retained default
#'   is 1000; explicitly request 30 for a short initial SGD trial.
#' @param eps Positive SMACOF tolerance for the change in normalized stress.
#'   Do not supply this argument for SGD, which runs its prescribed schedule.
#' @param seed Integer random seed, or `NULL` to use the current RNG stream.
#'   With a non-NULL seed, random starts do not change the caller's RNG state.
#' @param backend `"sgd"` (default) or `"smacof"`. No automatic fallback occurs.
#' @param pair_weights `"uniform"` (default) or `"inverse_squared"`. The latter
#'   applies the paper's inverse-squared shortest-path-distance weights in either
#'   backend. This is separate from the graph's `edge_weights`.
#' @param sgd_control Named list used only with `backend = "sgd"`:
#'   `scheduler` (`"hybrid"` or `"exponential"`), `learning_rate` (0.5),
#'   `final_rate` (0.01), `switch_ratio` (0.4), `checkpoint_every` (1), and
#'   `max_workspace_bytes` (256 MiB). These defaults are provisional calibration
#'   choices. Rates must be positive with `final_rate <= learning_rate`;
#'   `switch_ratio` is between zero and one inclusive. The memory allowance covers native workspace,
#'   not the R input matrices or total process memory.
#' @return A `"grip_gmds_layout"` object with method `"metric_mds"`.
#'   `metadata` records the objective, backend/version, achieved raw stress,
#'   weighted target-normalized RMSE, both weighted Stress-1 conventions, selected start,
#'   coordinate scale multiplier, and per-start losses and stopping information.
#'   SGD adds per-start elapsed seconds (native fitting including scoring),
#'   native seeds, pair-update counts, best epochs, and checkpoint histories in
#'   `metadata$sgd`. Histories use RMS-normalized target-distance units and include
#'   both raw and scale-profiled stress. Internally, inverse-squared weights are
#'   `1 / (distance / input_rms_distance)^2`; the corresponding stress is already
#'   dimensionless. Uniform raw stress is converted to squared input-distance
#'   units for public metadata. Weighted target-normalized RMSE takes the square
#'   root of the ratio of stress to the sum of weighted squared target distances.
#'   SGD reports `iteration_limit` and `converged = FALSE` when its
#'   schedule finishes, including when an earlier checkpoint is returned.
#' @section SGD behavior:
#' Every pass visits all unordered pairs in shuffled order. The pair step is
#' clipped at `min(rate * normalized_pair_weight, 1)`. Rates use RMS-normalized
#' distance units for both weighting choices; changing input distance units
#' therefore preserves the optimization trajectory up to scale. The default
#' schedule is shared across weighting choices, not a reproduction of the
#' graph-dependent annealing schedule used in Zheng et al. (2018). Symmetric clipped
#' updates follow Zheng et al. (2018); the hybrid exponential/harmonic schedule
#' follows Hangan et al. (2026). The terminal epoch is always scored. The best
#' independently scored checkpoint, including initialization, is retained using
#' scale-profiled stress; checkpoint scoring does not rescale the ongoing search.
#' The R wrapper independently rescales and recomputes its returned stress.
#' The exponential schedule's final rate is its boundary value after the planned
#' passes, not its last applied rate. The hybrid schedule switches after
#' `floor(switch_ratio * max_iter)` passes, at one tenth the initial rate
#' (or the initial rate if the switch is immediate). Its harmonic segment
#' approaches `final_rate` at the boundary only when the switch rate exceeds
#' `final_rate`; otherwise that segment stays constant at the switch rate.
#' With `switch_ratio = 1`, no harmonic segment occurs and `final_rate` is unused.
#' Coincident points with positive target distance are separated along a seeded
#' direction; zero targets are not floored. Without such collisions, updates
#' preserve a deficient initial affine span. Use full-dimensional random starts
#' to examine this sensitivity. This backend retains quadratic pair storage and
#' work per pass. It does not validate the supplied graph distances.
#' @references Zheng, J. X., Pawar, S. and Goodman, D. F. (2018).
#'   Graph Drawing by Stochastic Gradient Descent.
#'   \doi{10.1109/TVCG.2018.2859997}.
#'   Hangan, D., Kobourov, S. and Miller, J. (2026).
#'   Bridging Graph Drawing and Dimensionality Reduction with Stochastic Stress
#'   Optimization. <https://arxiv.org/abs/2605.00641>.
#' @examples
#' # Schedule completion is reported explicitly, not as convergence.
#' fit <- suppressWarnings(metric.mds(edges = edges.cycle(8), n = 8,
#'   dim = 3, init = "random", backend = "sgd", max_iter = 30,
#'   diagnostics = FALSE))
#' fit$metadata$engine
#' fit$metadata$target_normalized_rmse
#' @seealso [classical.mds()], [edge.kk()], [smacof::mds()]
#' @md
#' @export
#' @section Workflow guides:
#' Start with \code{vignette("function-guide", package = "grip")}
#' to choose a layout, diagnostic, or reference comparison. List installed
#' guides with \code{vignette(package = "grip")}.
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
                       band_quantiles = c(1 / 3, 2 / 3),
                       backend = c("sgd", "smacof"),
                       sgd_control = list(),
                       pair_weights = c("uniform", "inverse_squared")) {
  grip.validate.graph.arguments(edges, n, adj_list, weight_list, edge_weights, prepared)
  backend <- match.arg(backend)
  pair_weights <- match.arg(pair_weights)
  if (backend == "sgd" && !missing(eps)) {
    stop("eps is a SMACOF tolerance; SGD uses max_iter and sgd_control", call. = FALSE)
  }
  if (backend == "smacof" && length(sgd_control)) {
    stop("sgd_control requires backend = 'sgd'", call. = FALSE)
  }
  if (backend == "smacof" && !grip.mds.has.smacof()) {
    stop("metric.mds() requires the optional 'smacof' package; install it, ",
         "or use classical.mds() for classical scaling", call. = FALSE)
  }
  if (backend == "smacof" && utils::packageVersion("smacof") < "2.1-7") {
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
  if (backend == "sgd") {
    sgd_control <- grip.mds.sgd.control(sgd_control, max_iter)
    rates <- grip.mds.sgd.rates(sgd_control, max_iter)
  }
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
  if (pair_weights == "inverse_squared" && any(target <= 0)) {
    stop("inverse_squared pair weights require strictly positive off-diagonal distances; use uniform weighting for zero distances", call. = FALSE)
  }
  pair.stiffness <- if (pair_weights == "uniform") rep(1, length(target)) else
    (1 / target.normalized)^2
  if (any(!is.finite(pair.stiffness)) || any(pair.stiffness <= 0)) {
    stop("inverse_squared pair weights are outside the representable numeric range", call. = FALSE)
  }
  # Public inverse-squared stress is dimensionless; normalized-distance weights
  # absorb RMS^2. Uniform stress is converted back to squared input units.
  loss.units <- if (pair_weights == "uniform") target.rms^2 else 1
  target.energy <- sum(pair.stiffness * target.normalized^2)
  weight.matrix <- NULL
  if (backend == "smacof" && pair_weights != "uniform") {
    weight.dist <- stats::as.dist(delta)
    weight.dist[] <- pair.stiffness
    weight.matrix <- as.matrix(weight.dist)
  }
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
  if (!is.null(seed) && (backend == "sgd" || n_init > 1L || identical(init.name, "random"))) {
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
    denominator <- sum(pair.stiffness * d^2)
    if (!is.finite(denominator) || denominator <= 0) {
      stop("MDS configuration is collapsed or has nonfinite distances", call. = FALSE)
    }
    multiplier <- sum(pair.stiffness * d * target.normalized) / denominator
    list(coords = multiplier * x, multiplier = multiplier,
         loss = sum(pair.stiffness * (multiplier * d - target.normalized)^2))
  }
  best <- NULL
  best.loss <- Inf
  records <- vector("list", n_init)
  sgd.records <- vector("list", n_init)
  for (run in seq_len(n_init)) {
    start <- if (run == 1L && !is.null(first)) first else
      matrix(stats::rnorm(prepared$n * dim), nrow = prepared$n, ncol = dim)
    start <- rescale(start)
    notices <- character()
    if (backend == "sgd") native.seed <- sample.int(.Machine$integer.max, 1L)
    if (backend == "sgd") fit.started <- proc.time()[["elapsed"]]
    fit <- tryCatch(withCallingHandlers(
      if (backend == "sgd") {
        grip.mds.sgd.fit(start$coords, target.normalized, rates, sgd_control, native.seed,
                         if (pair_weights == "uniform") NULL else pair.stiffness)
      } else smacof::mds(delta.normalized, ndim = dim, type = "ratio",
                  init = start$coords, itmax = max_iter, eps = eps, weightmat = weight.matrix,
                  principal = FALSE, verbose = FALSE),
      warning = function(w) {
        notices <<- c(notices, conditionMessage(w))
        invokeRestart("muffleWarning")
      }), error = function(e) e)
    if (backend == "sgd") fit.elapsed <- proc.time()[["elapsed"]] - fit.started
    failed <- inherits(fit, "error")
    result <- if (failed) fit else tryCatch(rescale(fit$conf), error = function(e) e)
    failed <- inherits(result, "error")
    rejected <- !failed && result$loss > start$loss + 1e-10 * max(1, start$loss)
    if (rejected) result <- start
    reason <- if (failed) "backend_error" else if (rejected) "rejected_increase" else
      if (backend == "sgd" || fit$niter >= max_iter) "iteration_limit" else "stress_tolerance"
    records[[run]] <- data.frame(
      start = run, initialization = if (run == 1L) init.name else "random",
      initial_raw_stress = start$loss * loss.units,
      raw_stress = if (failed) NA_real_ else result$loss * loss.units,
      normalized_stress = if (failed) NA_real_ else result$loss / target.energy,
      iterations = if (inherits(fit, "error")) NA_integer_ else fit$niter,
      converged = identical(reason, "stress_tolerance"), termination = reason,
      backend_stress = if (inherits(fit, "error")) NA_real_ else fit$stress,
      warnings = paste(unique(notices), collapse = " | "),
      error = if (failed) conditionMessage(result) else "",
      stringsAsFactors = FALSE)
    if (backend == "sgd") {
      records[[run]]$elapsed_seconds <- fit.elapsed
      records[[run]]$native_seed <- native.seed
      records[[run]]$pair_updates <- if (inherits(fit, "error")) NA_real_ else fit$pair_updates
      records[[run]]$best_epoch <- if (inherits(fit, "error")) NA_integer_ else fit$best_epoch
      sgd.records[[run]] <- if (inherits(fit, "error")) list(error = conditionMessage(fit)) else
        list(trace = fit$trace, workspace_bytes = fit$workspace_bytes)
    }
    if (!failed && result$loss < best.loss) {
      best <- result
      best.loss <- result$loss
      selected <- run
    }
  }
  runs <- do.call(rbind, records)
  if (is.null(best)) {
    stop("All ", toupper(backend), " starts failed: ", paste(unique(runs$error), collapse = " | "),
         call. = FALSE)
  }
  if (!runs$converged[selected]) {
    warning("Selected metric.mds() start terminated with ", runs$termination[selected],
            "; inspect metadata$starts", call. = FALSE)
  }
  coords <- best$coords * target.rms
  d <- as.double(stats::dist(best$coords))
  target.scale <- sum(pair.stiffness * d * target.normalized) / target.energy
  diag <- if (diagnostics) score.gmds(coords = coords, prepared = prepared,
    scale_mode = scale_mode, distance_floor = distance_floor,
    edge_length_epsilon = edge_length_epsilon, band_quantiles = band_quantiles) else NULL
  output <- gmds.result(coords = coords, method = "metric_mds", prepared = prepared,
    trace = NULL, diagnostics = diag, metadata = list(
      engine = backend, grip_version = as.character(getNamespaceVersion("grip")),
      backend_version = if (backend == "smacof") as.character(utils::packageVersion("smacof")) else "grip-sgd-mds-v2",
      objective = "raw_distance_stress", pair_weights = pair_weights, type = "ratio",
      input_rms_distance = target.rms, coordinate_scale = best$multiplier * target.rms,
      raw_stress = best.loss * loss.units,
      target_normalized_rmse = sqrt(best.loss / target.energy),
      stress1_identity = sqrt(best.loss / sum(pair.stiffness * d^2)),
      stress1_profiled = sqrt(sum(pair.stiffness * (d - target.scale * target.normalized)^2) / sum(pair.stiffness * d^2)),
      selected_start = selected, converged = runs$converged[selected],
      termination = runs$termination[selected], starts = runs,
      settings = list(init = init.name, n_init = n_init, max_iter = max_iter, eps = eps,
                      seed = seed, dimension = dim, pair_weights = pair_weights)))
  if (backend == "sgd") {
    output$metadata$settings$eps <- NULL
    output$metadata$settings$sgd_control <- sgd_control
    output$metadata$sgd <- sgd.records
    output$metadata$settings$randomization <- "R per-start seed; mt19937_64 rejection-sampled Fisher-Yates"
  }
  output
}
