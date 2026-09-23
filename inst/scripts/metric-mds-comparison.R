# Public recipes for the metric-MDS comparison vignette. No work runs on source().
# All fitting and graph preparation use exported grip interfaces.
comparison_cases <- function(sides = c(8L, 12L)) {
  cases <- list()
  for (side in sides) {
    grid <- seq(-1, 1, length.out = side)
    uv <- as.matrix(expand.grid(x = grid, y = grid))
    n <- nrow(uv)
    # This is the explicit quadratic lift (x, y, a*x^2 + b*y^2).
    # No fitted surface geodesics or frozen external fixtures are used.
    clouds <- list(plane = cbind(uv, 0),
                   paraboloid = cbind(uv, .8 * (uv[, 1]^2 + uv[, 2]^2)),
                   saddle = cbind(uv, .8 * (uv[, 1]^2 - uv[, 2]^2)))
    clouds$gaussian <- withr::with_seed(701L + n, matrix(rnorm(n * 3L), n, 3L))
    theta <- seq(0, 3 * pi, length.out = n)
    clouds$helix <- cbind(cos(theta), sin(theta), theta / pi)
    complete <- t(utils::combn(n, 2L))
    index <- matrix(seq_len(n), side, side)
    edges <- faces <- list()
    for (j in seq_len(side - 1L)) for (i in seq_len(side - 1L)) {
      a <- index[i, j]; b <- index[i + 1L, j]
      c <- index[i, j + 1L]; d <- index[i + 1L, j + 1L]
      edges[[length(edges) + 1L]] <- rbind(c(a,b), c(a,c), c(b,d), c(c,d), c(a,d), c(b,c))
      faces[[length(faces) + 1L]] <- rbind(c(a,b,d), c(a,d,c))
    }
    mesh.edges <- unique(t(apply(do.call(rbind, edges), 1, sort)))
    triangles <- do.call(rbind, faces)
    for (kind in c('gaussian', 'plane', 'paraboloid', 'saddle',
                   'paraboloid_graph', 'saddle_graph', 'helix_graph')) {
      family <- sub('_graph$', '', kind)
      graph <- grepl('_graph$', kind)
      X <- unname(clouds[[family]])
      E <- if (!graph) complete else if (family == 'helix') {
        unique(rbind(cbind(seq_len(n - 1L), 2:n),
                     cbind(seq_len(n - 7L), 8:n)))
      } else mesh.edges
      weights <- sqrt(rowSums((X[E[,1], , drop = FALSE] - X[E[,2], , drop = FALSE])^2))
      labels <- c(gaussian = 'Gaussian cloud', plane = 'Plane', paraboloid = 'Paraboloid',
                  saddle = 'Saddle', helix = 'Spatial helix graph')
      id <- paste(kind, n, sep = '-')
      cases[[id]] <- list(id = id, label = unname(labels[[family]]), n = n,
        dimension = if (family == 'plane') 2L else 3L,
        target = if (graph) 'Local graph paths' else 'Ambient Euclidean',
        X = X, edges = E, weights = weights,
        draw_edges = if (graph) E else matrix(integer(), 0L, 2L),
        triangles = if (family %in% c('plane','paraboloid','saddle')) triangles else NULL)
    }
  }
  cases
}

comparison_score <- function(coords, distances) {
  residual <- as.double(stats::dist(coords)) - as.double(stats::as.dist(distances))
  raw <- sum(residual^2)
  c(raw_stress = raw, error = sqrt(raw / sum(stats::as.dist(distances)^2)))
}

comparison_align <- function(coords, reference) {
  # Orthogonal Procrustes: permit reflection, never change scale.
  if (ncol(coords) == 2L) reference <- reference[, 1:2, drop = FALSE]
  center <- colMeans(reference)
  A <- sweep(coords, 2L, colMeans(coords))
  B <- sweep(reference, 2L, center)
  s <- svd(crossprod(A, B))
  sweep(A %*% (s$u %*% t(s$v)), 2L, center, '+')
}

comparison_start <- function(case, seed) {
  withr::with_seed(1000L + seed, matrix(rnorm(case$n * case$dimension), case$n, case$dimension))
}

comparison_fit <- function(case, prepared, seed, budget, backend,
                           fit_fun = grip::metric.mds) {
  warnings <- character()
  args <- list(prepared = prepared, dim = case$dimension,
               init = comparison_start(case, seed), n_init = 1L, seed = seed,
               max_iter = budget, backend = backend, diagnostics = FALSE)
  if (backend == 'smacof') args$eps <- 1e-8
  if (backend == 'sgd') args$sgd_control <- list(scheduler = 'hybrid', learning_rate = .5,
    final_rate = .01, switch_ratio = .4, checkpoint_every = 1L, max_workspace_bytes = 256 * 1024^2)
  started <- proc.time()[['elapsed']]
  fit <- tryCatch(withCallingHandlers(do.call(fit_fun, args), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart('muffleWarning')
  }), error = function(e) e)
  seconds <- proc.time()[['elapsed']] - started
  row <- data.frame(case = case$id, seed = seed, budget = budget, backend = backend,
    status = 'error', seconds = seconds, raw_stress = NA_real_, error = NA_real_,
    iterations = NA_integer_, termination = NA_character_, converged = FALSE,
    warnings = paste(warnings, collapse = ' | '), failure = '', stringsAsFactors = FALSE)
  if (inherits(fit, 'error')) {
    row$failure <- conditionMessage(fit)
    return(list(row = row, coords = NULL))
  }
  score <- comparison_score(fit$coords, prepared$distance_matrix)
  stopifnot(isTRUE(all.equal(unname(score[['raw_stress']]), fit$metadata$raw_stress,
                            tolerance = 1e-8)), all(is.finite(fit$coords)))
  row$status <- 'ok'; row$raw_stress <- score[['raw_stress']]; row$error <- score[['error']]
  row$iterations <- fit$metadata$starts$iterations[[1]]
  row$termination <- fit$metadata$termination; row$converged <- fit$metadata$converged
  list(row = row, coords = fit$coords)
}

comparison_rows <- function(bundle) do.call(rbind, lapply(bundle$fits, `[[`, 'row'))

comparison_summary <- function(bundle) {
  rows <- comparison_rows(bundle)
  keys <- unique(rows[, c('case','budget','backend')])
  do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    key <- keys[i, ]
    z <- rows[rows$case == key$case & rows$budget == key$budget & rows$backend == key$backend, ]
    ok <- z$status == 'ok'
    q <- function(x, p) if (any(ok)) unname(stats::quantile(x[ok], p)) else NA_real_
    data.frame(key, successful = sum(ok), attempted = nrow(z),
      error = q(z$error, .5), error_q25 = q(z$error, .25), error_q75 = q(z$error, .75),
      seconds = q(z$seconds, .5), seconds_q25 = q(z$seconds, .25), seconds_q75 = q(z$seconds, .75))
  }))
}

comparison_representative <- function(bundle, id, seed = 1L, budget = 100L) {
  result <- lapply(c('sgd', 'smacof'), function(backend) {
    # Exact fixed settings, without selecting a winner or replacing a failed run.
    matches <- which(vapply(bundle$fits, function(f) f$row$case == id &&
      f$row$seed == seed && f$row$budget == budget && f$row$backend == backend, logical(1)))
    stopifnot(length(matches) == 1L)
    bundle$fits[[matches]]
  })
  stats::setNames(result, c('sgd', 'smacof'))
}

comparison_view <- function(case, fits, width = 900L, height = 460L,
                            show = c('overlay', 'reference', 'sgd', 'smacof'),
                            legend = TRUE, controls = TRUE) {
  show <- match.arg(show)
  stopifnot(case$dimension == 3L)
  coordinates <- list(Reference = case$X)
  for (backend in c('sgd','smacof')) if (!is.null(fits[[backend]]$coords)) {
    coordinates[[toupper(backend)]] <- comparison_align(fits[[backend]]$coords, case$X)
  }
  palette <- c(Reference = '#888888', SGD = '#1769AA', SMACOF = '#D66A19')
  bounds <- t(apply(do.call(rbind, coordinates), 2L, range))
  padding <- pmax(bounds[,2] - bounds[,1], 1e-8) * .04
  bounds <- bounds + cbind(-padding, padding)
  if (show != 'overlay') coordinates <- coordinates[intersect(names(coordinates),
    if (show == 'reference') 'Reference' else toupper(show))]
  if (!length(coordinates)) stop('This fit is unavailable; inspect its recorded failure.')
  X <- do.call(rbind, coordinates)
  rownames(X) <- NULL
  layers <- list()
  for (j in seq_along(coordinates)) {
    name <- names(coordinates)[j]; offset <- (j - 1L) * case$n
    if (nrow(case$draw_edges)) layers[[length(layers) + 1L]] <-
      ivue::layer3D.edges(case$draw_edges + offset, col = grDevices::adjustcolor(palette[[name]], .38), width = .7)
    if (!is.null(case$triangles)) layers[[length(layers) + 1L]] <-
      ivue::layer3D.mesh(case$triangles + offset, col = palette[[name]], alpha = .09,
                         edges = FALSE)
  }
  layers[[length(layers) + 1L]] <- ivue::layer3D.axes(limits = bounds, padding = 0,
    width = 1, head.length = .035, cex = .8)
  groups <- rep(names(coordinates), each = case$n)
  widget <- ivue::plot3D.groups(X, groups = groups,
    scale = ivue::color.scale.groups(names(palette), colors = palette),
    point.type = 'point', point.size = 4, alpha = .8, axes = FALSE,
    xlab = '', ylab = '', zlab = '', aspect = 'equal',
    camera = ivue::camera.zup(elevation = 22, turn = -125,
                            zoom = .85),
    limits = bounds, legend.show = legend, controls = controls,
    layers = layers, width = width, height = height, legend.width = 150,
    description = if (controls) paste(case$label, case$n, 'points;', case$target,
      'targets. Gray: reference; blue: SGD; orange: SMACOF. Rigid alignment only.') else NULL)
  htmlwidgets::onRender(widget, "function(el) {
    el.querySelectorAll('.ivue-legend details').forEach(function(node) { node.remove(); });
    el.querySelectorAll('.ivue-legend').forEach(function(node) { node.style.width = 'max-content'; });
  }")
}

comparison_planar <- function(case, fits) {
  old <- graphics::par(mfrow = c(1,3), mar = c(3,3,3,1)); on.exit(graphics::par(old))
  Z <- list(Reference = case$X[,1:2])
  for (b in c('sgd','smacof')) Z[[toupper(b)]] <- if (is.null(fits[[b]]$coords)) NULL else
    comparison_align(fits[[b]]$coords, case$X)
  limits <- range(unlist(Z), finite = TRUE)
  for (name in c('Reference','SGD','SMACOF')) {
    if (is.null(Z[[name]])) { graphics::plot.new(); graphics::title(paste(name, 'unavailable')); next }
    graphics::plot(Z[[name]], asp = 1, xlim = limits, ylim = limits, pch = 16, cex = .65,
      col = c(Reference = '#888888', SGD = '#1769AA', SMACOF = '#D66A19')[[name]],
      xlab = 'x', ylab = 'y', main = name)
  }
}

comparison_tradeoff <- function(bundle, ids) {
  summary <- comparison_summary(bundle)
  old <- graphics::par(mfrow = c(1, length(ids)), mar = c(4,4,3,1)); on.exit(graphics::par(old))
  for (id in ids) {
    s <- summary[summary$case == id & summary$successful > 0L, ]
    if (!nrow(s)) { graphics::plot.new(); graphics::title(paste(id, 'no completed fits')); next }
    x <- pmax(s$seconds, .001); y <- pmax(s$error, 1e-12)
    graphics::plot(x, y, log = 'xy', type = 'n',
      xlim = range(c(pmax(s$seconds_q25, .001), pmax(s$seconds_q75, .001))) * c(.65,1.5),
      ylim = range(c(pmax(s$error_q25, 1e-12), pmax(s$error_q75, 1e-12))) * c(.5,2),
      xlab = 'Fit seconds (log scale)', ylab = 'Normalized error (log scale)',
      main = paste(bundle$cases[[id]]$n, 'points'))
    for (b in c('sgd','smacof')) {
      z <- s[s$backend == b, ]; col <- if (b == 'sgd') '#1769AA' else '#D66A19'
      xx <- pmax(z$seconds, .001); yy <- pmax(z$error, 1e-12)
      graphics::segments(pmax(z$seconds_q25,.001), yy, pmax(z$seconds_q75,.001), yy, col = col)
      graphics::segments(xx, pmax(z$error_q25,1e-12), xx, pmax(z$error_q75,1e-12), col = col)
      graphics::points(xx, yy, pch = if (b == 'sgd') 16 else 17, col = col)
      graphics::text(xx, yy, labels = z$budget, pos = if (b == 'sgd') 3 else 1, cex = .7, col = col)
    }
    graphics::legend('topright', c('SGD','SMACOF'), col = c('#1769AA','#D66A19'), pch = c(16,17), bty = 'n', cex = .8)
  }
}
