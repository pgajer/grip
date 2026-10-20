# Runnable helpers for vignette("metric-mds-animation", package = "grip").
# All fits use the public metric.mds() tracing and continuation interfaces.
surface.sample <- function(n = 250L, seed = 20261020L) {
  set.seed(seed)
  u <- runif(n)
  theta <- runif(n, 0, 2*pi)
  radius <- sqrt(((1 + u*(5^(3/2) - 1))^(2/3) - 1)/4)
  cbind(x = radius*cos(theta), y = radius*sin(theta), z = radius^2)
}

surface.relative.stress <- function(X, target) {
  mean((as.numeric(stats::dist(X))/target - 1)^2)
}

surface.align <- function(X, truth) {
  center <- colMeans(truth)
  centered <- sweep(X, 2, colMeans(X))
  decomposition <- svd(crossprod(centered, sweep(truth, 2, center)))
  sweep(centered %*% decomposition$u %*% t(decomposition$v), 2, center, "+")
}

# One-epoch calls intentionally stop before the rate schedule ends. Do not
# suppress unrelated warnings (invalid targets, numerical failures, etc.).
surface.epoch <- function(...) {
  withCallingHandlers(grip::metric.mds(...), warning = function(w) {
    if (grepl("^Selected metric.mds\\(\\) start terminated with iteration_limit;",
              conditionMessage(w)))
      invokeRestart("muffleWarning")
  })
}

surface.trace <- function(truth, every = 32L, epoch.cap = 10L,
                          max.bytes = 64*1024^2) {
  target <- as.numeric(stats::dist(truth))
  stopifnot(ncol(truth) == 3L, all(is.finite(target)), all(target > 0),
            every >= 1L, epoch.cap >= 1L, epoch.cap <= 500L)
  set.seed(20261021L)
  initial <- matrix(rnorm(nrow(truth)*3L), nrow(truth), 3L)
  initial <- sweep(initial, 2, colMeans(initial))
  ratio <- as.numeric(stats::dist(initial))/target
  initial <- initial * sum(ratio)/sum(ratio^2)
  control <- list(
    keep.state = TRUE, schedule.epochs = 500L, checkpoint.every = 1L,
    scheduler = "hybrid", learning.rate = 0.5, final.rate = 0.01,
    switch.ratio = 0.4,
    trace = list(every = every, max.bytes = max.bytes)
  )
  frames <- list()
  steps <- list()
  pair.count <- length(target)
  fit <- NULL
  for (epoch in seq_len(epoch.cap)) {
    fit <- if (epoch == 1L) surface.epoch(
      distance.matrix = as.matrix(stats::dist(truth)),
      init = initial, dim = 3L, backend = "sgd", n.init = 1L,
      max.iter = 1L, seed = 20261022L, pair.weights = "inverse_squared",
      diagnostics = FALSE, sgd.control = control
    ) else surface.epoch(resume = fit, max.iter = 1L)
    recorded <- fit$update_trace
    stopifnot(recorded$coordinate_units == "input_distance",
              recorded$alignment == "none", recorded$scope == "current_call")
    # Playback preserves finer detail near the start of each epoch. Discard
    # duplicate continuation starts, but always keep the final examined pair.
    within <- recorded$steps$pair_updates - (epoch - 1L)*pair.count
    keep <- (within > 0 & (within <= 32L |
                           (within <= 1024L & within %% 32L == 0L) |
                           within %% 256L == 0L |
                           within == pair.count)) | (epoch == 1L & within == 0)
    frames <- c(frames, recorded$coordinates[keep])
    steps[[epoch]] <- recorded$steps[keep, , drop = FALSE]
    terminal <- tail(recorded$coordinates, 1L)[[1L]]
    if (100*surface.relative.stress(terminal, target) <= 1e-12) break
  }
  errors <- vapply(frames, surface.relative.stress, numeric(1), target = target)
  aligned <- lapply(frames, surface.align, truth = truth)
  stopifnot(all(is.finite(errors)))
  # Check that display alignment preserves distances, including the start.
  for (i in seq_along(frames))
    stopifnot(max(abs(as.numeric(stats::dist(frames[[i]])) -
                       as.numeric(stats::dist(aligned[[i]])))) < 1e-12)
  steps <- do.call(rbind, steps)
  rownames(steps) <- NULL
  stopifnot(!anyDuplicated(steps$pair_updates),
            all(diff(steps$pair_updates) > 0))
  list(truth = truth, initial = initial, raw = frames, aligned = aligned,
       steps = steps, errors = errors, epochs = epoch, pair.count = pair.count,
       control = control, fit = fit, seed = 20261022L,
       stop.reason = if (100*tail(errors, 1) <= 1e-12) "threshold" else "epoch cap")
}

surface.animation <- function(cases) {
  if (!isTRUE(getOption("grip.vignette.interactive", TRUE)) ||
      !requireNamespace("ivue", quietly = TRUE) ||
      !requireNamespace("rgl", quietly = TRUE)) {
    return(knitr::asis_output(paste(
      "\n\nInteractive animation unavailable. Install **ivue** and **rgl**",
      "and rebuild with `options(grip.vignette.interactive = TRUE)`.",
      "Numerical examples and results remain available.\n\n")))
  }
  ids <- list()
  data <- lapply(names(cases), function(key) {
    case <- cases[[key]]
    set.seed(20261023L)
    pairs <- t(combn(nrow(case$truth), 2L))
    chords <- pairs[sample.int(nrow(pairs), min(1000L, nrow(pairs))), , drop = FALSE]
    bounds <- t(apply(rbind(case$truth, do.call(rbind, case$aligned)), 2, range))
    bounds[, 1] <- pmin(bounds[, 1], c(-1, -1, if (key == "saddle") -1 else 0))
    bounds[, 2] <- pmax(bounds[, 2], 1)
    list(key = key, truth = unname(case$truth),
         frames = lapply(case$aligned, unname),
         updates = case$steps$pair_updates, errors = case$errors,
         epochs = case$epochs, pair.count = case$pair.count, chords = chords,
         bounds = bounds, stop.reason = case$stop.reason)
  })
  names(data) <- names(cases)
  for (key in names(data)) {
    limits <- if (key == "saddle") c(-1, 1) else c(0, 1)
    data[[key]]$colors <- grDevices::hcl.colors(256, "Viridis")[
      1 + floor((data[[key]]$truth[, 3] - limits[1])/diff(limits)*255)]
  }
  mesh <- function(sign) {
    point <- function(r, t) c(r*cos(t), r*sin(t), r^2*(cos(t)^2 + sign*sin(t)^2))
    triangles <- list()
    for (j in 0:39) for (i in 0:127) {
      r0 <- j/40; r1 <- (j + 1)/40
      t0 <- i*2*pi/128; t1 <- (i + 1)*2*pi/128
      triangles[[length(triangles) + 1L]] <-
        rbind(point(r0, t0), point(r1, t0), point(r1, t1))
      if (j > 0) triangles[[length(triangles) + 1L]] <-
        rbind(point(r0, t0), point(r1, t1), point(r0, t1))
    }
    vertices <- do.call(rbind, triangles)
    normals <- cbind(-2*vertices[, 1], -2*sign*vertices[, 2], 1)
    list(vertices = vertices, normals = normals/sqrt(rowSums(normals^2)))
  }
  layer <- ivue::layer3D.callback(function(ctx) {
    for (key in names(data)) {
      item <- data[[key]]
      shape <- mesh(if (key == "saddle") -1 else 1)
      ids[[key]] <<- list(
        surface = as.integer(rgl::triangles3d(
          shape$vertices, normals = shape$normals, col = "lightblue",
          alpha = .3, lit = TRUE, specular = "white", shininess = 40,
          front = "filled", back = "filled", polygon_offset = c(1, 1))),
        chords = as.integer(rgl::segments3d(
          item$truth[as.vector(t(item$chords)), ], col = "#375f75",
          alpha = .075, lit = FALSE)),
        truth = as.integer(rgl::spheres3d(item$truth, radius = .01,
                                        col = "#314455", alpha = .65)))
      if (key != names(data)[1L])
        ids[[key]]$points <<- as.integer(rgl::spheres3d(
          item$truth, radius = .02, col = item$colors, lit = TRUE))
    }
  })
  first <- data[[1L]]
  widget <- ivue::plot3D.plain(
    first$truth, col = first$colors, point.type = "sphere", sphere.radius = .02,
    axes = FALSE, layers = list(layer), limits = first$bounds,
    camera = ivue::camera.zup(elevation = 22, turn = -38, fov = 0, zoom = .62),
    width = NULL, height = 490, controls = FALSE,
    description = "Drag to rotate; scroll to zoom. Colors retain original point height."
  )
  ids[[1L]]$points <- attr(widget, "ivue")$draw.ids$object[1L]
  for (key in names(data)) data[[key]]$ids <- ids[[key]]
  script <- system.file("scripts", "metric-mds-animation.js", package = "grip")
  htmlwidgets::onRender(widget, paste(readLines(script), collapse = "\n"),
                        data = list(surfaces = data))
}
