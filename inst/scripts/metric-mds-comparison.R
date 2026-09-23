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

comparison_score <- function(coords, distances, pair_weights = 'uniform') {
  target <- as.double(stats::as.dist(distances))
  residual <- as.double(stats::dist(coords)) - target
  if (pair_weights == 'inverse_squared') {
    if (any(target <= 0)) stop('Inverse-squared scoring needs positive distances.')
    raw <- sum((residual / target)^2)
    return(c(raw_stress = raw, error = sqrt(raw / length(target))))
  }
  stopifnot(pair_weights == 'uniform')
  raw <- sum(residual^2)
  c(raw_stress = raw, error = sqrt(raw / sum(target^2)))
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
                           fit_fun = grip::metric.mds, pair_weights = 'uniform') {
  warnings <- character()
  args <- list(prepared = prepared, dim = case$dimension,
               init = comparison_start(case, seed), n_init = 1L, seed = seed,
               max_iter = budget, backend = backend, pair_weights = pair_weights, diagnostics = FALSE)
  if (backend == 'smacof') args$eps <- 1e-8
  if (backend == 'sgd') args$sgd_control <- list(scheduler = 'hybrid', learning_rate = .5,
    final_rate = .01, switch_ratio = .4, checkpoint_every = 1L, max_workspace_bytes = 256 * 1024^2)
  started <- proc.time()[['elapsed']]
  fit <- tryCatch(withCallingHandlers(do.call(fit_fun, args), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart('muffleWarning')
  }), error = function(e) e)
  seconds <- proc.time()[['elapsed']] - started
  row <- data.frame(case = case$id, seed = seed, budget = budget, backend = backend,
    pair_weights = pair_weights, uniform_error = NA_real_, relative_error = NA_real_,
    status = 'error', seconds = seconds, raw_stress = NA_real_, error = NA_real_,
    iterations = NA_integer_, termination = NA_character_, converged = FALSE,
    warnings = paste(warnings, collapse = ' | '), failure = '', stringsAsFactors = FALSE)
  if (inherits(fit, 'error')) {
    row$failure <- conditionMessage(fit)
    return(list(row = row, coords = NULL))
  }
  score <- comparison_score(fit$coords, prepared$distance_matrix, pair_weights)
  stopifnot(isTRUE(all.equal(unname(score[['raw_stress']]), fit$metadata$raw_stress,
                            tolerance = 1e-8)), all(is.finite(fit$coords)))
  row$uniform_error <- comparison_score(fit$coords, prepared$distance_matrix)[['error']]
  row$relative_error <- if (all(prepared$distance_matrix[lower.tri(prepared$distance_matrix)] > 0))
    comparison_score(fit$coords, prepared$distance_matrix, 'inverse_squared')[['error']] else NA_real_
  row$status <- 'ok'; row$raw_stress <- score[['raw_stress']]; row$error <- score[['error']]
  row$iterations <- fit$metadata$starts$iterations[[1]]
  row$termination <- fit$metadata$termination; row$converged <- fit$metadata$converged
  list(row = row, coords = fit$coords, metadata = fit$metadata)
}

comparison_rows <- function(bundle) {
  do.call(rbind, lapply(bundle$fits, function(fit) {
    row <- fit$row
    if (is.null(row$pair_weights)) row$pair_weights <- 'uniform'
    row
  }))
}

comparison_summary <- function(bundle) {
  rows <- comparison_rows(bundle)
  keys <- unique(rows[, c('case','budget','backend','pair_weights')])
  do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    key <- keys[i, ]
    z <- rows[rows$case == key$case & rows$budget == key$budget & rows$backend == key$backend & rows$pair_weights == key$pair_weights, ]
    ok <- z$status == 'ok'
    q <- function(x, p) if (any(ok)) unname(stats::quantile(x[ok], p)) else NA_real_
    data.frame(key, successful = sum(ok), attempted = nrow(z),
      error = q(z$error, .5), error_q25 = q(z$error, .25), error_q75 = q(z$error, .75),
      uniform_error = if ('uniform_error' %in% names(z)) q(z$uniform_error, .5) else NA_real_,
      relative_error = if ('relative_error' %in% names(z)) q(z$relative_error, .5) else NA_real_,
      seconds = q(z$seconds, .5), seconds_q25 = q(z$seconds, .25), seconds_q75 = q(z$seconds, .75))
  }))
}

comparison_representative <- function(bundle, id, seed = 1L, budget = 100L, pair_weights = 'uniform') {
  result <- lapply(c('sgd', 'smacof'), function(backend) {
    # Exact fixed settings, without selecting a winner or replacing a failed run.
    matches <- which(vapply(bundle$fits, function(f) f$row$case == id &&
      f$row$seed == seed && f$row$budget == budget && f$row$backend == backend &&
      (if (is.null(f$row$pair_weights)) 'uniform' else f$row$pair_weights) == pair_weights, logical(1)))
    stopifnot(length(matches) == 1L)
    bundle$fits[[matches]]
  })
  stats::setNames(result, c('sgd', 'smacof'))
}

comparison_view <- function(case, fits, width = 900L, height = 460L,
                            show = c('overlay', 'reference', 'sgd', 'smacof'),
                            legend = TRUE, controls = TRUE, weighting = FALSE, limits = NULL,
                            series_labels = NULL, series_colors = NULL,
                            reference = TRUE, description = NULL, layout_selector = FALSE) {
  show <- match.arg(show)
  stopifnot(case$dimension == 3L)
  coordinates <- if (reference) list(Reference = case$X) else list()
  labels <- if (!is.null(series_labels)) series_labels else if (weighting) c(sgd_uniform = 'SGD: uniform',
    sgd_inverse_squared = 'SGD: inverse-squared', smacof_uniform = 'SMACOF: uniform',
    smacof_inverse_squared = 'SMACOF: inverse-squared') else c(sgd = 'SGD', smacof = 'SMACOF')
  for (backend in names(labels)) if (!is.null(fits[[backend]]$coords)) {
    coordinates[[labels[[backend]]]] <- comparison_align(fits[[backend]]$coords, case$X)
  }
  palette <- if (!is.null(series_colors)) series_colors else if (weighting) c(Reference = '#888888', 'SGD: uniform' = '#1769AA',
    'SGD: inverse-squared' = '#009E73', 'SMACOF: uniform' = '#D66A19',
    'SMACOF: inverse-squared' = '#AA3377') else
      c(Reference = '#888888', SGD = '#1769AA', SMACOF = '#D66A19')
  bounds <- t(apply(do.call(rbind, coordinates), 2L, range))
  padding <- pmax(bounds[,2] - bounds[,1], 1e-8) * .04
  bounds <- bounds + cbind(-padding, padding)
  if (!is.null(limits)) bounds <- limits
  if (show != 'overlay') coordinates <- coordinates[intersect(names(coordinates),
    if (show == 'reference') 'Reference' else toupper(show))]
  if (!length(coordinates)) stop('This fit is unavailable; inspect its recorded failure.')
  X <- do.call(rbind, coordinates)
  rownames(X) <- NULL
  object.ids <- list()
  if (layout_selector) {
    # Track points, edges and surfaces separately for each selectable layout.
    layers <- list(ivue::layer3D.callback(function(ctx) {
      rgl::pop3d(id = unique(ctx$draw.ids$object))
      for (name in names(coordinates)) {
        coords <- coordinates[[name]]
        ids <- rgl::points3d(coords, col = palette[[name]], size = 4,
                             alpha = .8, lit = FALSE)
        if (nrow(case$draw_edges)) ids <- c(ids,
          rgl::segments3d(coords[as.vector(t(case$draw_edges)), , drop = FALSE],
            col = palette[[name]], alpha = .38, lwd = .7, lit = FALSE))
        if (!is.null(case$triangles)) ids <- c(ids,
          rgl::triangles3d(coords[as.vector(t(case$triangles)), , drop = FALSE],
            col = palette[[name]], alpha = .09, lit = FALSE))
        object.ids[[name]] <<- as.integer(ids)
      }
    }))
  } else {
    layers <- list()
    for (j in seq_along(coordinates)) {
      name <- names(coordinates)[j]; offset <- (j - 1L) * case$n
      if (nrow(case$draw_edges)) layers[[length(layers) + 1L]] <-
        ivue::layer3D.edges(case$draw_edges + offset, col = grDevices::adjustcolor(palette[[name]], .38), width = .7)
      if (!is.null(case$triangles)) layers[[length(layers) + 1L]] <-
        ivue::layer3D.mesh(case$triangles + offset, col = palette[[name]], alpha = .09,
                           edges = FALSE)
    }
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
    limits = bounds, legend.show = legend, controls = controls && !layout_selector,
    layers = layers, width = width, height = height, legend.width = 150,
    description = if (!is.null(description)) description else if (controls && weighting) paste(case$label, case$n,
      'vertices. Gray: reference; blue/green: SGD with uniform/inverse-squared weights;',
      'orange/purple: SMACOF with uniform/inverse-squared weights. Alignment preserves scale.') else if (controls) paste(case$label, case$n, 'points;', case$target,
      'targets. Gray: reference; blue: SGD; orange: SMACOF. Rigid alignment only.') else NULL)
  widget$sizingPolicy$browser$fill <- FALSE
  widget <- htmlwidgets::onRender(widget, "function(el) {
    // Keep the WebGL buffer and displayed canvas within the widget dimensions.
    if (el.gripResize) el.gripResize.disconnect();
    function fitCanvas() {
      var scene = el.rglinstance, width = el.clientWidth, height = el.clientHeight;
      if (width > 0 && height > 0 && (scene.canvas.width !== width || scene.canvas.height !== height)) {
        el.width = width; el.height = height; scene.resize(el); scene.drawScene();
      }
    }
    el.gripResize = new ResizeObserver(fitCanvas);
    el.gripResize.observe(el); el.gripResize.observe(el.querySelector('canvas'));
    fitCanvas();
    el.querySelectorAll('.ivue-legend details').forEach(function(node) { node.remove(); });
    el.querySelectorAll('.ivue-legend').forEach(function(node) { node.style.width = 'max-content'; });
    var oldPanel = el.querySelector('.ivue-tools');
    if (oldPanel) {
      var toolbar = document.createElement('div'); toolbar.className = 'ivue-tools';
      toolbar.setAttribute('role', 'group'); toolbar.setAttribute('aria-label', '3D view controls');
      toolbar.style.border = '0'; toolbar.style.maxHeight = 'none'; toolbar.style.overflow = 'visible';
      Array.from(oldPanel.children).forEach(function(child) {
        if (child.tagName !== 'SUMMARY') toolbar.appendChild(child);
      });
      var hint = toolbar.querySelector('p');
      if (hint) { hint.textContent = 'Drag to rotate; scroll to zoom.'; hint.style.margin = '2px 0'; }
      ['pointerdown', 'mousedown', 'touchstart', 'wheel'].forEach(function(event) {
        toolbar.addEventListener(event, function(e) { e.stopPropagation(); });
      });
      toolbar.style.position = 'static'; toolbar.style.boxSizing = 'border-box';
      toolbar.style.width = el.style.width; toolbar.style.maxWidth = '100%';
      toolbar.style.margin = '0 auto 8px';
      oldPanel.remove();
      if (el.gripControls) el.gripControls.remove();
      el.gripControls = toolbar; el.after(toolbar);
    }
  }")
  if (layout_selector && controls)
    widget <- comparison_layout_selector(widget, object.ids)
  widget
}

# Match the graph gallery's always-visible selector; changing the selection
# changes only object visibility, leaving alignment, scale and camera intact.
comparison_layout_selector <- function(widget, object.ids) {
  htmlwidgets::onRender(widget, "function(el, x, data) {
    var scene = el.rglinstance, root = scene.scene.rootSubscene;
    var legend = el.querySelector('.ivue-legend'), rows = {};
    if (legend) Object.keys(data.ids).forEach(function(name) {
      rows[name] = Array.from(legend.children).find(function(node) {
        return node.tagName === 'DIV' && node.textContent.trim().startsWith(name + ' (');
      });
    });
    var panel = document.createElement('div'); panel.className = 'ivue-tools';
    panel.setAttribute('role', 'group'); panel.setAttribute('aria-label', 'Layout controls');
    panel.style.border = '0'; panel.style.maxHeight = 'none'; panel.style.overflow = 'visible';
    var label = document.createElement('label'); label.textContent = 'Layout: ';
    var select = document.createElement('select');
    select.className = 'grip-layout-selector'; select.style.font = 'inherit';
    Object.keys(data.ids).filter(function(name) { return name !== 'Reference'; }).forEach(function(name) {
      var option = document.createElement('option'); option.value = name;
      option.textContent = name; select.appendChild(option);
    });
    var all = document.createElement('option'); all.value = 'all'; all.textContent = 'all layouts';
    select.appendChild(all); select.value = 'all';
    label.appendChild(select); panel.appendChild(label);
    var hint = document.createElement('span'); hint.textContent = ' Drag to rotate; scroll to zoom.';
    panel.appendChild(hint);
    function update() {
      Object.keys(data.ids).forEach(function(name) {
        var visible = name === 'Reference' || select.value === 'all' || select.value === name;
        data.ids[name].forEach(function(id) {
          if (visible) scene.addToSubscene(id, root); else scene.delFromSubscene(id, root);
        });
        if (rows[name]) rows[name].style.display = visible ? 'flex' : 'none';
      });
      scene.drawScene();
    }
    select.addEventListener('change', update);
    ['pointerdown', 'mousedown', 'touchstart', 'wheel'].forEach(function(event) {
      panel.addEventListener(event, function(e) { e.stopPropagation(); });
    });
    panel.style.position = 'static'; panel.style.boxSizing = 'border-box';
    panel.style.width = el.style.width; panel.style.maxWidth = '100%';
    panel.style.margin = '0 auto 8px';
    if (el.gripControls) el.gripControls.remove();
    el.gripControls = panel; el.after(panel); update();
  }", data = list(ids = object.ids))
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


comparison_weighting_fits <- function(bundle, id, seed = 1L) {
  fits <- list()
  for (weight in c('uniform','inverse_squared')) {
    pair <- comparison_representative(bundle, id, seed = seed, pair_weights = weight)
    for (backend in names(pair)) fits[[paste(backend,weight,sep='_')]] <- pair[[backend]]
  }
  fits[c('sgd_uniform','sgd_inverse_squared','smacof_uniform','smacof_inverse_squared')]
}

comparison_weighting_panels <- function(case, fits) {
  overlay <- comparison_view(case, fits, weighting = TRUE, controls = FALSE)
  limits <- attr(overlay,'ivue')$limits
  panels <- lapply(names(fits), function(key) {
    fit <- fits[[key]]; backend <- strsplit(key,'_',fixed=TRUE)[[1]][1]
    label <- paste(toupper(backend), if (fit$row$pair_weights=='uniform') 'uniform' else 'inverse-squared')
    htmltools::tags$section(htmltools::tags$h3(label),
      htmltools::tags$p(sprintf('Uniform error %.3g; relative error %.3g', fit$row$uniform_error, fit$row$relative_error)),
      if (is.null(fit$coords)) htmltools::tags$p('Fit unavailable') else
        comparison_view(case, stats::setNames(list(fit),key), weighting=TRUE,
          width=450L, height=320L, legend=FALSE, controls=FALSE, limits=limits))
  })
  htmltools::tagList(htmltools::tags$style(htmltools::HTML(
    'body {margin:8px;font:15px system-ui;width:900px;} .panels {display:grid;grid-template-columns:450px 450px;} h3,p {text-align:center;margin:4px;} section {width:450px;} .panels .ivue-description {display:none;}')),
    htmltools::tags$div(class='panels',panels))
}
