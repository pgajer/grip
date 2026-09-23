# Public, installed recipes for the Graph examples vignette. No fitting on source().
graph_example_cases <- function() {
  data <- new.env()
  utils::data('zheng.graphs', package = 'grip', envir = data)
  generated <- list(
    tree = grip::kary.tree.weighted.graph(k = 2, depth = 5, depth_decay = .85,
                                         branch_spread = .2, normalize = 'none'),
    carpet = grip::recursive.mask.grid.surface.graph(grip::mask.border(3), level = 2,
                surface = 'saddle', amplitude = .5, normalize = 'none'),
    saddle = grip::mesh.surface.graph(10, 10, surface = 'saddle', connectivity = 'diagonal', normalize = 'none'),
    paraboloid = grip::mesh.surface.graph(10, 10, surface = 'paraboloid', connectivity = 'diagonal', normalize = 'none'),
    torus = grip::torus.surface.graph(12, 8, normalize = 'none'),
    porous = grip::cube.channel.network.surface.graph(side = 5, level = 1,
                surface = 'twisted', twist = .2, normalize = 'none'))
  labels <- c(tree = 'A branching tree', carpet = 'A recursive carpet on a saddle',
    saddle = 'A saddle', paraboloid = 'A paraboloid', torus = 'A torus',
    porous = 'A porous cube', karate = 'Karate club', dwt_66 = 'dwt_66: a small structural graph',
    lesmis = 'Les Miserables', dwt_307 = 'dwt_307: a structural graph',
    '494_bus' = '494_bus: a power network', dwt_1005 = 'dwt_1005: a structural mesh',
    '1138_bus' = '1138_bus: a power network')
  karate.edges <- as.matrix(utils::read.csv(system.file('extdata', 'karate-club-edges.csv', package = 'grip')))
  cases <- c(generated, list(karate = list(n = 34L, edges = karate.edges,
                  edge_weights = rep(1, nrow(karate.edges)))), data$zheng.graphs)
  for (id in names(cases)) {
    cases[[id]]$id <- id; cases[[id]]$label <- unname(labels[id])
    cases[[id]]$group <- if (id %in% names(generated)) 'Generated graphs' else 'Graphs from applications'
    cases[[id]]$lengths <- if (id %in% names(generated)) 'Supplied edge lengths' else 'Unit edge lengths'
  }
  cases
}

graph_example_scores <- function(coords, graph, distances) {
  edges <- graph$edges
  drawn <- sqrt(rowSums((coords[edges[,1], , drop = FALSE] - coords[edges[,2], , drop = FALSE])^2))
  target <- as.double(stats::as.dist(distances))
  c(graph_error = sqrt(sum((as.double(stats::dist(coords)) - target)^2) / sum(target^2)),
    edge_error = sqrt(sum((drawn - graph$edge_weights)^2) / sum(graph$edge_weights^2)))
}

graph_example_fit <- function(graph, seed = 2026L) {
  prepared <- grip::prepare.graph.geodesic.mds(edges = graph$edges, n = graph$n,
                                               edge_weights = graph$edge_weights)
  stopifnot(all(is.finite(prepared$distance_matrix)))
  records <- list(); warnings <- character()
  capture <- function(expr) withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart('muffleWarning')
  })
  started <- proc.time()[['elapsed']]
  mds <- capture(grip::metric.mds(prepared = prepared, dim = 3, init = 'random',
    n_init = 3, max_iter = 80, seed = seed, diagnostics = FALSE))
  records$mds <- list(coords = mds$coords, metadata = mds$metadata,
    seconds = proc.time()[['elapsed']] - started, warnings = warnings,
    scores = graph_example_scores(mds$coords, graph, prepared$distance_matrix))
  warnings <- character(); started <- proc.time()[['elapsed']]
  refined <- tryCatch(capture(grip::edge.kk(coords = mds$coords, prepared = prepared,
    dim = 3, stiffness_method = 'uniform', density_mix_schedule = 1,
    scale_mode = 'identity', max_iter = 100, return_trace = FALSE, diagnostics = FALSE)),
    error = function(e) e)
  if (inherits(refined, 'error')) {
    records$refined <- list(error = conditionMessage(refined), warnings = warnings)
  } else {
    records$refined <- list(coords = refined$coords, metadata = refined$metadata,
      seconds = proc.time()[['elapsed']] - started, warnings = warnings,
      scores = graph_example_scores(refined$coords, graph, prepared$distance_matrix))
  }
  records
}

graph_example_align <- function(coords, reference) {
  A <- sweep(coords, 2, colMeans(coords)); B <- sweep(reference, 2, colMeans(reference))
  s <- svd(crossprod(A, B))
  sweep(A %*% (s$u %*% t(s$v)), 2, colMeans(reference), '+')
}

graph_example_view <- function(graph, fits, show = c('overlay', 'mds', 'refined', 'reference'),
                               width = 900L, height = 460L, controls = TRUE) {
  show <- match.arg(show)
  reference <- graph$coords_surface
  if (is.null(reference)) reference <- fits$mds$coords
  if (is.null(reference)) stop('No available coordinates for this graph.')
  configurations <- list()
  if (!is.null(graph$coords_surface)) configurations$reference <- reference
  for (method in c('mds','refined')) if (!is.null(fits[[method]]$coords)) {
    configurations[[method]] <- graph_example_align(fits[[method]]$coords, reference)
  }
  bounds <- t(apply(do.call(rbind, configurations), 2, range))
  pad <- pmax(bounds[,2] - bounds[,1], 1e-6) * .06
  bounds <- bounds + cbind(-pad, pad)
  if (show != 'overlay') configurations <- configurations[intersect(names(configurations), show)]
  if (!length(configurations)) stop('Requested layout is unavailable.')
  palette <- c(reference = '#888888', mds = '#1769AA', refined = '#D66A19')
  labels <- c(reference = 'Generating coordinates', mds = 'metric-MDS', refined = 'metric-MDS + edge-KK')
  # Separate objects let the selector hide both vertices and edges without
  # rebuilding the scene or changing its camera, bounds, or alignment.
  object.ids <- list()
  layers <- list(ivue::layer3D.callback(function(ctx) {
    rgl::pop3d(id = unique(ctx$draw.ids$object))
    for (name in names(configurations)) {
      coords <- configurations[[name]]
      points <- rgl::points3d(coords, col = palette[[name]],
                             size = if (graph$n > 500) 2.5 else 4, lit = FALSE)
      edges <- rgl::segments3d(coords[as.vector(t(graph$edges)), , drop = FALSE],
                              col = palette[[name]], alpha = .5, lwd = .8, lit = FALSE)
      object.ids[[name]] <<- as.integer(c(points, edges))
    }
  }))
  layers[[length(layers) + 1L]] <- ivue::layer3D.axes(limits = bounds, padding = 0,
    width = 1, head.length = .035, cex = .8)
  names <- names(configurations)
  widget <- ivue::plot3D.groups(do.call(rbind, configurations),
    groups = rep(unname(labels[names]), each = graph$n),
    scale = ivue::color.scale.groups(unname(labels[names]),
      colors = stats::setNames(unname(palette[names]), unname(labels[names]))),
    point.type = 'point', point.size = if (graph$n > 500) 2.5 else 4,
    axes = FALSE, xlab = '', ylab = '', zlab = '', aspect = 'equal',
    camera = ivue::camera.zup(elevation = 22, turn = -125, zoom = .6),
    limits = bounds, layers = layers, legend.show = controls, controls = FALSE,
    width = width, height = height, legend.width = 180,
    description = if (controls) paste(graph$label,
      'Gray: generating coordinates (when available); blue: MDS; orange: MDS + edge-KK.') else NULL)
  if (!controls) return(widget)
  htmlwidgets::onRender(widget, "function(el, x, data) {
    var scene = el.rglinstance, root = scene.scene.rootSubscene;
    el.querySelectorAll('.ivue-legend details').forEach(function(node) { node.remove(); });
    var legend = el.querySelector('.ivue-legend');
    var rows = {};
    if (legend) Object.keys(data.labels).forEach(function(name) {
      rows[name] = Array.from(legend.children).find(function(node) {
        return node.tagName === 'DIV' && node.textContent.trim().startsWith(data.labels[name] + ' (');
      });
    });
    var panel = document.createElement('details');
    panel.className = 'ivue-tools'; panel.open = true;
    var summary = document.createElement('summary'); summary.textContent = 'View controls';
    panel.appendChild(summary);
    var label = document.createElement('label'); label.textContent = 'Layout: ';
    var select = document.createElement('select');
    select.className = 'grip-layout-selector'; select.style.font = 'inherit';
    [['mds', 'metric-MDS'], ['refined', 'metric-MDS + edge-KK'], ['both', 'both']].forEach(function(choice) {
      var option = document.createElement('option');
      option.value = choice[0]; option.textContent = choice[1];
      option.disabled = choice[0] === 'both' ? !(data.ids.mds && data.ids.refined) : !data.ids[choice[0]];
      select.appendChild(option);
    });
    select.value = data.initial;
    label.appendChild(select); panel.appendChild(label);
    var hint = document.createElement('span'); hint.textContent = ' Drag to rotate; scroll to zoom.';
    panel.appendChild(hint);
    function update() {
      Object.keys(data.ids).forEach(function(name) {
        var visible = name === 'reference' || select.value === 'both' || select.value === name;
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
    el.appendChild(panel); update();
  }", data = list(ids = object.ids, labels = as.list(labels[names]),
    initial = if (all(c('mds', 'refined') %in% names)) 'both' else
      if ('mds' %in% names) 'mds' else if ('refined' %in% names) 'refined' else 'both'))
}
