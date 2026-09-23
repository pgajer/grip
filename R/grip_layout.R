grip.build.adj.from.edges <- function(edges, n, edge.weights = NULL) {
  n <- grip.validate.vertex.count(n)
  edges <- as.matrix(edges)
  if (!is.numeric(edges) || ncol(edges) != 2) {
    stop("edges must be a two-column integer matrix of 1-based vertex ids")
  }
  grip.validate.vertex.ids(edges, "edges", n)
  edges <- matrix(as.integer(edges), ncol = 2)
  if (any(!is.finite(edges))) {
    stop("edges must contain finite integer vertex ids")
  }
  if (any(edges <= 0L | edges > n)) {
    stop("edges must be 1-based and within [1, n]")
  }

  use.weights <- !is.null(edge.weights)
  if (use.weights) {
    if (length(edge.weights) != nrow(edges)) {
      stop("edge.weights length must match number of edges")
    }
    if (!is.numeric(edge.weights)) {
      stop("edge.weights must be a numeric vector")
    }
    bad <- which(!is.finite(edge.weights) | edge.weights <= 0)
    if (length(bad) > 0L) {
      i <- bad[[1L]]
      stop(sprintf(
        "edge.weights must contain finite values > 0; first invalid at edge.weights[%d] = %s",
        i,
        format(edge.weights[i], digits = 16)
      ))
    }
    edge.weights <- as.double(edge.weights)
  }

  out.adj <- vector("list", n)
  out.w <- if (use.weights) vector("list", n) else NULL
  if (nrow(edges) == 0L) {
    return(list(adj_list = out.adj, weight_list = out.w))
  }

  for (i in seq_len(nrow(edges))) {
    u <- edges[i, 1L]
    v <- edges[i, 2L]
    if (u == v) next
    out.adj[[u]] <- c(out.adj[[u]], v)
    out.adj[[v]] <- c(out.adj[[v]], u)
    if (use.weights) {
      w <- edge.weights[[i]]
      out.w[[u]] <- c(out.w[[u]], w)
      out.w[[v]] <- c(out.w[[v]], w)
    }
  }

  list(adj_list = out.adj, weight_list = out.w)
}

grip.connected.components <- function(adj.list, n) {
  comp <- integer(n)
  cid <- 0L
  for (v in seq_len(n)) {
    if (comp[[v]] != 0L) next
    cid <- cid + 1L
    q <- integer(n)
    head <- 1L
    tail <- 1L
    q[[tail]] <- v
    comp[[v]] <- cid
    while (head <= tail) {
      x <- q[[head]]
      head <- head + 1L
      nb <- adj.list[[x]]
      if (length(nb) == 0L) next
      for (u in nb) {
        if (comp[[u]] == 0L) {
          tail <- tail + 1L
          q[[tail]] <- u
          comp[[u]] <- cid
        }
      }
    }
  }
  comp
}

grip.induce.subgraph <- function(adj.list, weight.list, vertices, n) {
  vertices <- as.integer(vertices)
  in.comp <- rep(FALSE, n)
  in.comp[vertices] <- TRUE
  map <- integer(n)
  map[vertices] <- seq_along(vertices)
  sub.adj <- vector("list", length(vertices))
  use.weights <- !is.null(weight.list)
  sub.w <- if (use.weights) vector("list", length(vertices)) else NULL

  for (i in seq_along(vertices)) {
    v <- vertices[[i]]
    nb <- as.integer(adj.list[[v]])
    if (length(nb) == 0L) {
      sub.adj[[i]] <- integer(0)
      if (use.weights) sub.w[[i]] <- numeric(0)
      next
    }
    keep <- in.comp[nb]
    nb.keep <- nb[keep]
    sub.adj[[i]] <- as.integer(map[nb.keep])
    if (use.weights) {
      sub.w[[i]] <- as.double(weight.list[[v]][keep])
    }
  }

  list(adj_list = sub.adj, weight_list = sub.w)
}

grip.pack.component.layouts <- function(layouts, comp, n, dim) {
  comp.ids <- sort(unique(comp))
  out <- matrix(NA_real_, nrow = n, ncol = dim)
  if (length(comp.ids) == 1L) {
    out[,] <- as.matrix(layouts[[1L]])
    return(out)
  }

  centered <- vector("list", length(layouts))
  radii <- rep(1.0, length(layouts))
  for (k in seq_along(layouts)) {
    z <- as.matrix(layouts[[k]])
    ctr <- colMeans(z)
    zc <- sweep(z, 2L, ctr, "-")
    centered[[k]] <- zc
    if (nrow(zc) > 0L) {
      rr <- sqrt(rowSums(zc * zc))
      rmax <- max(rr, na.rm = TRUE)
      if (is.finite(rmax) && rmax > 0) radii[[k]] <- as.double(rmax)
    }
  }

  ord <- order(vapply(comp.ids, function(id) sum(comp == id), integer(1L)), decreasing = TRUE)
  gap <- stats::median(radii, na.rm = TRUE)
  if (!is.finite(gap) || gap <= 0) gap <- 1.0
  gap <- gap * 1.5

  x.offset <- 0.0
  prev.r <- 0.0
  for (kk in seq_along(ord)) {
    idx <- ord[[kk]]
    cid <- comp.ids[[idx]]
    rows <- which(comp == cid)
    if (kk == 1L) {
      x.offset <- 0.0
    } else {
      x.offset <- x.offset + prev.r + radii[[idx]] + gap
    }
    shift <- rep(0.0, dim)
    shift[[1L]] <- x.offset
    out[rows, ] <- sweep(centered[[idx]], 2L, shift, "+")
    prev.r <- radii[[idx]]
  }
  out
}

grip.normalize.preset <- function(preset, fn = "grip") {
  if (is.null(preset)) {
    return(NULL)
  }
  if (!is.character(preset) || length(preset) != 1L || is.na(preset)) {
    stop("preset must be NULL, 'carpet', 'mesh', 'torus', or 'tree'")
  }
  if (preset %in% c("carpet", "mesh", "torus", "tree")) {
    return(preset)
  }
  stop("preset must be NULL, 'carpet', 'mesh', 'torus', or 'tree'")
}

grip.carpet.preset.defaults <- function() {
  list(
    placement = "barycenter",
    rounds = 160L,
    final_rounds = 288L,
    num_init = 28L,
    num_nbrs = 24L,
    r = 0.03,
    s = 6.0,
    repulsion_factor = 2.5
  )
}

grip.mesh.preset.defaults <- function() {
  list(
    placement = "barycenter",
    rounds = 128L,
    final_rounds = 128L,
    num_init = 12L,
    num_nbrs = 20L,
    r = 0.10,
    s = 4.5,
    repulsion_factor = 1.5
  )
}

grip.torus.preset.defaults <- function() {
  list(
    placement = "barycenter",
    rounds = 192L,
    final_rounds = 288L,
    num_init = 12L,
    num_nbrs = 28L,
    r = 0.05,
    s = 7.5,
    repulsion_factor = 0.75
  )
}

grip.tree.preset.defaults <- function(dim = 2L) {
  list(
    placement = if (isTRUE(as.integer(dim) == 2L)) "circle" else "barycenter",
    rounds = 64L,
    final_rounds = 160L,
    num_init = 28L,
    num_nbrs = 8L,
    r = 0.05,
    s = 7.5,
    repulsion_factor = 0.0
  )
}

grip.globalrep.default.final.rounds <- function(n) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n)) {
    stop("n must be a single finite numeric value")
  }
  n <- as.integer(n)
  if (is.na(n) || n <= 0L) {
    stop("n must be a positive integer")
  }
  if (n <= 1000L) {
    return(384L)
  }
  if (n <= 5000L) {
    return(320L)
  }
  if (n <= 20000L) {
    return(256L)
  }
  if (n <= 50000L) {
    return(200L)
  }
  128L
}

grip.globalrep.base.defaults <- function(n = NULL) {
  list(
    placement = "barycenter",
    rounds = 160L,
    final_rounds = if (is.null(n)) 384L else grip.globalrep.default.final.rounds(n),
    num_init = 24L,
    num_nbrs = 20L,
    r = 0.03,
    s = 7.5,
    repulsion_factor = 2.5,
    coarse_repulsion_factor = 1.5,
    coarse_repulsion_sample = 16L,
    coarse_repulsion_exact_below = 64L
  )
}

grip.forward.call <- function(target, call, env = parent.frame()) {
  args <- as.list(call)[-1L]
  do.call(target, args, envir = env)
}

grip.resolve.preset <- function(preset,
                                dim = 2L,
                                placement,
                                placement.missing,
                                rounds,
                                rounds.missing,
                                final.rounds,
                                final.rounds.missing,
                                num.init,
                                num.init.missing,
                                num.nbrs,
                                num.nbrs.missing,
                                r,
                                r.missing,
                                s,
                                s.missing,
                                repulsion.factor,
                                repulsion.factor.missing) {
  if (is.null(preset)) {
    return(list(
      placement = placement,
      rounds = rounds,
      final_rounds = final.rounds,
      num_init = num.init,
      num_nbrs = num.nbrs,
      r = r,
      s = s,
      repulsion_factor = repulsion.factor
    ))
  }

  defaults <- switch(
    preset,
    carpet = grip.carpet.preset.defaults(),
    mesh = grip.mesh.preset.defaults(),
    torus = grip.torus.preset.defaults(),
    tree = grip.tree.preset.defaults(dim = dim),
    stop("unknown preset")
  )
  if (placement.missing) placement <- defaults$placement
  if (rounds.missing) rounds <- defaults$rounds
  if (final.rounds.missing) final.rounds <- defaults$final_rounds
  if (num.init.missing) num.init <- defaults$num_init
  if (num.nbrs.missing) num.nbrs <- defaults$num_nbrs
  if (r.missing) r <- defaults$r
  if (s.missing) s <- defaults$s
  if (repulsion.factor.missing) repulsion.factor <- defaults$repulsion_factor

  list(
    placement = placement,
    rounds = rounds,
    final_rounds = final.rounds,
    num_init = num.init,
    num_nbrs = num.nbrs,
    r = r,
    s = s,
    repulsion_factor = repulsion.factor
  )
}

grip.validate.tuning.inputs <- function(num.nbrs, r, s, repulsion.factor) {
  if (!is.numeric(num.nbrs) || length(num.nbrs) != 1L || !is.finite(num.nbrs)) {
    stop("num.nbrs must be a single finite numeric value")
  }
  if (abs(num.nbrs - round(num.nbrs)) > sqrt(.Machine$double.eps)) {
    stop("num.nbrs must be a positive integer")
  }
  num.nbrs <- as.integer(round(num.nbrs))
  if (is.na(num.nbrs) || num.nbrs <= 0L) {
    stop("num.nbrs must be a positive integer")
  }

  if (!is.numeric(r) || length(r) != 1L || !is.finite(r)) {
    stop("r must be a single finite numeric value")
  }
  r <- as.double(r)
  if (r < 0 || r > 1) {
    stop("r must be in [0, 1]")
  }

  if (!is.numeric(s) || length(s) != 1L || !is.finite(s)) {
    stop("s must be a single finite numeric value")
  }
  s <- as.double(s)
  if (s < 0) {
    stop("s must be >= 0")
  }

  if (!is.numeric(repulsion.factor) || length(repulsion.factor) != 1L ||
      !is.finite(repulsion.factor)) {
    stop("repulsion.factor must be a single finite numeric value")
  }
  repulsion.factor <- as.double(repulsion.factor)
  if (repulsion.factor < 0) {
    stop("repulsion.factor must be >= 0")
  }

  list(
    num_nbrs = num.nbrs,
    r = r,
    s = s,
    repulsion_factor = repulsion.factor
  )
}

grip.validate.globalrep.tuning.inputs <- function(num.nbrs,
                                                  r,
                                                  s,
                                                  repulsion.factor,
                                                  coarse.repulsion.factor,
                                                  coarse.repulsion.sample,
                                                  coarse.repulsion.exact.below,
                                                  final.anchor.factor = 0,
                                                  final.move.scale.after.first = 1,
                                                  insertion.anchor.count = 3,
                                                  insertion.anchor.scope = "any_higher",
                                                  insertion.anchor.strategy = "first",
                                                  level0.insertion.mode = "inherit",
                                                  level0.anchor.count = insertion.anchor.count,
                                                  level0.local.kk.steps = 3) {
  tuning <- grip.validate.tuning.inputs(
    num.nbrs = num.nbrs,
    r = r,
    s = s,
    repulsion.factor = repulsion.factor
  )

  if (!is.numeric(coarse.repulsion.factor) ||
      length(coarse.repulsion.factor) != 1L ||
      !is.finite(coarse.repulsion.factor)) {
    stop("coarse.repulsion.factor must be a single finite numeric value")
  }
  coarse.repulsion.factor <- as.double(coarse.repulsion.factor)
  if (coarse.repulsion.factor < 0) {
    stop("coarse.repulsion.factor must be >= 0")
  }

  if (!is.numeric(coarse.repulsion.sample) ||
      length(coarse.repulsion.sample) != 1L ||
      !is.finite(coarse.repulsion.sample)) {
    stop("coarse.repulsion.sample must be a single finite numeric value")
  }
  if (abs(coarse.repulsion.sample - round(coarse.repulsion.sample)) >
      sqrt(.Machine$double.eps)) {
    stop("coarse.repulsion.sample must be a positive integer")
  }
  coarse.repulsion.sample <- as.integer(round(coarse.repulsion.sample))
  if (is.na(coarse.repulsion.sample) || coarse.repulsion.sample <= 0L) {
    stop("coarse.repulsion.sample must be a positive integer")
  }

  if (!is.numeric(coarse.repulsion.exact.below) ||
      length(coarse.repulsion.exact.below) != 1L ||
      !is.finite(coarse.repulsion.exact.below)) {
    stop("coarse.repulsion.exact.below must be a single finite numeric value")
  }
  if (abs(coarse.repulsion.exact.below - round(coarse.repulsion.exact.below)) >
      sqrt(.Machine$double.eps)) {
    stop("coarse.repulsion.exact.below must be a positive integer")
  }
  coarse.repulsion.exact.below <- as.integer(round(coarse.repulsion.exact.below))
  if (is.na(coarse.repulsion.exact.below) || coarse.repulsion.exact.below <= 0L) {
    stop("coarse.repulsion.exact.below must be a positive integer")
  }

  if (!is.numeric(final.anchor.factor) ||
      length(final.anchor.factor) != 1L ||
      !is.finite(final.anchor.factor)) {
    stop("final.anchor.factor must be a single finite numeric value")
  }
  final.anchor.factor <- as.double(final.anchor.factor)
  if (final.anchor.factor < 0) {
    stop("final.anchor.factor must be >= 0")
  }

  if (!is.numeric(final.move.scale.after.first) ||
      length(final.move.scale.after.first) != 1L ||
      !is.finite(final.move.scale.after.first)) {
    stop("final.move.scale.after.first must be a single finite numeric value")
  }
  final.move.scale.after.first <- as.double(final.move.scale.after.first)
  if (final.move.scale.after.first < 0 || final.move.scale.after.first > 1) {
    stop("final.move.scale.after.first must be in [0, 1]")
  }

  if (!is.numeric(insertion.anchor.count) ||
      length(insertion.anchor.count) != 1L ||
      !is.finite(insertion.anchor.count)) {
    stop("insertion.anchor.count must be a single finite numeric value")
  }
  if (abs(insertion.anchor.count - round(insertion.anchor.count)) >
      sqrt(.Machine$double.eps)) {
    stop("insertion.anchor.count must be a positive integer")
  }
  insertion.anchor.count <- as.integer(round(insertion.anchor.count))
  if (is.na(insertion.anchor.count) || insertion.anchor.count <= 0L) {
    stop("insertion.anchor.count must be a positive integer")
  }

  insertion.anchor.scope <- match.arg(
    insertion.anchor.scope,
    choices = c("any_higher", "prev_misf")
  )

  insertion.anchor.strategy <- match.arg(
    insertion.anchor.strategy,
    choices = c("first", "distance_band", "balanced_band", "spread_prev")
  )

  level0.insertion.mode <- match.arg(
    level0.insertion.mode,
    choices = c("inherit", "barycenter", "least_squares")
  )

  if (!is.numeric(level0.anchor.count) ||
      length(level0.anchor.count) != 1L ||
      !is.finite(level0.anchor.count)) {
    stop("level0.anchor.count must be a single finite numeric value")
  }
  if (abs(level0.anchor.count - round(level0.anchor.count)) >
      sqrt(.Machine$double.eps)) {
    stop("level0.anchor.count must be a positive integer")
  }
  level0.anchor.count <- as.integer(round(level0.anchor.count))
  if (is.na(level0.anchor.count) || level0.anchor.count <= 0L) {
    stop("level0.anchor.count must be a positive integer")
  }

  if (!is.numeric(level0.local.kk.steps) ||
      length(level0.local.kk.steps) != 1L ||
      !is.finite(level0.local.kk.steps)) {
    stop("level0.local.kk.steps must be a single finite numeric value")
  }
  if (abs(level0.local.kk.steps - round(level0.local.kk.steps)) >
      sqrt(.Machine$double.eps)) {
    stop("level0.local.kk.steps must be a non-negative integer")
  }
  level0.local.kk.steps <- as.integer(round(level0.local.kk.steps))
  if (is.na(level0.local.kk.steps) || level0.local.kk.steps < 0L) {
    stop("level0.local.kk.steps must be a non-negative integer")
  }

  c(
    tuning,
    list(
      coarse_repulsion_factor = coarse.repulsion.factor,
      coarse_repulsion_sample = coarse.repulsion.sample,
      coarse_repulsion_exact_below = coarse.repulsion.exact.below,
      final_anchor_factor = final.anchor.factor,
      final_move_scale_after_first = final.move.scale.after.first,
      insertion_anchor_count = insertion.anchor.count,
      insertion_anchor_scope = insertion.anchor.scope,
      insertion_anchor_strategy = insertion.anchor.strategy,
      level0_insertion_mode = level0.insertion.mode,
      level0_anchor_count = level0.anchor.count,
      level0_local_kk_steps = level0.local.kk.steps
    )
  )
}

grip.validate.lgkk.polish.inputs <- function(lgkk.polish.rounds = 0L,
                                             lgkk.multiscale.rounds = 0L,
                                             lgkk.rounds.coarse = NULL,
                                             lgkk.rounds.pre.final = NULL,
                                             lgkk.rounds.final = NULL,
                                             lgkk.local.nbrs = 20L,
                                             lgkk.landmark.count = 8L,
                                             lgkk.multiscale.scope = "all",
                                             lgkk.active.limit = 4096L) {
  if (!is.numeric(lgkk.polish.rounds) ||
      length(lgkk.polish.rounds) != 1L ||
      !is.finite(lgkk.polish.rounds)) {
    stop("lgkk.polish.rounds must be a single finite numeric value")
  }
  if (abs(lgkk.polish.rounds - round(lgkk.polish.rounds)) >
      sqrt(.Machine$double.eps)) {
    stop("lgkk.polish.rounds must be a non-negative integer")
  }
  lgkk.polish.rounds <- as.integer(round(lgkk.polish.rounds))
  if (is.na(lgkk.polish.rounds) || lgkk.polish.rounds < 0L) {
    stop("lgkk.polish.rounds must be a non-negative integer")
  }

  if (!is.numeric(lgkk.multiscale.rounds) ||
      length(lgkk.multiscale.rounds) != 1L ||
      !is.finite(lgkk.multiscale.rounds)) {
    stop("lgkk.multiscale.rounds must be a single finite numeric value")
  }
  if (abs(lgkk.multiscale.rounds - round(lgkk.multiscale.rounds)) >
      sqrt(.Machine$double.eps)) {
    stop("lgkk.multiscale.rounds must be a non-negative integer")
  }
  lgkk.multiscale.rounds <- as.integer(round(lgkk.multiscale.rounds))
  if (is.na(lgkk.multiscale.rounds) || lgkk.multiscale.rounds < 0L) {
    stop("lgkk.multiscale.rounds must be a non-negative integer")
  }

  validate.optional.rounds <- function(x, name) {
    if (is.null(x)) {
      return(NULL)
    }
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
      stop(sprintf("%s must be NULL or a single finite numeric value", name))
    }
    if (abs(x - round(x)) > sqrt(.Machine$double.eps)) {
      stop(sprintf("%s must be a non-negative integer", name))
    }
    x <- as.integer(round(x))
    if (is.na(x) || x < 0L) {
      stop(sprintf("%s must be a non-negative integer", name))
    }
    x
  }

  lgkk.rounds.coarse <- validate.optional.rounds(lgkk.rounds.coarse, "lgkk.rounds.coarse")
  lgkk.rounds.pre.final <- validate.optional.rounds(lgkk.rounds.pre.final, "lgkk.rounds.pre.final")
  lgkk.rounds.final <- validate.optional.rounds(lgkk.rounds.final, "lgkk.rounds.final")

  if (is.null(lgkk.rounds.coarse)) {
    lgkk.rounds.coarse <- lgkk.multiscale.rounds
  }
  if (is.null(lgkk.rounds.pre.final)) {
    lgkk.rounds.pre.final <- lgkk.multiscale.rounds
  }
  if (is.null(lgkk.rounds.final)) {
    lgkk.rounds.final <- lgkk.multiscale.rounds
  }

  if (!is.numeric(lgkk.local.nbrs) ||
      length(lgkk.local.nbrs) != 1L ||
      !is.finite(lgkk.local.nbrs)) {
    stop("lgkk.local.nbrs must be a single finite numeric value")
  }
  if (abs(lgkk.local.nbrs - round(lgkk.local.nbrs)) >
      sqrt(.Machine$double.eps)) {
    stop("lgkk.local.nbrs must be a non-negative integer")
  }
  lgkk.local.nbrs <- as.integer(round(lgkk.local.nbrs))
  if (is.na(lgkk.local.nbrs) || lgkk.local.nbrs < 0L) {
    stop("lgkk.local.nbrs must be a non-negative integer")
  }

  if (!is.numeric(lgkk.landmark.count) ||
      length(lgkk.landmark.count) != 1L ||
      !is.finite(lgkk.landmark.count)) {
    stop("lgkk.landmark.count must be a single finite numeric value")
  }
  if (abs(lgkk.landmark.count - round(lgkk.landmark.count)) >
      sqrt(.Machine$double.eps)) {
    stop("lgkk.landmark.count must be a non-negative integer")
  }
  lgkk.landmark.count <- as.integer(round(lgkk.landmark.count))
  if (is.na(lgkk.landmark.count) || lgkk.landmark.count < 0L) {
    stop("lgkk.landmark.count must be a non-negative integer")
  }

  lgkk.multiscale.scope <- match.arg(
    lgkk.multiscale.scope,
    choices = c("all", "coarse")
  )

  if (!is.numeric(lgkk.active.limit) ||
      length(lgkk.active.limit) != 1L ||
      !is.finite(lgkk.active.limit)) {
    stop("lgkk.active.limit must be a single finite numeric value")
  }
  if (abs(lgkk.active.limit - round(lgkk.active.limit)) >
      sqrt(.Machine$double.eps)) {
    stop("lgkk.active.limit must be a positive integer")
  }
  lgkk.active.limit <- as.integer(round(lgkk.active.limit))
  if (is.na(lgkk.active.limit) || lgkk.active.limit <= 0L) {
    stop("lgkk.active.limit must be a positive integer")
  }

  list(
    lgkk_polish_rounds = lgkk.polish.rounds,
    lgkk_multiscale_rounds = lgkk.multiscale.rounds,
    lgkk_rounds_coarse = lgkk.rounds.coarse,
    lgkk_rounds_pre_final = lgkk.rounds.pre.final,
    lgkk_rounds_final = lgkk.rounds.final,
    lgkk_local_nbrs = lgkk.local.nbrs,
    lgkk_landmark_count = lgkk.landmark.count,
    lgkk_multiscale_scope = lgkk.multiscale.scope,
    lgkk_active_limit = lgkk.active.limit
  )
}

grip.apply.lgkk.polish <- function(coords,
                                   adj.list,
                                   weight.list,
                                   rounds,
                                   lgkk.local.nbrs,
                                   lgkk.landmark.count,
                                   return.trace = FALSE) {
  if (is.null(rounds) || rounds <= 0L) {
    return(list(
      coords = coords,
      trace = data.frame(),
      frames = list(coords)
    ))
  }
  prepared <- prepare.landmark.geodesic.kk(
    adj.list = adj.list,
    weight.list = weight.list,
    n = nrow(coords),
    local.nbrs = lgkk.local.nbrs,
    landmark.count = lgkk.landmark.count
  )
  landmark.geodesic.kk(
    coords = coords,
    prepared = prepared,
    max.iter = rounds,
    local.nbrs = lgkk.local.nbrs,
    landmark.count = lgkk.landmark.count,
    return.trace = return.trace
  )
}

#' Compute a GRIP layout with coarse global repulsion
#'
#' This compatibility entry point is equivalent to
#' \code{\link{grip}(..., metric = "hop")}. It keeps the old
#' global-repulsion name available while using the same quality-first
#' multiscale GRIP engine and adaptive \code{final.rounds} schedule as the
#' primary layout API.
#'
#' @inheritParams grip
#' @param weight.list Optional parallel list of positive edge lengths for
#'   \code{adj.list}. \code{NULL} treats every edge as length 1.
#' @param edge.weights Optional positive edge lengths for \code{edges}, in row
#'   order. \code{NULL} treats every edge as length 1.
#' @param coarse.repulsion.factor Non-negative multiplier applied to the extra
#'   coarse-level active-set repulsion term. \code{0} disables that extra term;
#'   when the remaining tuning arguments also match
#'   \code{\link{legacy.grip}()}, the result matches the legacy
#'   layout behavior.
#' @param coarse.repulsion.sample Positive integer sample size used to
#'   approximate active-set-wide repulsion on larger coarse levels.
#' @param coarse.repulsion.exact.below Positive integer threshold. When the
#'   active set size is at most this value, the coarse repulsion is computed
#'   exactly against all currently active vertices instead of being sampled.
#' @param final.anchor.factor Non-negative multiplier for an anchor term that
#'   pulls the final FR stage back toward the pre-final full-graph layout.
#'   `0` disables the anchor and preserves the current behavior.
#' @param final.move.scale.after.first Scalar in `[0, 1]` applied to the final
#'   FR displacement after the first finest-level round. Values below `1`
#'   damp later full-graph movement while keeping the first FR round unchanged.
#' @param final.mode Final full-graph refinement mode. \code{"fr"} keeps the
#'   current Fruchterman-Reingold-style final stage. \code{"kk_repulse"} uses a
#'   KK-style local distance-matching update with explicit active-set
#'   repulsion instead of the final FR phase.
#' @param insertion.anchor.count Positive integer number of anchor vertices used
#'   during multiscale insertion on non-initial MISF refinement levels. This is
#'   the closest current implementation to a global \code{K_mish} parameter.
#' @param insertion.anchor.scope Anchor-eligibility rule used during multiscale
#'   insertion. \code{"any_higher"} matches the historical GRIP behavior and
#'   allows anchors from any already placed higher MISF level.
#'   \code{"prev_misf"} restricts anchors to the immediately previous MISF
#'   level only.
#' @param insertion.anchor.strategy Anchor-selection rule used during
#'   multiscale insertion. \code{"first"} keeps the historical
#'   first-anchors-found BFS behavior. \code{"distance_band"} keeps exploring
#'   until the \code{K_mish}-th anchor distance band is exhausted, then places
#'   the new vertex from that less order-sensitive anchor pool.
#'   \code{"balanced_band"} uses the same band expansion, then explicitly
#'   selects a subset whose centroid stays centered in the candidate cloud while
#'   remaining geometrically spread out. \code{"spread_prev"} is a
#'   symmetry-oriented band strategy intended to be paired with
#'   \code{insertion.anchor.scope = "prev_misf"}; it selects anchors with broad
#'   angular and geometric coverage before placement.
#' @param level0.insertion.mode Level-0 insertion placement override used only
#'   when the finest filtration level is first populated. \code{"inherit"}
#'   keeps the current GRIP behavior. \code{"barycenter"} disables the 2D
#'   circle heuristic at level 0 and uses barycentric anchor placement.
#'   \code{"least_squares"} uses a multi-anchor least-squares distance fit at
#'   level 0 before any local micro-polish.
#' @param level0.anchor.count Positive integer number of already placed anchors
#'   to collect for level-0 insertion experiments. By default this inherits
#'   \code{insertion.anchor.count}. The legacy behavior uses 3.
#' @param level0.local.kk.steps Non-negative integer number of tiny local KK
#'   micro-polish steps applied immediately after each level-0 insertion. The
#'   legacy behavior uses 3.
#' @param lgkk.polish.rounds Non-negative integer number of experimental
#'   landmark-geodesic KK polish iterations applied after the main GRIP solve.
#'   \code{0} disables the polish.
#' @param lgkk.multiscale.rounds Non-negative integer number of compiled
#'   landmark-geodesic KK refinement rounds applied inside the multiscale solver
#'   after each eligible MISF level completes its standard GRIP rounds.
#'   This legacy shared budget is used as a fallback when any of the more
#'   specific per-stage budgets below are left \code{NULL}.
#' @param lgkk.rounds.coarse Optional non-negative integer number of compiled
#'   LGKK rounds applied on coarse MISF levels with \code{misf_level > 1}.
#'   When \code{NULL}, this falls back to \code{lgkk.multiscale.rounds}.
#' @param lgkk.rounds.pre.final Optional non-negative integer number of compiled
#'   LGKK rounds applied on the last coarse level just before the full graph is
#'   opened (\code{misf_level == 1}). When \code{NULL}, this falls back to
#'   \code{lgkk.multiscale.rounds}.
#' @param lgkk.rounds.final Optional non-negative integer number of compiled
#'   LGKK rounds applied after the full graph level completes its standard GRIP
#'   rounds (\code{misf_level == 0}). When \code{NULL}, this falls back to
#'   \code{lgkk.multiscale.rounds}.
#' @param lgkk.local.nbrs Number of nearest graph-metric neighbors retained per
#'   vertex in the LGKK sparse local set when either LGKK stage is enabled.
#' @param lgkk.landmark.count Number of farthest-point landmarks retained per
#'   vertex in the LGKK sparse long-range set when either LGKK stage is
#'   enabled.
#' @param lgkk.multiscale.scope Scope for the compiled multiscale LGKK stage.
#'   \code{"all"} applies it after every eligible MISF level, including the
#'   final full-graph level. \code{"coarse"} applies it only on coarse levels.
#' @param lgkk.active.limit Positive integer upper bound on the active-set size
#'   for compiled multiscale LGKK cache construction. Levels larger than this
#'   skip the multiscale LGKK stage.
#' @return A numeric matrix with `n` rows and `dim` columns.
#' @examples
#' edges <- edges.mesh(4, 4)
#' coords <- globalrep.grip(edges, n = max(edges), dim = 2,
#'                                 rounds = 8, final.rounds = 8,
#'                                 num.init = 6, num.nbrs = 8,
#'                                 coarse.repulsion.factor = 0.2,
#'                                 coarse.repulsion.sample = 8,
#'                                 coarse.repulsion.exact.below = 32,
#'                                 seed = 1)
#' round(coords, 3)
#' @export
#' @md
globalrep.grip <- function(edges = NULL,
                                  n = NULL,
                                  adj.list = NULL,
                                  weight.list = NULL,
                                  edge.weights = NULL,
                                  dim = 3,
                                  placement = c("barycenter", "circle"),
                                  preset = NULL,
                                  rounds = 160,
                                  final.rounds = 384,
                                  num.init = 24,
                                  num.nbrs = 20,
                                  r = 0.03,
                                  s = 7.5,
                                  repulsion.factor = 2.5,
                                  coarse.repulsion.factor = 1.5,
                                  coarse.repulsion.sample = 16,
                                  coarse.repulsion.exact.below = 64,
                                  final.anchor.factor = 0,
                                  final.move.scale.after.first = 1,
                                  final.mode = c("fr", "kk_repulse"),
                                  insertion.anchor.count = 3,
                                  insertion.anchor.scope = c("any_higher", "prev_misf"),
                                  insertion.anchor.strategy = c("first", "distance_band", "balanced_band", "spread_prev"),
                                  level0.insertion.mode = c("inherit", "barycenter", "least_squares"),
                                  level0.anchor.count = insertion.anchor.count,
                                  level0.local.kk.steps = 3,
                                  lgkk.polish.rounds = 0L,
                                  lgkk.multiscale.rounds = 0L,
                                  lgkk.rounds.coarse = NULL,
                                  lgkk.rounds.pre.final = NULL,
                                  lgkk.rounds.final = NULL,
                                  lgkk.local.nbrs = 20L,
                                  lgkk.landmark.count = 8L,
                                  lgkk.multiscale.scope = c("all", "coarse"),
                                  lgkk.active.limit = 4096L,
                                  tinit.factor = 6,
                                  seed = 6,
                                  disconnected = c("components", "error")) {
  grip.validate.graph.arguments(edges, n, adj.list, weight.list, edge.weights)
  placement_missing <- missing(placement)
  rounds_missing <- missing(rounds)
  final_rounds_missing <- missing(final.rounds)
  num_init_missing <- missing(num.init)
  num_nbrs_missing <- missing(num.nbrs)
  r_missing <- missing(r)
  s_missing <- missing(s)
  repulsion_factor_missing <- missing(repulsion.factor)

  preset <- grip.normalize.preset(preset, fn = "globalrep.grip")
  resolved <- grip.resolve.preset(
    preset = preset,
    dim = dim,
    placement = placement,
    placement.missing = placement_missing,
    rounds = rounds,
    rounds.missing = rounds_missing,
    final.rounds = final.rounds,
    final.rounds.missing = final_rounds_missing,
    num.init = num.init,
    num.init.missing = num_init_missing,
    num.nbrs = num.nbrs,
    num.nbrs.missing = num_nbrs_missing,
    r = r,
    r.missing = r_missing,
    s = s,
    s.missing = s_missing,
    repulsion.factor = repulsion.factor,
    repulsion.factor.missing = repulsion_factor_missing
  )
  placement <- resolved$placement
  rounds <- resolved$rounds
  final.rounds <- resolved$final_rounds
  num.init <- resolved$num_init
  num.nbrs <- resolved$num_nbrs
  r <- resolved$r
  s <- resolved$s
  repulsion.factor <- resolved$repulsion_factor
  placement <- match.arg(placement)
  final.mode <- match.arg(final.mode)
  insertion.anchor.scope <- match.arg(insertion.anchor.scope)
  insertion.anchor.strategy <- match.arg(insertion.anchor.strategy)
  level0.insertion.mode <- match.arg(level0.insertion.mode)
  lgkk.multiscale.scope <- match.arg(lgkk.multiscale.scope)
  disconnected <- match.arg(disconnected)

  validated <- grip.validate.layout.inputs(
    edges = edges,
    n = n,
    adj.list = adj.list,
    weight.list = weight.list,
    edge.weights = edge.weights,
    dim = dim,
    placement = placement,
    seed = seed
  )
  adj.list <- validated$adj_list
  weight.list <- validated$weight_list
  n <- validated$n
  dim <- validated$dim
  seed <- validated$seed
  if (is.null(preset) && final_rounds_missing) {
    final.rounds <- grip.globalrep.default.final.rounds(n)
  }
  tuning <- grip.validate.globalrep.tuning.inputs(
    num.nbrs = num.nbrs,
    r = r,
    s = s,
    repulsion.factor = repulsion.factor,
    coarse.repulsion.factor = coarse.repulsion.factor,
    coarse.repulsion.sample = coarse.repulsion.sample,
    coarse.repulsion.exact.below = coarse.repulsion.exact.below,
    final.anchor.factor = final.anchor.factor,
    final.move.scale.after.first = final.move.scale.after.first,
    insertion.anchor.count = insertion.anchor.count,
    insertion.anchor.scope = insertion.anchor.scope,
    insertion.anchor.strategy = insertion.anchor.strategy,
    level0.insertion.mode = level0.insertion.mode,
    level0.anchor.count = level0.anchor.count,
    level0.local.kk.steps = level0.local.kk.steps
  )
  num.nbrs <- tuning$num_nbrs
  r <- tuning$r
  s <- tuning$s
  repulsion.factor <- tuning$repulsion_factor
  coarse.repulsion.factor <- tuning$coarse_repulsion_factor
  coarse.repulsion.sample <- tuning$coarse_repulsion_sample
  coarse.repulsion.exact.below <- tuning$coarse_repulsion_exact_below
  final.anchor.factor <- tuning$final_anchor_factor
  final.move.scale.after.first <- tuning$final_move_scale_after_first
  insertion.anchor.count <- tuning$insertion_anchor_count
  insertion.anchor.scope <- tuning$insertion_anchor_scope
  insertion.anchor.strategy <- tuning$insertion_anchor_strategy
  level0.insertion.mode <- tuning$level0_insertion_mode
  level0.anchor.count <- tuning$level0_anchor_count
  level0.local.kk.steps <- tuning$level0_local_kk_steps
  lgkk <- grip.validate.lgkk.polish.inputs(
    lgkk.polish.rounds = lgkk.polish.rounds,
    lgkk.multiscale.rounds = lgkk.multiscale.rounds,
    lgkk.rounds.coarse = lgkk.rounds.coarse,
    lgkk.rounds.pre.final = lgkk.rounds.pre.final,
    lgkk.rounds.final = lgkk.rounds.final,
    lgkk.local.nbrs = lgkk.local.nbrs,
    lgkk.landmark.count = lgkk.landmark.count,
    lgkk.multiscale.scope = lgkk.multiscale.scope,
    lgkk.active.limit = lgkk.active.limit
  )
  lgkk.polish.rounds <- lgkk$lgkk_polish_rounds
  lgkk.multiscale.rounds <- lgkk$lgkk_multiscale_rounds
  lgkk.rounds.coarse <- lgkk$lgkk_rounds_coarse
  lgkk.rounds.pre.final <- lgkk$lgkk_rounds_pre_final
  lgkk.rounds.final <- lgkk$lgkk_rounds_final
  lgkk.local.nbrs <- lgkk$lgkk_local_nbrs
  lgkk.landmark.count <- lgkk$lgkk_landmark_count
  lgkk.multiscale.scope <- lgkk$lgkk_multiscale_scope
  lgkk.active.limit <- lgkk$lgkk_active_limit

  layout.adj <- function(adj.list, weight.list, n) {
    coords <- grip_layout_globalrep_adj_cpp(
      adj_list = adj.list,
      weight_list = weight.list,
      n = n,
      dim = dim,
      placement = placement,
      rounds = as.integer(rounds),
      final_rounds = as.integer(final.rounds),
      num_init = as.integer(num.init),
      num_nbrs = num.nbrs,
      r = r,
      s = s,
      repulsion_factor = repulsion.factor,
      coarse_repulsion_factor = coarse.repulsion.factor,
      coarse_repulsion_sample = coarse.repulsion.sample,
      coarse_repulsion_exact_below = coarse.repulsion.exact.below,
      final_anchor_factor = final.anchor.factor,
      final_move_scale_after_first = final.move.scale.after.first,
      insertion_anchor_count = insertion.anchor.count,
      insertion_anchor_scope = insertion.anchor.scope,
      insertion_anchor_strategy = insertion.anchor.strategy,
      level0_insertion_mode = level0.insertion.mode,
      level0_anchor_count = level0.anchor.count,
      level0_local_kk_steps = level0.local.kk.steps,
      lgkk_multiscale_rounds = lgkk.multiscale.rounds,
      lgkk_rounds_coarse = lgkk.rounds.coarse,
      lgkk_rounds_pre_final = lgkk.rounds.pre.final,
      lgkk_rounds_final = lgkk.rounds.final,
      lgkk_local_nbrs = lgkk.local.nbrs,
      lgkk_landmark_count = lgkk.landmark.count,
      lgkk_multiscale_scope = lgkk.multiscale.scope,
      lgkk_active_limit = lgkk.active.limit,
      final_mode = final.mode,
      tinit_factor = as.integer(tinit.factor),
      seed = seed
    )
    polished <- grip.apply.lgkk.polish(
      coords = coords,
      adj.list = adj.list,
      weight.list = weight.list,
      rounds = lgkk.polish.rounds,
      lgkk.local.nbrs = lgkk.local.nbrs,
      lgkk.landmark.count = lgkk.landmark.count,
      return.trace = FALSE
    )
    polished$coords
  }

  comp <- grip.connected.components(adj.list = adj.list, n = n)
  n.comp <- length(unique(comp))

  if (n.comp == 1L) {
    return(layout.adj(adj.list = adj.list, weight.list = weight.list, n = n))
  }

  if (identical(disconnected, "error")) {
    stop(sprintf(
      "Input graph has %d connected components; the GRIP layout core assumes connected graphs. Use disconnected = 'components' to lay out each component safely.",
      n.comp
    ))
  }

  warning(
    sprintf(
      "Input graph has %d connected components; laying out components separately to avoid disconnected-graph instability.",
      n.comp
    ),
    call. = FALSE
  )

  comp.ids <- sort(unique(comp))
  layouts <- vector("list", length(comp.ids))
  for (k in seq_along(comp.ids)) {
    rows <- which(comp == comp.ids[[k]])
    sub <- grip.induce.subgraph(
      adj.list = adj.list,
      weight.list = weight.list,
      vertices = rows,
      n = n
    )
    layouts[[k]] <- layout.adj(
      adj.list = sub$adj_list,
      weight.list = sub$weight_list,
      n = length(rows)
    )
  }
  grip.pack.component.layouts(layouts = layouts, comp = comp, n = n, dim = dim)
}

#' Compute a GRIP layout
#'
#' This is the primary layout API. It uses the quality-first multiscale GRIP
#' engine with extra coarse-level global repulsion to reduce foldovers while
#' preserving the usual GRIP refinement structure. With \code{preset = NULL},
#' the default profile is tuned for higher-quality layouts and automatically
#' tapers \code{final.rounds} on larger graphs.
#'
#' @param edges Two-column integer matrix of edges (1-based vertex ids).
#'   Supply either \code{edges}/\code{edge.weights} or
#'   \code{adj.list}/\code{weight.list}, not both. Fractional, nonfinite,
#'   and out-of-range vertex ids are rejected before integer conversion.
#' @param n Number of vertices, a finite positive integer.
#' @param adj.list Adjacency list (1-based) for undirected graphs.
#' @param weight.list Parallel list of edge lengths for \code{adj.list}.
#'   The edge-length-metric engine requires it; in the hop-metric engine,
#'   \code{NULL} treats all edges as length 1. All supplied lengths must be
#'   finite and strictly positive.
#' @param edge.weights Vector of edge lengths for \code{edges}, in the same
#'   order as its rows. The edge-length-metric engine requires it; in the
#'   hop-metric engine, \code{NULL} treats all edges as length 1. All supplied
#'   lengths must be finite and strictly positive. The selected engine
#'   determines whether lengths also define graph distances.
#' @param metric Graph metric used by the multiscale GRIP engine.
#'   \code{"hop"} (default) uses unweighted shortest-path hop counts to build
#'   the MISF hierarchy and graph neighborhoods. \code{"edge_length"} uses
#'   shortest-path distances obtained by summing the supplied positive edge
#'   lengths. See \bold{Edge-length semantics} for the exact role of edge
#'   lengths in each mode.
#' @param dim Layout dimension (2 or 3). Default is 3.
#' @param placement Initial placement strategy. "circle" is only used for 2D.
#' @param preset Optional tuning preset. \code{NULL} uses the quality-first
#'   defaults. \code{"carpet"} applies a preset tuned for
#'   Sierpinski-carpet-like graphs and validated on carpet levels 3 and 4.
#'   \code{"mesh"} applies a preset tuned for rectangular lattice graphs and
#'   validated on 8x8 and 12x12 mesh layouts. \code{"torus"} applies a preset
#'   tuned for 3D torus layouts and validated on torus sizes from 8x8 through
#'   20x20. \code{"tree"} applies a preset tuned for symmetric force-directed
#'   layouts of tree-like graphs and validated on binary trees of depths 5 and
#'   6. Presets only fill in tuning arguments that you did not supply
#'   explicitly.
#' @param rounds Initial rounds for refinement.
#' @param final.rounds Final rounds for refinement.
#' @param num.init Number of initial vertices in the coarsest level.
#' @param num.nbrs Maximum number of graph-distance neighbors retained for local
#'   refinement at each filtration level.
#' @param r Main local temperature adaptation rate in \code{[0, 1]}.
#' @param s Non-negative boost factor applied when successive displacements have
#'   a consistent direction.
#' @param repulsion.factor Non-negative multiplier applied to GRIP's
#'   finest-level repulsive force scale.
#' @param coarse.repulsion.factor Non-negative multiplier applied to the extra
#'   coarse-level active-set repulsion term. \code{0} disables that extra term.
#' @param coarse.repulsion.sample Positive integer sample size used to
#'   approximate active-set-wide repulsion on larger coarse levels.
#' @param coarse.repulsion.exact.below Positive integer threshold. When the
#'   active set size is at most this value, the coarse repulsion is computed
#'   exactly against all currently active vertices instead of being sampled.
#' @param final.anchor.factor Non-negative multiplier for an anchor term that
#'   pulls the final FR stage back toward the pre-final full-graph layout.
#'   `0` disables the anchor and preserves the current behavior.
#' @param final.move.scale.after.first Scalar in `[0, 1]` applied to the final
#'   FR displacement after the first finest-level round. Values below `1`
#'   damp later full-graph movement while keeping the first FR round unchanged.
#' @param final.mode Final full-graph refinement mode. \code{"fr"} keeps the
#'   current Fruchterman-Reingold-style final stage. \code{"kk_repulse"} uses a
#'   KK-style local distance-matching update with explicit active-set
#'   repulsion instead of the final FR phase.
#' @param insertion.anchor.count Positive integer number of anchor vertices used
#'   during multiscale insertion on non-initial MISF refinement levels. This is
#'   the closest current implementation to a global \code{K_mish} parameter.
#' @param insertion.anchor.scope Anchor-eligibility rule used during multiscale
#'   insertion. \code{"any_higher"} matches the historical GRIP behavior and
#'   allows anchors from any already placed higher MISF level.
#'   \code{"prev_misf"} restricts anchors to the immediately previous MISF
#'   level only.
#' @param insertion.anchor.strategy Anchor-selection rule used during
#'   multiscale insertion. \code{"first"} keeps the historical
#'   first-anchors-found BFS behavior. \code{"distance_band"} keeps exploring
#'   until the \code{K_mish}-th anchor distance band is exhausted, then places
#'   the new vertex from that less order-sensitive anchor pool.
#'   \code{"balanced_band"} uses the same band expansion, then explicitly
#'   selects a subset whose centroid stays centered in the candidate cloud while
#'   remaining geometrically spread out. \code{"spread_prev"} is a
#'   symmetry-oriented band strategy intended to be paired with
#'   \code{insertion.anchor.scope = "prev_misf"}; it selects anchors with broad
#'   angular and geometric coverage before placement.
#' @param level0.insertion.mode Level-0 insertion placement override used only
#'   when the finest filtration level is first populated. \code{"inherit"}
#'   keeps the current GRIP behavior. \code{"barycenter"} disables the 2D
#'   circle heuristic at level 0 and uses barycentric anchor placement.
#'   \code{"least_squares"} uses a multi-anchor least-squares distance fit at
#'   level 0 before any local micro-polish.
#' @param level0.anchor.count Positive integer number of already placed anchors
#'   to collect for level-0 insertion experiments. By default this inherits
#'   \code{insertion.anchor.count}. The legacy behavior uses 3.
#' @param level0.local.kk.steps Non-negative integer number of tiny local KK
#'   micro-polish steps applied immediately after each level-0 insertion. The
#'   legacy behavior uses 3.
#' @param lgkk.polish.rounds Non-negative integer number of experimental
#'   landmark-geodesic KK polish iterations applied after the main GRIP solve.
#'   \code{0} disables the polish.
#' @param lgkk.multiscale.rounds Non-negative integer number of compiled
#'   landmark-geodesic KK refinement rounds applied inside the multiscale solver
#'   after each eligible MISF level completes its standard GRIP rounds.
#'   This legacy shared budget is used as a fallback when any of the more
#'   specific per-stage budgets below are left \code{NULL}.
#' @param lgkk.rounds.coarse Optional non-negative integer number of compiled
#'   LGKK rounds applied on coarse MISF levels with \code{misf_level > 1}.
#'   When \code{NULL}, this falls back to \code{lgkk.multiscale.rounds}.
#' @param lgkk.rounds.pre.final Optional non-negative integer number of compiled
#'   LGKK rounds applied on the last coarse level just before the full graph is
#'   opened (\code{misf_level == 1}). When \code{NULL}, this falls back to
#'   \code{lgkk.multiscale.rounds}.
#' @param lgkk.rounds.final Optional non-negative integer number of compiled
#'   LGKK rounds applied after the full graph level completes its standard GRIP
#'   rounds (\code{misf_level == 0}). When \code{NULL}, this falls back to
#'   \code{lgkk.multiscale.rounds}.
#' @param lgkk.local.nbrs Number of nearest graph-metric neighbors retained per
#'   vertex in the LGKK sparse local set when either LGKK stage is enabled.
#' @param lgkk.landmark.count Number of farthest-point landmarks retained per
#'   vertex in the LGKK sparse long-range set when either LGKK stage is
#'   enabled.
#' @param lgkk.multiscale.scope Scope for the compiled multiscale LGKK stage.
#'   \code{"all"} applies it after every eligible MISF level, including the
#'   final full-graph level. \code{"coarse"} applies it only on coarse levels.
#' @param lgkk.active.limit Positive integer upper bound on the active-set size
#'   for compiled multiscale LGKK cache construction. Levels larger than this
#'   skip the multiscale LGKK stage.
#' @param metric.neighbor.cap Weighted-metric search limit used only when
#'   \code{metric = "edge_length"}. \code{NULL} performs the exact weighted
#'   neighborhood search and stops once the required neighbors and anchors are
#'   filled. A positive integer enables an approximate search by limiting the
#'   number of settled vertices per search. It is an error to supply this
#'   argument with \code{metric = "hop"}.
#' @param length.normalization Global normalization applied only when
#'   \code{metric = "edge_length"}: \code{"median"} (default) divides every
#'   edge length by their median, \code{"mean"} divides by their mean, and
#'   \code{"none"} preserves the supplied numerical scale. It is an error to
#'   supply this argument with \code{metric = "hop"}.
#' @param tinit.factor Initial temperature factor.
#' @param seed Optional RNG seed for reproducibility. If NULL, uses current time.
#' @param disconnected How to handle disconnected graphs:
#'   \code{"components"} (default) lays out each connected component separately
#'   and packs them into one coordinate matrix; \code{"error"} stops with an
#'   error.
#' @details
#' \strong{Edge-length semantics}
#'
#' The arguments \code{edge.weights} and \code{weight.list} are historically
#' named but represent positive edge \emph{lengths} or traversal costs, not
#' connection strengths, capacities, or similarities. A larger value requests
#' a longer geometric edge. If the available values are strengths for which a
#' larger value means a closer connection, convert them to positive lengths
#' before calling \code{grip()}, for example with a scientifically appropriate
#' reciprocal or other monotone decreasing transformation.
#'
#' With \code{metric = "hop"}, the standard GRIP hierarchy, insertion anchors,
#' and retained local neighborhoods use combinatorial shortest-path distance:
#' every traversed edge contributes one hop. When edge lengths are supplied,
#' they are nevertheless used as the desired lengths of adjacent vertex pairs
#' in the attractive force calculation. Thus this mode is useful when topology
#' should determine the multiscale organization but adjacent edges should have
#' unequal target lengths. No global length normalization is performed in this
#' mode. If no lengths are supplied, every edge has desired length one. Optional
#' LGKK stages use the supplied edge lengths for their geodesic distances even
#' though the standard GRIP stages remain hop based.
#'
#' With \code{metric = "edge_length"}, positive edge lengths are required. The
#' same lengths determine desired adjacent-edge lengths and the shortest-path
#' metric used for the MISF hierarchy, insertion anchors, retained graph
#' neighborhoods, and LGKK stages. Dijkstra-style weighted shortest paths are
#' used instead of hop-count breadth-first searches. By default, all lengths
#' are divided by their median before layout; this preserves relative geometry
#' while placing the numerical scale near the solver's unit scale. Use
#' \code{length.normalization = "mean"} for mean scaling or \code{"none"} when
#' the absolute supplied scale is intentional. Multiplying all input lengths
#' by the same positive constant therefore leaves the default normalized solve
#' unchanged.
#'
#' For \code{edges} input, provide one value per row through
#' \code{edge.weights}. For \code{adj.list} input, provide a parallel
#' \code{weight.list}: \code{weight.list[[i]][j]} is the length of the edge from
#' vertex \code{i} to \code{adj.list[[i]][j]}. For an undirected graph, the
#' adjacency and length entries should be symmetric.
#' @references
#' Gajer, P. and Kobourov, S.G. (2002). GRIP: Graph dRawing with Intelligent
#' Placement. \emph{Journal of Graph Algorithms and Applications}, 6(3),
#' 203--224. doi:10.7155/jgaa.00052.
#'
#' Gajer, P., Goodrich, M.T. and Kobourov, S.G. (2004). A multi-dimensional
#' approach to force-directed layouts of large graphs.
#' \emph{Computational Geometry}, 29(1), 3--18.
#' doi:10.1016/j.comgeo.2004.03.014.
#' @return A numeric matrix with `n` rows and `dim` columns.
#' @examples
#' edges <- edges.mesh(4, 4)
#' coords <- grip(edges, n = max(edges), dim = 2,
#'                       coarse.repulsion.factor = 0.2,
#'                       coarse.repulsion.sample = 8,
#'                       coarse.repulsion.exact.below = 32,
#'                       seed = 1)
#' round(coords, 3)
#'
#' # Use edge lengths throughout the multiscale graph metric.
#' path <- cbind(1:5, 2:6)
#' lengths <- c(1, 1, 2, 1, 1)
#' weighted.coords <- grip(
#'   path, n = 6, edge.weights = lengths,
#'   metric = "edge_length", dim = 2,
#'   rounds = 4, final.rounds = 4, num.init = 3, seed = 1
#' )
#' @export
#' @md
#' @section Workflow guides:
#' Start with \code{vignette("function-guide", package = "grip")}
#' to choose a layout, diagnostic, or reference comparison. List installed
#' guides with \code{vignette(package = "grip")}.
grip <- function(edges = NULL,
                        n = NULL,
                        adj.list = NULL,
                        weight.list = NULL,
                        edge.weights = NULL,
                        dim = 3,
                        placement = c("barycenter", "circle"),
                        preset = NULL,
                        rounds = 160,
                        final.rounds = 384,
                        num.init = 24,
                        num.nbrs = 20,
                        r = 0.03,
                        s = 7.5,
                        repulsion.factor = 2.5,
                        coarse.repulsion.factor = 1.5,
                        coarse.repulsion.sample = 16,
                        coarse.repulsion.exact.below = 64,
                        final.anchor.factor = 0,
                        final.move.scale.after.first = 1,
                        final.mode = c("fr", "kk_repulse"),
                        insertion.anchor.count = 3,
                        insertion.anchor.scope = c("any_higher", "prev_misf"),
                        insertion.anchor.strategy = c("first", "distance_band", "balanced_band", "spread_prev"),
                        level0.insertion.mode = c("inherit", "barycenter", "least_squares"),
                        level0.anchor.count = insertion.anchor.count,
                        level0.local.kk.steps = 3,
                        lgkk.polish.rounds = 0L,
                        lgkk.multiscale.rounds = 0L,
                        lgkk.rounds.coarse = NULL,
                        lgkk.rounds.pre.final = NULL,
                        lgkk.rounds.final = NULL,
                        lgkk.local.nbrs = 20L,
                        lgkk.landmark.count = 8L,
                        lgkk.multiscale.scope = c("all", "coarse"),
                        lgkk.active.limit = 4096L,
                        tinit.factor = 6,
                        seed = 6,
                        disconnected = c("components", "error"),
                        metric = c("hop", "edge_length"),
                        metric.neighbor.cap = NULL,
                        length.normalization = c("median", "mean", "none")) {
  grip.validate.graph.arguments(edges, n, adj.list, weight.list, edge.weights)
  metric_neighbor_cap_missing <- missing(metric.neighbor.cap)
  length_normalization_missing <- missing(length.normalization)
  metric <- match.arg(metric)

  call <- match.call(expand.dots = FALSE)
  call$metric <- NULL

  if (identical(metric, "hop")) {
    if (!metric_neighbor_cap_missing) {
      stop("metric_neighbor_cap is only available when metric = 'edge_length'", call. = FALSE)
    }
    if (!length_normalization_missing) {
      stop("length_normalization is only available when metric = 'edge_length'", call. = FALSE)
    }
    call$metric_neighbor_cap <- NULL
    call$length_normalization <- NULL
    return(grip.forward.call(globalrep.grip, call, env = parent.frame()))
  }

  grip.forward.call(globalrep.weighted.grip, call, env = parent.frame())
}

#' Compute the legacy GRIP layout
#'
#' This function preserves the original local-force wrapper and historical
#' default values that were previously exposed as \code{\link{grip}()}.
#' Use it for backwards-compatible comparisons or when you explicitly want the
#' pre-global-repulsion behavior.
#'
#' @param edges Two-column integer matrix of edges (1-based vertex ids).
#'   Supply either \code{edges}/\code{edge.weights} or
#'   \code{adj.list}/\code{weight.list}, not both. Fractional, nonfinite,
#'   and out-of-range vertex ids are rejected before integer conversion.
#' @param n Number of vertices, a finite positive integer.
#' @param adj.list Adjacency list (1-based) for undirected graphs.
#' @param weight.list Optional parallel list of edge weights (edge lengths).
#'   If NULL, all edges are treated as weight 1. All weights must be finite
#'   and strictly positive.
#' @param edge.weights Optional vector of edge weights for \code{edges}. All
#'   weights must be finite and strictly positive.
#' @param dim Layout dimension (2 or 3). Default is 3.
#' @param placement Initial placement strategy. "circle" is only used for 2D.
#' @param preset Optional tuning preset. \code{NULL} uses the historical
#'   defaults. \code{"carpet"} applies a preset tuned for
#'   Sierpinski-carpet-like graphs and validated on carpet levels 3 and 4.
#'   \code{"mesh"} applies a preset tuned for rectangular lattice graphs and
#'   validated on 8x8 and 12x12 mesh layouts. \code{"torus"} applies a preset
#'   tuned for 3D torus layouts and validated on torus sizes from 8x8 through
#'   20x20. \code{"tree"} applies a preset tuned for symmetric force-directed
#'   layouts of tree-like graphs and validated on binary trees of depths 5 and
#'   6. Presets only fill in tuning arguments that you did not supply
#'   explicitly.
#' @param rounds Initial rounds for refinement.
#' @param final.rounds Final rounds for refinement.
#' @param num.init Number of initial vertices in the coarsest level.
#' @param num.nbrs Maximum number of graph-distance neighbors retained for local
#'   refinement at each filtration level.
#' @param r Main local temperature adaptation rate in \code{[0, 1]}.
#' @param s Non-negative boost factor applied when successive displacements have
#'   a consistent direction.
#' @param repulsion.factor Non-negative multiplier applied to GRIP's
#'   finest-level repulsive force scale. \code{1} keeps the historical
#'   repulsion strength; \code{0} disables that repulsive term.
#' @param tinit.factor Initial temperature factor.
#' @param seed Optional RNG seed for reproducibility. If NULL, uses current time.
#' @param disconnected How to handle disconnected graphs:
#'   \code{"components"} (default) lays out each connected component separately
#'   and packs them into one coordinate matrix; \code{"error"} stops with an error.
#' @references
#' Gajer, P. and Kobourov, S.G. (2002). GRIP: Graph dRawing with Intelligent
#' Placement. \emph{Journal of Graph Algorithms and Applications}, 6(3),
#' 203--224. doi:10.7155/jgaa.00052.
#'
#' Gajer, P., Goodrich, M.T. and Kobourov, S.G. (2004). A multi-dimensional
#' approach to force-directed layouts of large graphs.
#' \emph{Computational Geometry}, 29(1), 3--18.
#' doi:10.1016/j.comgeo.2004.03.014.
#' @return A numeric matrix with `n` rows and `dim` columns.
#' @examples
#' edges <- cbind(1:5, 2:6)
#' coords <- legacy.grip(edges, n = 6, dim = 2,
#'                              placement = "barycenter",
#'                              rounds = 5, final.rounds = 5,
#'                              num.init = 3, num.nbrs = 4,
#'                              seed = 1)
#' round(coords, 3)
#' @export
#' @md
legacy.grip <- function(edges = NULL,
                               n = NULL,
                               adj.list = NULL,
                               weight.list = NULL,
                               edge.weights = NULL,
                               dim = 3,
                               placement = c("barycenter", "circle"),
                               preset = NULL,
                               rounds = 20,
                               final.rounds = 25,
                               num.init = 36,
                               num.nbrs = 10,
                               r = 0.15,
                               s = 3.0,
                               repulsion.factor = 1.0,
                               tinit.factor = 6,
                               seed = 6,
                               disconnected = c("components", "error")) {
  grip.validate.graph.arguments(edges, n, adj.list, weight.list, edge.weights)
  placement_missing <- missing(placement)
  rounds_missing <- missing(rounds)
  final_rounds_missing <- missing(final.rounds)
  num_init_missing <- missing(num.init)
  num_nbrs_missing <- missing(num.nbrs)
  r_missing <- missing(r)
  s_missing <- missing(s)
  repulsion_factor_missing <- missing(repulsion.factor)

  preset <- grip.normalize.preset(preset, fn = "legacy.grip")
  resolved <- grip.resolve.preset(
    preset = preset,
    dim = dim,
    placement = placement,
    placement.missing = placement_missing,
    rounds = rounds,
    rounds.missing = rounds_missing,
    final.rounds = final.rounds,
    final.rounds.missing = final_rounds_missing,
    num.init = num.init,
    num.init.missing = num_init_missing,
    num.nbrs = num.nbrs,
    num.nbrs.missing = num_nbrs_missing,
    r = r,
    r.missing = r_missing,
    s = s,
    s.missing = s_missing,
    repulsion.factor = repulsion.factor,
    repulsion.factor.missing = repulsion_factor_missing
  )
  placement <- resolved$placement
  rounds <- resolved$rounds
  final.rounds <- resolved$final_rounds
  num.init <- resolved$num_init
  num.nbrs <- resolved$num_nbrs
  r <- resolved$r
  s <- resolved$s
  repulsion.factor <- resolved$repulsion_factor
  placement <- match.arg(placement)
  disconnected <- match.arg(disconnected)

  validated <- grip.validate.layout.inputs(
    edges = edges,
    n = n,
    adj.list = adj.list,
    weight.list = weight.list,
    edge.weights = edge.weights,
    dim = dim,
    placement = placement,
    seed = seed
  )
  adj.list <- validated$adj_list
  weight.list <- validated$weight_list
  n <- validated$n
  dim <- validated$dim
  seed <- validated$seed
  tuning <- grip.validate.tuning.inputs(
    num.nbrs = num.nbrs,
    r = r,
    s = s,
    repulsion.factor = repulsion.factor
  )
  num.nbrs <- tuning$num_nbrs
  r <- tuning$r
  s <- tuning$s
  repulsion.factor <- tuning$repulsion_factor

  layout.adj <- function(adj.list, weight.list, n) {
    grip_layout_adj_cpp(adj_list = adj.list,
                        weight_list = weight.list,
                        n = n,
                        dim = dim,
                        placement = placement,
                        rounds = as.integer(rounds),
                        final_rounds = as.integer(final.rounds),
                        num_init = as.integer(num.init),
                        num_nbrs = num.nbrs,
                        r = r,
                        s = s,
                        repulsion_factor = repulsion.factor,
                        tinit_factor = as.integer(tinit.factor),
                        seed = seed)
  }

  comp <- grip.connected.components(adj.list = adj.list, n = n)
  n.comp <- length(unique(comp))

  if (n.comp == 1L) {
    return(layout.adj(adj.list = adj.list, weight.list = weight.list, n = n))
  }

  if (identical(disconnected, "error")) {
    stop(sprintf(
      "Input graph has %d connected components; the GRIP layout core assumes connected graphs. Use disconnected = 'components' to lay out each component safely.",
      n.comp
    ))
  }

  warning(
    sprintf(
      "Input graph has %d connected components; laying out components separately to avoid disconnected-graph instability.",
      n.comp
    ),
    call. = FALSE
  )

  comp.ids <- sort(unique(comp))
  layouts <- vector("list", length(comp.ids))
  for (k in seq_along(comp.ids)) {
    rows <- which(comp == comp.ids[[k]])
    sub <- grip.induce.subgraph(
      adj.list = adj.list,
      weight.list = weight.list,
      vertices = rows,
      n = n
    )
    layouts[[k]] <- layout.adj(
      adj.list = sub$adj_list,
      weight.list = sub$weight_list,
      n = length(rows)
    )
  }
  grip.pack.component.layouts(layouts = layouts, comp = comp, n = n, dim = dim)
}

#' Compute a GRIP layout trace
#'
#' This traces the primary \code{\link{grip}()} engine, including the
#' coarse-level global-repulsion term used by the quality-first default layout.
#' For backwards-compatible traces of the historical local-force wrapper, use
#' \code{\link{trace.legacy.grip}()}.
#'
#' @inheritParams grip
#' @param trace Snapshot granularity. \code{"round"} records the coarsest
#'   initialization, each level start, every \code{trace.every} completed rounds,
#'   and the final layout. \code{"level"} records the coarsest initialization,
#'   every \code{trace.every}th level start, and the final layout.
#' @param trace.every Positive integer thinning factor for recorded rounds or
#'   levels. Initial and final snapshots are always included.
#' @param diagnostics Optional per-frame diagnostic mode. \code{"none"} skips
#'   extra scoring, \code{"light"} appends lightweight shape diagnostics, and
#'   \code{"full"} also computes sampled stress on each traced frame.
#' @param target.coords Optional numeric target coordinate matrix used to append
#'   per-frame Procrustes RMSE diagnostics. It must have `n` rows and `dim`
#'   columns.
#' @param diagnostic.sample.size.nonedge Positive integer sample size used for
#'   per-frame non-edge separation diagnostics when \code{diagnostics != "none"}.
#' @param diagnostic.sample.size.stress Positive integer sample size used for
#'   per-frame sampled stress when \code{diagnostics = "full"}.
#' @param diagnostic.nonedge.seed RNG seed base used for per-frame non-edge
#'   separation diagnostics.
#' @param diagnostic.stress.seed RNG seed base used for per-frame sampled stress
#'   diagnostics.
#' @return A list with \code{final}, \code{frames}, \code{meta}, \code{trace},
#'   \code{trace.every}, canonical \code{stage_trace} and \code{stage_data},
#'   and optionally \code{diagnostics} and \code{lgkk.polish}. \code{final} is
#'   the final coordinate matrix. \code{frames} is a list of coordinate
#'   matrices with \code{NA} rows for vertices that have not yet been
#'   introduced by GRIP. \code{meta} is a data frame describing each frame with
#'   columns \code{frame}, \code{phase}, \code{level_index}, \code{misf_level},
#'   \code{round_in_level}, and \code{active_vertices}. \code{stage_trace} and
#'   \code{stage_data} lift the raw frames onto the shared MISF stage schema
#'   used by the geodesic multiscale methods, with canonical states such as
#'   \code{seed}, \code{initial_placement}, \code{top_level},
#'   \code{insertion}, \code{refinement}, and \code{final_polish}. When
#'   diagnostics are requested, \code{diagnostics} is a data frame parallel to
#'   \code{meta} that appends per-frame quality metrics such as
#'   \code{edge.length.cv}, \code{sampled.nonedge.sep.ratio}, and optional
#'   \code{procrustes.rmse}.
#' @examples
#' edges <- cbind(1:5, 2:6)
#' tr <- trace.grip(edges, n = 6, dim = 2,
#'                         placement = "barycenter",
#'                         rounds = 3, final.rounds = 2,
#'                         num.init = 3, num.nbrs = 4,
#'                         trace = "level",
#'                         trace.every = 1,
#'                         diagnostics = "light",
#'                         seed = 1)
#' tr$diagnostics
#' @noRd
grip.trace.hop <- function(edges = NULL,
                              n = NULL,
                              adj.list = NULL,
                              weight.list = NULL,
                              edge.weights = NULL,
                              dim = 3,
                              placement = c("barycenter", "circle"),
                              preset = NULL,
                              rounds = 160,
                              final.rounds = 384,
                              num.init = 24,
                              num.nbrs = 20,
                              r = 0.03,
                              s = 7.5,
                              repulsion.factor = 2.5,
                              coarse.repulsion.factor = 1.5,
                              coarse.repulsion.sample = 16,
                              coarse.repulsion.exact.below = 64,
                              final.anchor.factor = 0,
                              final.move.scale.after.first = 1,
                              final.mode = c("fr", "kk_repulse"),
                              insertion.anchor.count = 3,
                              insertion.anchor.scope = c("any_higher", "prev_misf"),
                              insertion.anchor.strategy = c("first", "distance_band", "balanced_band", "spread_prev"),
                              level0.insertion.mode = c("inherit", "barycenter", "least_squares"),
                              level0.anchor.count = insertion.anchor.count,
                              level0.local.kk.steps = 3,
                              lgkk.polish.rounds = 0L,
                              lgkk.multiscale.rounds = 0L,
                              lgkk.rounds.coarse = NULL,
                              lgkk.rounds.pre.final = NULL,
                              lgkk.rounds.final = NULL,
                              lgkk.local.nbrs = 20L,
                              lgkk.landmark.count = 8L,
                              lgkk.multiscale.scope = c("all", "coarse"),
                              lgkk.active.limit = 4096L,
                              tinit.factor = 6,
                              seed = 6,
                              trace = c("round", "level"),
                              trace.every = 1,
                              diagnostics = c("none", "light", "full"),
                              target.coords = NULL,
                              diagnostic.sample.size.nonedge = 1000L,
                              diagnostic.sample.size.stress = 500L,
                              diagnostic.nonedge.seed = 1L,
                              diagnostic.stress.seed = 1L) {
  grip.validate.graph.arguments(edges, n, adj.list, weight.list, edge.weights)
  placement_missing <- missing(placement)
  rounds_missing <- missing(rounds)
  final_rounds_missing <- missing(final.rounds)
  num_init_missing <- missing(num.init)
  num_nbrs_missing <- missing(num.nbrs)
  r_missing <- missing(r)
  s_missing <- missing(s)
  repulsion_factor_missing <- missing(repulsion.factor)

  preset <- grip.normalize.preset(preset, fn = "trace.grip")
  resolved <- grip.resolve.preset(
    preset = preset,
    dim = dim,
    placement = placement,
    placement.missing = placement_missing,
    rounds = rounds,
    rounds.missing = rounds_missing,
    final.rounds = final.rounds,
    final.rounds.missing = final_rounds_missing,
    num.init = num.init,
    num.init.missing = num_init_missing,
    num.nbrs = num.nbrs,
    num.nbrs.missing = num_nbrs_missing,
    r = r,
    r.missing = r_missing,
    s = s,
    s.missing = s_missing,
    repulsion.factor = repulsion.factor,
    repulsion.factor.missing = repulsion_factor_missing
  )
  placement <- resolved$placement
  rounds <- resolved$rounds
  final.rounds <- resolved$final_rounds
  num.init <- resolved$num_init
  num.nbrs <- resolved$num_nbrs
  r <- resolved$r
  s <- resolved$s
  repulsion.factor <- resolved$repulsion_factor
  placement <- match.arg(placement)
  final.mode <- match.arg(final.mode)
  insertion.anchor.scope <- match.arg(insertion.anchor.scope)
  insertion.anchor.strategy <- match.arg(insertion.anchor.strategy)
  level0.insertion.mode <- match.arg(level0.insertion.mode)
  lgkk.multiscale.scope <- match.arg(lgkk.multiscale.scope)
  trace <- match.arg(trace)
  diagnostics <- match.arg(diagnostics)

  if (!is.numeric(trace.every) || length(trace.every) != 1L || !is.finite(trace.every)) {
    stop("trace.every must be a single finite numeric value")
  }
  trace.every <- as.integer(trace.every)
  if (is.na(trace.every) || trace.every <= 0L) {
    stop("trace.every must be a positive integer")
  }

  validated <- grip.validate.layout.inputs(
    edges = edges,
    n = n,
    adj.list = adj.list,
    weight.list = weight.list,
    edge.weights = edge.weights,
    dim = dim,
    placement = placement,
    seed = seed
  )
  adj.list <- validated$adj_list
  weight.list <- validated$weight_list
  n <- validated$n
  dim <- validated$dim
  seed <- validated$seed
  trace.edges <- grip.edges.from.adj.list(adj.list)
  if (is.null(preset) && final_rounds_missing) {
    final.rounds <- grip.globalrep.default.final.rounds(n)
  }
  tuning <- grip.validate.globalrep.tuning.inputs(
    num.nbrs = num.nbrs,
    r = r,
    s = s,
    repulsion.factor = repulsion.factor,
    coarse.repulsion.factor = coarse.repulsion.factor,
    coarse.repulsion.sample = coarse.repulsion.sample,
    coarse.repulsion.exact.below = coarse.repulsion.exact.below,
    final.anchor.factor = final.anchor.factor,
    final.move.scale.after.first = final.move.scale.after.first,
    insertion.anchor.count = insertion.anchor.count,
    insertion.anchor.scope = insertion.anchor.scope,
    insertion.anchor.strategy = insertion.anchor.strategy,
    level0.insertion.mode = level0.insertion.mode,
    level0.anchor.count = level0.anchor.count,
    level0.local.kk.steps = level0.local.kk.steps
  )
  num.nbrs <- tuning$num_nbrs
  r <- tuning$r
  s <- tuning$s
  repulsion.factor <- tuning$repulsion_factor
  coarse.repulsion.factor <- tuning$coarse_repulsion_factor
  coarse.repulsion.sample <- tuning$coarse_repulsion_sample
  coarse.repulsion.exact.below <- tuning$coarse_repulsion_exact_below
  final.anchor.factor <- tuning$final_anchor_factor
  final.move.scale.after.first <- tuning$final_move_scale_after_first
  insertion.anchor.count <- tuning$insertion_anchor_count
  insertion.anchor.scope <- tuning$insertion_anchor_scope
  insertion.anchor.strategy <- tuning$insertion_anchor_strategy
  level0.insertion.mode <- tuning$level0_insertion_mode
  level0.anchor.count <- tuning$level0_anchor_count
  level0.local.kk.steps <- tuning$level0_local_kk_steps
  lgkk <- grip.validate.lgkk.polish.inputs(
    lgkk.polish.rounds = lgkk.polish.rounds,
    lgkk.multiscale.rounds = lgkk.multiscale.rounds,
    lgkk.rounds.coarse = lgkk.rounds.coarse,
    lgkk.rounds.pre.final = lgkk.rounds.pre.final,
    lgkk.rounds.final = lgkk.rounds.final,
    lgkk.local.nbrs = lgkk.local.nbrs,
    lgkk.landmark.count = lgkk.landmark.count,
    lgkk.multiscale.scope = lgkk.multiscale.scope,
    lgkk.active.limit = lgkk.active.limit
  )
  lgkk.polish.rounds <- lgkk$lgkk_polish_rounds
  lgkk.multiscale.rounds <- lgkk$lgkk_multiscale_rounds
  lgkk.rounds.coarse <- lgkk$lgkk_rounds_coarse
  lgkk.rounds.pre.final <- lgkk$lgkk_rounds_pre_final
  lgkk.rounds.final <- lgkk$lgkk_rounds_final
  lgkk.local.nbrs <- lgkk$lgkk_local_nbrs
  lgkk.landmark.count <- lgkk$lgkk_landmark_count
  lgkk.multiscale.scope <- lgkk$lgkk_multiscale_scope
  lgkk.active.limit <- lgkk$lgkk_active_limit

  comp <- grip.connected.components(adj.list = adj.list, n = n)
  n.comp <- length(unique(comp))
  if (n.comp != 1L) {
    stop(sprintf(
      "trace.grip() currently supports only connected graphs; input graph has %d connected components.",
      n.comp
    ))
  }

  out <- grip_layout_globalrep_trace_adj_cpp(
    adj_list = adj.list,
    weight_list = weight.list,
    n = n,
    dim = dim,
    placement = placement,
    rounds = as.integer(rounds),
    final_rounds = as.integer(final.rounds),
    num_init = as.integer(num.init),
    num_nbrs = num.nbrs,
    r = r,
    s = s,
    repulsion_factor = repulsion.factor,
    coarse_repulsion_factor = coarse.repulsion.factor,
    coarse_repulsion_sample = coarse.repulsion.sample,
    coarse_repulsion_exact_below = coarse.repulsion.exact.below,
    final_anchor_factor = final.anchor.factor,
    final_move_scale_after_first = final.move.scale.after.first,
    insertion_anchor_count = insertion.anchor.count,
    insertion_anchor_scope = insertion.anchor.scope,
    insertion_anchor_strategy = insertion.anchor.strategy,
    level0_insertion_mode = level0.insertion.mode,
    level0_anchor_count = level0.anchor.count,
    level0_local_kk_steps = level0.local.kk.steps,
    lgkk_multiscale_rounds = lgkk.multiscale.rounds,
    lgkk_rounds_coarse = lgkk.rounds.coarse,
    lgkk_rounds_pre_final = lgkk.rounds.pre.final,
    lgkk_rounds_final = lgkk.rounds.final,
    lgkk_local_nbrs = lgkk.local.nbrs,
    lgkk_landmark_count = lgkk.landmark.count,
    lgkk_multiscale_scope = lgkk.multiscale.scope,
    lgkk_active_limit = lgkk.active.limit,
    final_mode = final.mode,
    tinit_factor = as.integer(tinit.factor),
    seed = seed,
    trace = trace,
    trace_every = trace.every
  )
  if (lgkk.polish.rounds > 0L) {
    polished <- grip.apply.lgkk.polish(
      coords = out$final,
      adj.list = adj.list,
      weight.list = weight.list,
      rounds = lgkk.polish.rounds,
      lgkk.local.nbrs = lgkk.local.nbrs,
      lgkk.landmark.count = lgkk.landmark.count,
      return.trace = TRUE
    )
    if (length(polished$frames) > 1L) {
      add.frames <- polished$frames[-1L]
      add.meta <- data.frame(
        frame = seq.int(nrow(out$meta) + 1L, nrow(out$meta) + length(add.frames)),
        phase = rep("lgkk", length(add.frames)),
        level_index = rep(utils::tail(out$meta$level_index, 1L), length(add.frames)),
        misf_level = rep(utils::tail(out$meta$misf_level, 1L), length(add.frames)),
        round_in_level = seq_len(length(add.frames)),
        active_vertices = rep(n, length(add.frames)),
        stringsAsFactors = FALSE
      )
      out$frames <- c(out$frames, add.frames)
      out$meta <- rbind(out$meta, add.meta)
    } else {
      out$frames[[length(out$frames)]] <- polished$coords
    }
    out$final <- polished$coords
    out$lgkk.polish <- polished$trace
  } else {
    out$lgkk.polish <- data.frame()
  }
  out$trace <- trace
  out$trace.every <- trace.every
  out$diagnostics <- grip.trace.compute.diagnostics(
    frames = out$frames,
    meta = out$meta,
    adj.list = adj.list,
    weight.list = weight.list,
    diagnostics = diagnostics,
    target.coords = target.coords,
    sample.size.nonedge = diagnostic.sample.size.nonedge,
    sample.size.stress = diagnostic.sample.size.stress,
    nonedge.seed = diagnostic.nonedge.seed,
    stress.seed = diagnostic.stress.seed
  )
  stage.bundle <- grip.layout.trace.as.stage.bundle(
    trace = out,
    edges = trace.edges
  )
  out$stage_trace <- stage.bundle$stage_trace
  out$stage_data <- stage.bundle$stage_data
  class(out) <- c("grip_layout_trace", class(out))
  out
}

#' Trace a GRIP layout
#'
#' This is the tracing counterpart to \code{\link{grip}()}. The
#' \code{metric} argument selects the hop-based or edge-length-based compiled
#' engine while the remaining trace and diagnostic options keep the same
#' meaning in both cases. See \code{\link{grip}()} for detailed edge-length
#' semantics.
#'
#' @inheritParams grip
#' @param metric Graph metric traced by the multiscale engine:
#'   \code{"hop"} (default) or \code{"edge_length"}. See
#'   \code{\link{grip}()} for the exact role and normalization of supplied edge
#'   lengths in each mode.
#' @param trace Snapshot granularity. \code{"round"} records the coarsest
#'   initialization, each level start, every \code{trace.every} completed
#'   rounds, and the final layout. \code{"level"} records the coarsest
#'   initialization, every \code{trace.every}th level start, and the final
#'   layout.
#' @param trace.every Positive integer thinning factor for recorded rounds or
#'   levels. Initial and final snapshots are always included.
#' @param diagnostics Optional per-frame diagnostic mode. \code{"none"} skips
#'   extra scoring, \code{"light"} appends lightweight shape diagnostics, and
#'   \code{"full"} also computes sampled stress on each traced frame.
#' @param target.coords Optional numeric target coordinate matrix used to append
#'   per-frame Procrustes RMSE diagnostics. It must have \code{n} rows and
#'   \code{dim} columns.
#' @param diagnostic.sample.size.nonedge Positive integer sample size used for
#'   per-frame non-edge separation diagnostics when
#'   \code{diagnostics != "none"}.
#' @param diagnostic.sample.size.stress Positive integer sample size used for
#'   per-frame sampled stress when \code{diagnostics = "full"}.
#' @param diagnostic.nonedge.seed RNG seed base used for per-frame non-edge
#'   separation diagnostics.
#' @param diagnostic.stress.seed RNG seed base used for per-frame sampled stress
#'   diagnostics.
#' @return A list containing the final layout, recorded coordinate frames,
#'   frame metadata, trace settings, canonical stage data, and any requested
#'   diagnostics.
#' @examples
#' edges <- cbind(1:5, 2:6)
#' tr <- trace.grip(
#'   edges, n = 6, metric = "hop", dim = 2,
#'   rounds = 3, final.rounds = 2, num.init = 3,
#'   trace = "level", diagnostics = "light", seed = 1
#' )
#' tr$meta
#' @export
#' @section Workflow guides:
#' Start with \code{vignette("function-guide", package = "grip")}
#' to choose a layout, diagnostic, or reference comparison. List installed
#' guides with \code{vignette(package = "grip")}.
trace.grip <- function(edges = NULL,
                       n = NULL,
                       adj.list = NULL,
                       weight.list = NULL,
                       edge.weights = NULL,
                       dim = 3,
                       placement = c("barycenter", "circle"),
                       preset = NULL,
                       rounds = 160,
                       final.rounds = 384,
                       num.init = 24,
                       num.nbrs = 20,
                       r = 0.03,
                       s = 7.5,
                       repulsion.factor = 2.5,
                       coarse.repulsion.factor = 1.5,
                       coarse.repulsion.sample = 16,
                       coarse.repulsion.exact.below = 64,
                       final.anchor.factor = 0,
                       final.move.scale.after.first = 1,
                       final.mode = c("fr", "kk_repulse"),
                       insertion.anchor.count = 3,
                       insertion.anchor.scope = c("any_higher", "prev_misf"),
                       insertion.anchor.strategy = c("first", "distance_band", "balanced_band", "spread_prev"),
                       level0.insertion.mode = c("inherit", "barycenter", "least_squares"),
                       level0.anchor.count = insertion.anchor.count,
                       level0.local.kk.steps = 3,
                       lgkk.polish.rounds = 0L,
                       lgkk.multiscale.rounds = 0L,
                       lgkk.rounds.coarse = NULL,
                       lgkk.rounds.pre.final = NULL,
                       lgkk.rounds.final = NULL,
                       lgkk.local.nbrs = 20L,
                       lgkk.landmark.count = 8L,
                       lgkk.multiscale.scope = c("all", "coarse"),
                       lgkk.active.limit = 4096L,
                       tinit.factor = 6,
                       seed = 6,
                       trace = c("round", "level"),
                       trace.every = 1,
                       diagnostics = c("none", "light", "full"),
                       target.coords = NULL,
                       diagnostic.sample.size.nonedge = 1000L,
                       diagnostic.sample.size.stress = 500L,
                       diagnostic.nonedge.seed = 1L,
                       diagnostic.stress.seed = 1L,
                       metric = c("hop", "edge_length"),
                       metric.neighbor.cap = NULL,
                       length.normalization = c("median", "mean", "none")) {
  grip.validate.graph.arguments(edges, n, adj.list, weight.list, edge.weights)
  metric_neighbor_cap_missing <- missing(metric.neighbor.cap)
  length_normalization_missing <- missing(length.normalization)
  metric <- match.arg(metric)

  call <- match.call(expand.dots = FALSE)
  call$metric <- NULL

  if (identical(metric, "hop")) {
    if (!metric_neighbor_cap_missing) {
      stop("metric_neighbor_cap is only available when metric = 'edge_length'", call. = FALSE)
    }
    if (!length_normalization_missing) {
      stop("length_normalization is only available when metric = 'edge_length'", call. = FALSE)
    }
    call$metric_neighbor_cap <- NULL
    call$length_normalization <- NULL
    return(grip.forward.call(grip.trace.hop, call, env = parent.frame()))
  }

  grip.forward.call(grip.trace.edge.length, call, env = parent.frame())
}

#' Compute a trace for the legacy GRIP layout
#'
#' @inheritParams legacy.grip
#' @param trace Snapshot granularity. \code{"round"} records the coarsest
#'   initialization, each level start, every \code{trace.every} completed rounds,
#'   and the final layout. \code{"level"} records the coarsest initialization,
#'   every \code{trace.every}th level start, and the final layout.
#' @param trace.every Positive integer thinning factor for recorded rounds or
#'   levels. Initial and final snapshots are always included.
#' @return A list with \code{final}, \code{frames}, \code{meta}, \code{trace},
#'   and \code{trace.every}. \code{final} is the final coordinate matrix.
#'   \code{frames} is a list of coordinate matrices with \code{NA} rows for
#'   vertices that have not yet been introduced by GRIP. \code{meta} is a data
#'   frame describing each frame with columns \code{frame}, \code{phase},
#'   \code{level_index}, \code{misf_level}, \code{round_in_level}, and
#'   \code{active_vertices}.
#' @examples
#' edges <- cbind(1:5, 2:6)
#' tr <- trace.legacy.grip(edges, n = 6, dim = 2,
#'                                placement = "barycenter",
#'                                rounds = 3, final.rounds = 2,
#'                                num.init = 3, num.nbrs = 4,
#'                                trace = "level",
#'                                trace.every = 1,
#'                                seed = 1)
#' tr$meta
#' @export
trace.legacy.grip <- function(edges = NULL,
                                     n = NULL,
                                     adj.list = NULL,
                                     weight.list = NULL,
                                     edge.weights = NULL,
                                     dim = 3,
                                     placement = c("barycenter", "circle"),
                                     preset = NULL,
                                     rounds = 20,
                                     final.rounds = 25,
                                     num.init = 36,
                                     num.nbrs = 10,
                                     r = 0.15,
                                     s = 3.0,
                                     repulsion.factor = 1.0,
                                     tinit.factor = 6,
                                     seed = 6,
                                     trace = c("round", "level"),
                                     trace.every = 1) {
  grip.validate.graph.arguments(edges, n, adj.list, weight.list, edge.weights)
  placement_missing <- missing(placement)
  rounds_missing <- missing(rounds)
  final_rounds_missing <- missing(final.rounds)
  num_init_missing <- missing(num.init)
  num_nbrs_missing <- missing(num.nbrs)
  r_missing <- missing(r)
  s_missing <- missing(s)
  repulsion_factor_missing <- missing(repulsion.factor)

  preset <- grip.normalize.preset(preset, fn = "trace.legacy.grip")
  resolved <- grip.resolve.preset(
    preset = preset,
    dim = dim,
    placement = placement,
    placement.missing = placement_missing,
    rounds = rounds,
    rounds.missing = rounds_missing,
    final.rounds = final.rounds,
    final.rounds.missing = final_rounds_missing,
    num.init = num.init,
    num.init.missing = num_init_missing,
    num.nbrs = num.nbrs,
    num.nbrs.missing = num_nbrs_missing,
    r = r,
    r.missing = r_missing,
    s = s,
    s.missing = s_missing,
    repulsion.factor = repulsion.factor,
    repulsion.factor.missing = repulsion_factor_missing
  )
  placement <- resolved$placement
  rounds <- resolved$rounds
  final.rounds <- resolved$final_rounds
  num.init <- resolved$num_init
  num.nbrs <- resolved$num_nbrs
  r <- resolved$r
  s <- resolved$s
  repulsion.factor <- resolved$repulsion_factor
  placement <- match.arg(placement)
  trace <- match.arg(trace)

  if (!is.numeric(trace.every) || length(trace.every) != 1L || !is.finite(trace.every)) {
    stop("trace.every must be a single finite numeric value")
  }
  trace.every <- as.integer(trace.every)
  if (is.na(trace.every) || trace.every <= 0L) {
    stop("trace.every must be a positive integer")
  }

  validated <- grip.validate.layout.inputs(
    edges = edges,
    n = n,
    adj.list = adj.list,
    weight.list = weight.list,
    edge.weights = edge.weights,
    dim = dim,
    placement = placement,
    seed = seed
  )
  adj.list <- validated$adj_list
  weight.list <- validated$weight_list
  n <- validated$n
  dim <- validated$dim
  seed <- validated$seed
  tuning <- grip.validate.tuning.inputs(
    num.nbrs = num.nbrs,
    r = r,
    s = s,
    repulsion.factor = repulsion.factor
  )
  num.nbrs <- tuning$num_nbrs
  r <- tuning$r
  s <- tuning$s
  repulsion.factor <- tuning$repulsion_factor

  comp <- grip.connected.components(adj.list = adj.list, n = n)
  n.comp <- length(unique(comp))
  if (n.comp != 1L) {
    stop(sprintf(
      "trace.legacy.grip() currently supports only connected graphs; input graph has %d connected components.",
      n.comp
    ))
  }

  out <- grip_layout_trace_adj_cpp(
    adj_list = adj.list,
    weight_list = weight.list,
    n = n,
    dim = dim,
    placement = placement,
    rounds = as.integer(rounds),
    final_rounds = as.integer(final.rounds),
    num_init = as.integer(num.init),
    num_nbrs = num.nbrs,
    r = r,
    s = s,
    repulsion_factor = repulsion.factor,
    tinit_factor = as.integer(tinit.factor),
    seed = seed,
    trace = trace,
    trace_every = trace.every
  )
  out$trace <- trace
  out$trace.every <- trace.every
  class(out) <- c("grip_layout_trace", class(out))
  out
}
