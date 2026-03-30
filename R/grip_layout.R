grip.build.adj.from.edges <- function(edges, n, edge_weights = NULL) {
  edges <- as.matrix(edges)
  if (!is.numeric(edges) || ncol(edges) != 2) {
    stop("edges must be a two-column integer matrix of 1-based vertex ids")
  }
  edges <- matrix(as.integer(edges), ncol = 2)
  if (any(!is.finite(edges))) {
    stop("edges must contain finite integer vertex ids")
  }
  if (any(edges <= 0L | edges > n)) {
    stop("edges must be 1-based and within [1, n]")
  }

  use.weights <- !is.null(edge_weights)
  if (use.weights) {
    if (length(edge_weights) != nrow(edges)) {
      stop("edge_weights length must match number of edges")
    }
    if (!is.numeric(edge_weights)) {
      stop("edge_weights must be a numeric vector")
    }
    bad <- which(!is.finite(edge_weights) | edge_weights <= 0)
    if (length(bad) > 0L) {
      i <- bad[[1L]]
      stop(sprintf(
        "edge_weights must contain finite values > 0; first invalid at edge_weights[%d] = %s",
        i,
        format(edge_weights[i], digits = 16)
      ))
    }
    edge_weights <- as.double(edge_weights)
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
      w <- edge_weights[[i]]
      out.w[[u]] <- c(out.w[[u]], w)
      out.w[[v]] <- c(out.w[[v]], w)
    }
  }

  list(adj_list = out.adj, weight_list = out.w)
}

grip.connected.components <- function(adj_list, n) {
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
      nb <- adj_list[[x]]
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

grip.induce.subgraph <- function(adj_list, weight_list, vertices, n) {
  vertices <- as.integer(vertices)
  in.comp <- rep(FALSE, n)
  in.comp[vertices] <- TRUE
  map <- integer(n)
  map[vertices] <- seq_along(vertices)
  sub.adj <- vector("list", length(vertices))
  use.weights <- !is.null(weight_list)
  sub.w <- if (use.weights) vector("list", length(vertices)) else NULL

  for (i in seq_along(vertices)) {
    v <- vertices[[i]]
    nb <- as.integer(adj_list[[v]])
    if (length(nb) == 0L) {
      sub.adj[[i]] <- integer(0)
      if (use.weights) sub.w[[i]] <- numeric(0)
      next
    }
    keep <- in.comp[nb]
    nb.keep <- nb[keep]
    sub.adj[[i]] <- as.integer(map[nb.keep])
    if (use.weights) {
      sub.w[[i]] <- as.double(weight_list[[v]][keep])
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

grip.normalize.preset <- function(preset, fn = "grip.layout") {
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

grip.globalrep.default.final_rounds <- function(n) {
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
    final_rounds = if (is.null(n)) 384L else grip.globalrep.default.final_rounds(n),
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

grip.forward_call <- function(target, call, env = parent.frame()) {
  args <- as.list(call)[-1L]
  do.call(target, args, envir = env)
}

grip.resolve.preset <- function(preset,
                                dim = 2L,
                                placement,
                                placement_missing,
                                rounds,
                                rounds_missing,
                                final_rounds,
                                final_rounds_missing,
                                num_init,
                                num_init_missing,
                                num_nbrs,
                                num_nbrs_missing,
                                r,
                                r_missing,
                                s,
                                s_missing,
                                repulsion_factor,
                                repulsion_factor_missing) {
  if (is.null(preset)) {
    return(list(
      placement = placement,
      rounds = rounds,
      final_rounds = final_rounds,
      num_init = num_init,
      num_nbrs = num_nbrs,
      r = r,
      s = s,
      repulsion_factor = repulsion_factor
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
  if (placement_missing) placement <- defaults$placement
  if (rounds_missing) rounds <- defaults$rounds
  if (final_rounds_missing) final_rounds <- defaults$final_rounds
  if (num_init_missing) num_init <- defaults$num_init
  if (num_nbrs_missing) num_nbrs <- defaults$num_nbrs
  if (r_missing) r <- defaults$r
  if (s_missing) s <- defaults$s
  if (repulsion_factor_missing) repulsion_factor <- defaults$repulsion_factor

  list(
    placement = placement,
    rounds = rounds,
    final_rounds = final_rounds,
    num_init = num_init,
    num_nbrs = num_nbrs,
    r = r,
    s = s,
    repulsion_factor = repulsion_factor
  )
}

grip.validate.tuning.inputs <- function(num_nbrs, r, s, repulsion_factor) {
  if (!is.numeric(num_nbrs) || length(num_nbrs) != 1L || !is.finite(num_nbrs)) {
    stop("num_nbrs must be a single finite numeric value")
  }
  if (abs(num_nbrs - round(num_nbrs)) > sqrt(.Machine$double.eps)) {
    stop("num_nbrs must be a positive integer")
  }
  num_nbrs <- as.integer(round(num_nbrs))
  if (is.na(num_nbrs) || num_nbrs <= 0L) {
    stop("num_nbrs must be a positive integer")
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

  if (!is.numeric(repulsion_factor) || length(repulsion_factor) != 1L ||
      !is.finite(repulsion_factor)) {
    stop("repulsion_factor must be a single finite numeric value")
  }
  repulsion_factor <- as.double(repulsion_factor)
  if (repulsion_factor < 0) {
    stop("repulsion_factor must be >= 0")
  }

  list(
    num_nbrs = num_nbrs,
    r = r,
    s = s,
    repulsion_factor = repulsion_factor
  )
}

grip.validate.globalrep.tuning.inputs <- function(num_nbrs,
                                                  r,
                                                  s,
                                                  repulsion_factor,
                                                  coarse_repulsion_factor,
                                                  coarse_repulsion_sample,
                                                  coarse_repulsion_exact_below,
                                                  final_anchor_factor = 0,
                                                  final_move_scale_after_first = 1,
                                                  level0_insertion_mode = "inherit",
                                                  level0_anchor_count = 3,
                                                  level0_local_kk_steps = 3) {
  tuning <- grip.validate.tuning.inputs(
    num_nbrs = num_nbrs,
    r = r,
    s = s,
    repulsion_factor = repulsion_factor
  )

  if (!is.numeric(coarse_repulsion_factor) ||
      length(coarse_repulsion_factor) != 1L ||
      !is.finite(coarse_repulsion_factor)) {
    stop("coarse_repulsion_factor must be a single finite numeric value")
  }
  coarse_repulsion_factor <- as.double(coarse_repulsion_factor)
  if (coarse_repulsion_factor < 0) {
    stop("coarse_repulsion_factor must be >= 0")
  }

  if (!is.numeric(coarse_repulsion_sample) ||
      length(coarse_repulsion_sample) != 1L ||
      !is.finite(coarse_repulsion_sample)) {
    stop("coarse_repulsion_sample must be a single finite numeric value")
  }
  if (abs(coarse_repulsion_sample - round(coarse_repulsion_sample)) >
      sqrt(.Machine$double.eps)) {
    stop("coarse_repulsion_sample must be a positive integer")
  }
  coarse_repulsion_sample <- as.integer(round(coarse_repulsion_sample))
  if (is.na(coarse_repulsion_sample) || coarse_repulsion_sample <= 0L) {
    stop("coarse_repulsion_sample must be a positive integer")
  }

  if (!is.numeric(coarse_repulsion_exact_below) ||
      length(coarse_repulsion_exact_below) != 1L ||
      !is.finite(coarse_repulsion_exact_below)) {
    stop("coarse_repulsion_exact_below must be a single finite numeric value")
  }
  if (abs(coarse_repulsion_exact_below - round(coarse_repulsion_exact_below)) >
      sqrt(.Machine$double.eps)) {
    stop("coarse_repulsion_exact_below must be a positive integer")
  }
  coarse_repulsion_exact_below <- as.integer(round(coarse_repulsion_exact_below))
  if (is.na(coarse_repulsion_exact_below) || coarse_repulsion_exact_below <= 0L) {
    stop("coarse_repulsion_exact_below must be a positive integer")
  }

  if (!is.numeric(final_anchor_factor) ||
      length(final_anchor_factor) != 1L ||
      !is.finite(final_anchor_factor)) {
    stop("final_anchor_factor must be a single finite numeric value")
  }
  final_anchor_factor <- as.double(final_anchor_factor)
  if (final_anchor_factor < 0) {
    stop("final_anchor_factor must be >= 0")
  }

  if (!is.numeric(final_move_scale_after_first) ||
      length(final_move_scale_after_first) != 1L ||
      !is.finite(final_move_scale_after_first)) {
    stop("final_move_scale_after_first must be a single finite numeric value")
  }
  final_move_scale_after_first <- as.double(final_move_scale_after_first)
  if (final_move_scale_after_first < 0 || final_move_scale_after_first > 1) {
    stop("final_move_scale_after_first must be in [0, 1]")
  }

  level0_insertion_mode <- match.arg(
    level0_insertion_mode,
    choices = c("inherit", "barycenter", "least_squares")
  )

  if (!is.numeric(level0_anchor_count) ||
      length(level0_anchor_count) != 1L ||
      !is.finite(level0_anchor_count)) {
    stop("level0_anchor_count must be a single finite numeric value")
  }
  if (abs(level0_anchor_count - round(level0_anchor_count)) >
      sqrt(.Machine$double.eps)) {
    stop("level0_anchor_count must be a positive integer")
  }
  level0_anchor_count <- as.integer(round(level0_anchor_count))
  if (is.na(level0_anchor_count) || level0_anchor_count <= 0L) {
    stop("level0_anchor_count must be a positive integer")
  }

  if (!is.numeric(level0_local_kk_steps) ||
      length(level0_local_kk_steps) != 1L ||
      !is.finite(level0_local_kk_steps)) {
    stop("level0_local_kk_steps must be a single finite numeric value")
  }
  if (abs(level0_local_kk_steps - round(level0_local_kk_steps)) >
      sqrt(.Machine$double.eps)) {
    stop("level0_local_kk_steps must be a non-negative integer")
  }
  level0_local_kk_steps <- as.integer(round(level0_local_kk_steps))
  if (is.na(level0_local_kk_steps) || level0_local_kk_steps < 0L) {
    stop("level0_local_kk_steps must be a non-negative integer")
  }

  c(
    tuning,
    list(
      coarse_repulsion_factor = coarse_repulsion_factor,
      coarse_repulsion_sample = coarse_repulsion_sample,
      coarse_repulsion_exact_below = coarse_repulsion_exact_below,
      final_anchor_factor = final_anchor_factor,
      final_move_scale_after_first = final_move_scale_after_first,
      level0_insertion_mode = level0_insertion_mode,
      level0_anchor_count = level0_anchor_count,
      level0_local_kk_steps = level0_local_kk_steps
    )
  )
}

grip.validate.layout.inputs <- function(edges = NULL,
                                        n = NULL,
                                        adj_list = NULL,
                                        weight_list = NULL,
                                        edge_weights = NULL,
                                        dim = 3,
                                        placement = "barycenter",
                                        seed = 6) {
  if (!is.null(adj_list)) {
    if (!is.list(adj_list)) stop("adj_list must be a list of integer vectors")
    if (is.null(n)) n <- length(adj_list)
    if (!is.numeric(n) || length(n) != 1L) {
      stop("n must be a single numeric value")
    }
    n <- as.integer(n)
    if (n <= 0L) stop("n must be positive")
    if (length(adj_list) != n) {
      stop("adj_list length must match n")
    }

    adj_list <- lapply(adj_list, function(v) {
      if (!is.numeric(v)) stop("adj_list entries must be numeric/integer vectors")
      v <- as.integer(v)
      if (any(!is.finite(v))) stop("adj_list entries must be finite")
      if (any(v <= 0L | v > n)) stop("adj_list must be 1-based and within [1, n]")
      v
    })

    if (!is.null(weight_list)) {
      if (!is.list(weight_list) || length(weight_list) != n) {
        stop("weight_list must be a list parallel to adj_list")
      }
      for (i in seq_len(n)) {
        wi <- weight_list[[i]]
        if (!is.numeric(wi)) {
          stop(sprintf("weight_list[[%d]] must be a numeric vector", i))
        }
        if (length(wi) != length(adj_list[[i]])) {
          stop(sprintf("weight_list[[%d]] must be parallel to adj_list[[%d]]", i, i))
        }
        bad <- which(!is.finite(wi) | wi <= 0)
        if (length(bad) > 0L) {
          j <- bad[[1L]]
          stop(sprintf(
            "weight_list must contain finite values > 0; first invalid at weight_list[[%d]][%d] = %s",
            i,
            j,
            format(wi[j], digits = 16)
          ))
        }
        weight_list[[i]] <- as.double(wi)
      }
    }
  } else {
    if (is.null(edges)) {
      stop("provide either edges or adj_list/weight_list")
    }
    if (!is.numeric(n) || length(n) != 1L) {
      stop("n must be a single numeric value")
    }
    n <- as.integer(n)
    if (n <= 0L) stop("n must be positive")
    if (!is.null(edge_weights)) {
      if (length(edge_weights) != nrow(edges)) {
        stop("edge_weights length must match number of edges")
      }
      if (!is.numeric(edge_weights)) {
        stop("edge_weights must be a numeric vector")
      }
      bad <- which(!is.finite(edge_weights) | edge_weights <= 0)
      if (length(bad) > 0L) {
        i <- bad[[1L]]
        stop(sprintf(
          "edge_weights must contain finite values > 0; first invalid at edge_weights[%d] = %s",
          i,
          format(edge_weights[i], digits = 16)
        ))
      }
    }

    converted <- grip.build.adj.from.edges(edges = edges, n = n, edge_weights = edge_weights)
    adj_list <- converted$adj_list
    weight_list <- converted$weight_list
  }

  if (!is.numeric(dim) || !(dim %in% c(2, 3))) {
    stop("dim must be 2 or 3")
  }
  dim <- as.integer(dim)
  if (!is.null(seed)) seed <- as.integer(seed)
  if (placement == "circle" && dim != 2) {
    warning("circle placement is only used for 2D; falling back to barycenter")
  }

  list(
    adj_list = adj_list,
    weight_list = weight_list,
    n = n,
    dim = dim,
    seed = seed
  )
}

#' Compute a GRIP layout with coarse global repulsion
#'
#' This compatibility entry point is an alias of \code{\link{grip.layout}()}.
#' It keeps the old global-repulsion name available during the API migration,
#' while using the same quality-first multiscale GRIP engine and adaptive
#' \code{final_rounds} schedule as the primary layout API.
#'
#' @inheritParams grip.layout
#' @param coarse_repulsion_factor Non-negative multiplier applied to the extra
#'   coarse-level active-set repulsion term. \code{0} disables that extra term;
#'   when the remaining tuning arguments also match
#'   \code{\link{grip.layout.legacy}()}, the result matches the legacy
#'   layout behavior.
#' @param coarse_repulsion_sample Positive integer sample size used to
#'   approximate active-set-wide repulsion on larger coarse levels.
#' @param coarse_repulsion_exact_below Positive integer threshold. When the
#'   active set size is at most this value, the coarse repulsion is computed
#'   exactly against all currently active vertices instead of being sampled.
#' @param final_anchor_factor Non-negative multiplier for an anchor term that
#'   pulls the final FR stage back toward the pre-final full-graph layout.
#'   `0` disables the anchor and preserves the current behavior.
#' @param final_move_scale_after_first Scalar in `[0, 1]` applied to the final
#'   FR displacement after the first finest-level round. Values below `1`
#'   damp later full-graph movement while keeping the first FR round unchanged.
#' @param final_mode Final full-graph refinement mode. \code{"fr"} keeps the
#'   current Fruchterman-Reingold-style final stage. \code{"kk_repulse"} uses a
#'   KK-style local distance-matching update with explicit active-set
#'   repulsion instead of the final FR phase.
#' @param level0_insertion_mode Level-0 insertion placement override used only
#'   when the finest filtration level is first populated. \code{"inherit"}
#'   keeps the current GRIP behavior. \code{"barycenter"} disables the 2D
#'   circle heuristic at level 0 and uses barycentric anchor placement.
#'   \code{"least_squares"} uses a multi-anchor least-squares distance fit at
#'   level 0 before any local micro-polish.
#' @param level0_anchor_count Positive integer number of already placed anchors
#'   to collect for level-0 insertion experiments. The legacy behavior uses 3.
#' @param level0_local_kk_steps Non-negative integer number of tiny local KK
#'   micro-polish steps applied immediately after each level-0 insertion. The
#'   legacy behavior uses 3.
#' @return A numeric matrix with `n` rows and `dim` columns.
#' @examples
#' edges <- edges.mesh(4, 4)
#' coords <- grip.layout.globalrep(edges, n = max(edges), dim = 2,
#'                                 rounds = 8, final_rounds = 8,
#'                                 num_init = 6, num_nbrs = 8,
#'                                 coarse_repulsion_factor = 0.2,
#'                                 coarse_repulsion_sample = 8,
#'                                 coarse_repulsion_exact_below = 32,
#'                                 seed = 1)
#' round(coords, 3)
#' @export
grip.layout.globalrep <- function(edges = NULL,
                                  n = NULL,
                                  adj_list = NULL,
                                  weight_list = NULL,
                                  edge_weights = NULL,
                                  dim = 3,
                                  placement = c("barycenter", "circle"),
                                  preset = NULL,
                                  rounds = 160,
                                  final_rounds = 384,
                                  num_init = 24,
                                  num_nbrs = 20,
                                  r = 0.03,
                                  s = 7.5,
                                  repulsion_factor = 2.5,
                                  coarse_repulsion_factor = 1.5,
                                  coarse_repulsion_sample = 16,
                                  coarse_repulsion_exact_below = 64,
                                  final_anchor_factor = 0,
                                  final_move_scale_after_first = 1,
                                  final_mode = c("fr", "kk_repulse"),
                                  level0_insertion_mode = c("inherit", "barycenter", "least_squares"),
                                  level0_anchor_count = 3,
                                  level0_local_kk_steps = 3,
                                  tinit_factor = 6,
                                  seed = 6,
                                  disconnected = c("components", "error")) {
  placement_missing <- missing(placement)
  rounds_missing <- missing(rounds)
  final_rounds_missing <- missing(final_rounds)
  num_init_missing <- missing(num_init)
  num_nbrs_missing <- missing(num_nbrs)
  r_missing <- missing(r)
  s_missing <- missing(s)
  repulsion_factor_missing <- missing(repulsion_factor)

  preset <- grip.normalize.preset(preset, fn = "grip.layout.globalrep")
  resolved <- grip.resolve.preset(
    preset = preset,
    dim = dim,
    placement = placement,
    placement_missing = placement_missing,
    rounds = rounds,
    rounds_missing = rounds_missing,
    final_rounds = final_rounds,
    final_rounds_missing = final_rounds_missing,
    num_init = num_init,
    num_init_missing = num_init_missing,
    num_nbrs = num_nbrs,
    num_nbrs_missing = num_nbrs_missing,
    r = r,
    r_missing = r_missing,
    s = s,
    s_missing = s_missing,
    repulsion_factor = repulsion_factor,
    repulsion_factor_missing = repulsion_factor_missing
  )
  placement <- resolved$placement
  rounds <- resolved$rounds
  final_rounds <- resolved$final_rounds
  num_init <- resolved$num_init
  num_nbrs <- resolved$num_nbrs
  r <- resolved$r
  s <- resolved$s
  repulsion_factor <- resolved$repulsion_factor
  placement <- match.arg(placement)
  final_mode <- match.arg(final_mode)
  level0_insertion_mode <- match.arg(level0_insertion_mode)
  disconnected <- match.arg(disconnected)

  validated <- grip.validate.layout.inputs(
    edges = edges,
    n = n,
    adj_list = adj_list,
    weight_list = weight_list,
    edge_weights = edge_weights,
    dim = dim,
    placement = placement,
    seed = seed
  )
  adj_list <- validated$adj_list
  weight_list <- validated$weight_list
  n <- validated$n
  dim <- validated$dim
  seed <- validated$seed
  if (is.null(preset) && final_rounds_missing) {
    final_rounds <- grip.globalrep.default.final_rounds(n)
  }
  tuning <- grip.validate.globalrep.tuning.inputs(
    num_nbrs = num_nbrs,
    r = r,
    s = s,
    repulsion_factor = repulsion_factor,
    coarse_repulsion_factor = coarse_repulsion_factor,
    coarse_repulsion_sample = coarse_repulsion_sample,
    coarse_repulsion_exact_below = coarse_repulsion_exact_below,
    final_anchor_factor = final_anchor_factor,
    final_move_scale_after_first = final_move_scale_after_first,
    level0_insertion_mode = level0_insertion_mode,
    level0_anchor_count = level0_anchor_count,
    level0_local_kk_steps = level0_local_kk_steps
  )
  num_nbrs <- tuning$num_nbrs
  r <- tuning$r
  s <- tuning$s
  repulsion_factor <- tuning$repulsion_factor
  coarse_repulsion_factor <- tuning$coarse_repulsion_factor
  coarse_repulsion_sample <- tuning$coarse_repulsion_sample
  coarse_repulsion_exact_below <- tuning$coarse_repulsion_exact_below
  final_anchor_factor <- tuning$final_anchor_factor
  final_move_scale_after_first <- tuning$final_move_scale_after_first
  level0_insertion_mode <- tuning$level0_insertion_mode
  level0_anchor_count <- tuning$level0_anchor_count
  level0_local_kk_steps <- tuning$level0_local_kk_steps

  layout.adj <- function(adj_list, weight_list, n) {
    grip_layout_globalrep_adj_cpp(
      adj_list = adj_list,
      weight_list = weight_list,
      n = n,
      dim = dim,
      placement = placement,
      rounds = as.integer(rounds),
      final_rounds = as.integer(final_rounds),
      num_init = as.integer(num_init),
      num_nbrs = num_nbrs,
      r = r,
      s = s,
      repulsion_factor = repulsion_factor,
      coarse_repulsion_factor = coarse_repulsion_factor,
      coarse_repulsion_sample = coarse_repulsion_sample,
      coarse_repulsion_exact_below = coarse_repulsion_exact_below,
      final_anchor_factor = final_anchor_factor,
      final_move_scale_after_first = final_move_scale_after_first,
      level0_insertion_mode = level0_insertion_mode,
      level0_anchor_count = level0_anchor_count,
      level0_local_kk_steps = level0_local_kk_steps,
      final_mode = final_mode,
      tinit_factor = as.integer(tinit_factor),
      seed = seed
    )
  }

  comp <- grip.connected.components(adj_list = adj_list, n = n)
  n.comp <- length(unique(comp))

  if (n.comp == 1L) {
    return(layout.adj(adj_list = adj_list, weight_list = weight_list, n = n))
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
      adj_list = adj_list,
      weight_list = weight_list,
      vertices = rows,
      n = n
    )
    layouts[[k]] <- layout.adj(
      adj_list = sub$adj_list,
      weight_list = sub$weight_list,
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
#' tapers \code{final_rounds} on larger graphs.
#'
#' @param edges Two-column integer matrix of edges (1-based vertex ids).
#' @param n Number of vertices.
#' @param adj_list Adjacency list (1-based) for undirected graphs.
#' @param weight_list Optional parallel list of edge weights (edge lengths).
#'   If NULL, all edges are treated as weight 1. All weights must be finite
#'   and strictly positive.
#' @param edge_weights Optional vector of edge weights for \code{edges}. All
#'   weights must be finite and strictly positive.
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
#' @param final_rounds Final rounds for refinement.
#' @param num_init Number of initial vertices in the coarsest level.
#' @param num_nbrs Maximum number of graph-distance neighbors retained for local
#'   refinement at each filtration level.
#' @param r Main local temperature adaptation rate in \code{[0, 1]}.
#' @param s Non-negative boost factor applied when successive displacements have
#'   a consistent direction.
#' @param repulsion_factor Non-negative multiplier applied to GRIP's
#'   finest-level repulsive force scale.
#' @param coarse_repulsion_factor Non-negative multiplier applied to the extra
#'   coarse-level active-set repulsion term. \code{0} disables that extra term.
#' @param coarse_repulsion_sample Positive integer sample size used to
#'   approximate active-set-wide repulsion on larger coarse levels.
#' @param coarse_repulsion_exact_below Positive integer threshold. When the
#'   active set size is at most this value, the coarse repulsion is computed
#'   exactly against all currently active vertices instead of being sampled.
#' @param final_anchor_factor Non-negative multiplier for an anchor term that
#'   pulls the final FR stage back toward the pre-final full-graph layout.
#'   `0` disables the anchor and preserves the current behavior.
#' @param final_move_scale_after_first Scalar in `[0, 1]` applied to the final
#'   FR displacement after the first finest-level round. Values below `1`
#'   damp later full-graph movement while keeping the first FR round unchanged.
#' @param final_mode Final full-graph refinement mode. \code{"fr"} keeps the
#'   current Fruchterman-Reingold-style final stage. \code{"kk_repulse"} uses a
#'   KK-style local distance-matching update with explicit active-set
#'   repulsion instead of the final FR phase.
#' @param level0_insertion_mode Level-0 insertion placement override used only
#'   when the finest filtration level is first populated. \code{"inherit"}
#'   keeps the current GRIP behavior. \code{"barycenter"} disables the 2D
#'   circle heuristic at level 0 and uses barycentric anchor placement.
#'   \code{"least_squares"} uses a multi-anchor least-squares distance fit at
#'   level 0 before any local micro-polish.
#' @param level0_anchor_count Positive integer number of already placed anchors
#'   to collect for level-0 insertion experiments. The legacy behavior uses 3.
#' @param level0_local_kk_steps Non-negative integer number of tiny local KK
#'   micro-polish steps applied immediately after each level-0 insertion. The
#'   legacy behavior uses 3.
#' @param tinit_factor Initial temperature factor.
#' @param seed Optional RNG seed for reproducibility. If NULL, uses current time.
#' @param disconnected How to handle disconnected graphs:
#'   \code{"components"} (default) lays out each connected component separately
#'   and packs them into one coordinate matrix; \code{"error"} stops with an
#'   error.
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
#' coords <- grip.layout(edges, n = max(edges), dim = 2,
#'                       coarse_repulsion_factor = 0.2,
#'                       coarse_repulsion_sample = 8,
#'                       coarse_repulsion_exact_below = 32,
#'                       seed = 1)
#' round(coords, 3)
#' @export
grip.layout <- function(edges = NULL,
                        n = NULL,
                        adj_list = NULL,
                        weight_list = NULL,
                        edge_weights = NULL,
                        dim = 3,
                        placement = c("barycenter", "circle"),
                        preset = NULL,
                        rounds = 160,
                        final_rounds = 384,
                        num_init = 24,
                        num_nbrs = 20,
                        r = 0.03,
                        s = 7.5,
                        repulsion_factor = 2.5,
                        coarse_repulsion_factor = 1.5,
                        coarse_repulsion_sample = 16,
                        coarse_repulsion_exact_below = 64,
                        final_anchor_factor = 0,
                        final_move_scale_after_first = 1,
                        final_mode = c("fr", "kk_repulse"),
                        level0_insertion_mode = c("inherit", "barycenter", "least_squares"),
                        level0_anchor_count = 3,
                        level0_local_kk_steps = 3,
                        tinit_factor = 6,
                        seed = 6,
                        disconnected = c("components", "error")) {
  grip.forward_call(grip.layout.globalrep, match.call(expand.dots = FALSE), env = parent.frame())
}

#' Compute the legacy GRIP layout
#'
#' This function preserves the original local-force wrapper and historical
#' default values that were previously exposed as \code{\link{grip.layout}()}.
#' Use it for backwards-compatible comparisons or when you explicitly want the
#' pre-global-repulsion behavior.
#'
#' @param edges Two-column integer matrix of edges (1-based vertex ids).
#' @param n Number of vertices.
#' @param adj_list Adjacency list (1-based) for undirected graphs.
#' @param weight_list Optional parallel list of edge weights (edge lengths).
#'   If NULL, all edges are treated as weight 1. All weights must be finite
#'   and strictly positive.
#' @param edge_weights Optional vector of edge weights for \code{edges}. All
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
#' @param final_rounds Final rounds for refinement.
#' @param num_init Number of initial vertices in the coarsest level.
#' @param num_nbrs Maximum number of graph-distance neighbors retained for local
#'   refinement at each filtration level.
#' @param r Main local temperature adaptation rate in \code{[0, 1]}.
#' @param s Non-negative boost factor applied when successive displacements have
#'   a consistent direction.
#' @param repulsion_factor Non-negative multiplier applied to GRIP's
#'   finest-level repulsive force scale. \code{1} keeps the historical
#'   repulsion strength; \code{0} disables that repulsive term.
#' @param tinit_factor Initial temperature factor.
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
#' coords <- grip.layout.legacy(edges, n = 6, dim = 2,
#'                              placement = "barycenter",
#'                              rounds = 5, final_rounds = 5,
#'                              num_init = 3, num_nbrs = 4,
#'                              seed = 1)
#' round(coords, 3)
#' @export
grip.layout.legacy <- function(edges = NULL,
                               n = NULL,
                               adj_list = NULL,
                               weight_list = NULL,
                               edge_weights = NULL,
                               dim = 3,
                               placement = c("barycenter", "circle"),
                               preset = NULL,
                               rounds = 20,
                               final_rounds = 25,
                               num_init = 36,
                               num_nbrs = 10,
                               r = 0.15,
                               s = 3.0,
                               repulsion_factor = 1.0,
                               tinit_factor = 6,
                               seed = 6,
                               disconnected = c("components", "error")) {
  placement_missing <- missing(placement)
  rounds_missing <- missing(rounds)
  final_rounds_missing <- missing(final_rounds)
  num_init_missing <- missing(num_init)
  num_nbrs_missing <- missing(num_nbrs)
  r_missing <- missing(r)
  s_missing <- missing(s)
  repulsion_factor_missing <- missing(repulsion_factor)

  preset <- grip.normalize.preset(preset, fn = "grip.layout.legacy")
  resolved <- grip.resolve.preset(
    preset = preset,
    dim = dim,
    placement = placement,
    placement_missing = placement_missing,
    rounds = rounds,
    rounds_missing = rounds_missing,
    final_rounds = final_rounds,
    final_rounds_missing = final_rounds_missing,
    num_init = num_init,
    num_init_missing = num_init_missing,
    num_nbrs = num_nbrs,
    num_nbrs_missing = num_nbrs_missing,
    r = r,
    r_missing = r_missing,
    s = s,
    s_missing = s_missing,
    repulsion_factor = repulsion_factor,
    repulsion_factor_missing = repulsion_factor_missing
  )
  placement <- resolved$placement
  rounds <- resolved$rounds
  final_rounds <- resolved$final_rounds
  num_init <- resolved$num_init
  num_nbrs <- resolved$num_nbrs
  r <- resolved$r
  s <- resolved$s
  repulsion_factor <- resolved$repulsion_factor
  placement <- match.arg(placement)
  disconnected <- match.arg(disconnected)

  validated <- grip.validate.layout.inputs(
    edges = edges,
    n = n,
    adj_list = adj_list,
    weight_list = weight_list,
    edge_weights = edge_weights,
    dim = dim,
    placement = placement,
    seed = seed
  )
  adj_list <- validated$adj_list
  weight_list <- validated$weight_list
  n <- validated$n
  dim <- validated$dim
  seed <- validated$seed
  tuning <- grip.validate.tuning.inputs(
    num_nbrs = num_nbrs,
    r = r,
    s = s,
    repulsion_factor = repulsion_factor
  )
  num_nbrs <- tuning$num_nbrs
  r <- tuning$r
  s <- tuning$s
  repulsion_factor <- tuning$repulsion_factor

  layout.adj <- function(adj_list, weight_list, n) {
    grip_layout_adj_cpp(adj_list = adj_list,
                        weight_list = weight_list,
                        n = n,
                        dim = dim,
                        placement = placement,
                        rounds = as.integer(rounds),
                        final_rounds = as.integer(final_rounds),
                        num_init = as.integer(num_init),
                        num_nbrs = num_nbrs,
                        r = r,
                        s = s,
                        repulsion_factor = repulsion_factor,
                        tinit_factor = as.integer(tinit_factor),
                        seed = seed)
  }

  comp <- grip.connected.components(adj_list = adj_list, n = n)
  n.comp <- length(unique(comp))

  if (n.comp == 1L) {
    return(layout.adj(adj_list = adj_list, weight_list = weight_list, n = n))
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
      adj_list = adj_list,
      weight_list = weight_list,
      vertices = rows,
      n = n
    )
    layouts[[k]] <- layout.adj(
      adj_list = sub$adj_list,
      weight_list = sub$weight_list,
      n = length(rows)
    )
  }
  grip.pack.component.layouts(layouts = layouts, comp = comp, n = n, dim = dim)
}

#' Compute a GRIP layout trace
#'
#' This traces the primary \code{\link{grip.layout}()} engine, including the
#' coarse-level global-repulsion term used by the quality-first default layout.
#' For backwards-compatible traces of the historical local-force wrapper, use
#' \code{\link{grip.layout.trace.legacy}()}.
#'
#' @inheritParams grip.layout
#' @param trace Snapshot granularity. \code{"round"} records the coarsest
#'   initialization, each level start, every \code{trace.every} completed rounds,
#'   and the final layout. \code{"level"} records the coarsest initialization,
#'   every \code{trace.every}th level start, and the final layout.
#' @param trace.every Positive integer thinning factor for recorded rounds or
#'   levels. Initial and final snapshots are always included.
#' @param diagnostics Optional per-frame diagnostic mode. \code{"none"} skips
#'   extra scoring, \code{"light"} appends lightweight shape diagnostics, and
#'   \code{"full"} also computes sampled stress on each traced frame.
#' @param target_coords Optional numeric target coordinate matrix used to append
#'   per-frame Procrustes RMSE diagnostics. It must have `n` rows and `dim`
#'   columns.
#' @param diagnostic_sample_size_nonedge Positive integer sample size used for
#'   per-frame non-edge separation diagnostics when \code{diagnostics != "none"}.
#' @param diagnostic_sample_size_stress Positive integer sample size used for
#'   per-frame sampled stress when \code{diagnostics = "full"}.
#' @param diagnostic_nonedge_seed RNG seed base used for per-frame non-edge
#'   separation diagnostics.
#' @param diagnostic_stress_seed RNG seed base used for per-frame sampled stress
#'   diagnostics.
#' @return A list with \code{final}, \code{frames}, \code{meta}, \code{trace},
#'   \code{trace.every}, and optionally \code{diagnostics}. \code{final} is the
#'   final coordinate matrix.
#'   \code{frames} is a list of coordinate matrices with \code{NA} rows for
#'   vertices that have not yet been introduced by GRIP. \code{meta} is a data
#'   frame describing each frame with columns \code{frame}, \code{phase},
#'   \code{level_index}, \code{misf_level}, \code{round_in_level}, and
#'   \code{active_vertices}. When diagnostics are requested, \code{diagnostics}
#'   is a data frame parallel to \code{meta} that appends per-frame quality
#'   metrics such as \code{edge.length.cv}, \code{sampled.nonedge.sep.ratio},
#'   and optional \code{procrustes.rmse}.
#' @examples
#' edges <- cbind(1:5, 2:6)
#' tr <- grip.layout.trace(edges, n = 6, dim = 2,
#'                         placement = "barycenter",
#'                         rounds = 3, final_rounds = 2,
#'                         num_init = 3, num_nbrs = 4,
#'                         trace = "level",
#'                         trace.every = 1,
#'                         diagnostics = "light",
#'                         seed = 1)
#' tr$diagnostics
#' @export
grip.layout.trace <- function(edges = NULL,
                              n = NULL,
                              adj_list = NULL,
                              weight_list = NULL,
                              edge_weights = NULL,
                              dim = 3,
                              placement = c("barycenter", "circle"),
                              preset = NULL,
                              rounds = 160,
                              final_rounds = 384,
                              num_init = 24,
                              num_nbrs = 20,
                              r = 0.03,
                              s = 7.5,
                              repulsion_factor = 2.5,
                              coarse_repulsion_factor = 1.5,
                              coarse_repulsion_sample = 16,
                              coarse_repulsion_exact_below = 64,
                              final_anchor_factor = 0,
                              final_move_scale_after_first = 1,
                              final_mode = c("fr", "kk_repulse"),
                              level0_insertion_mode = c("inherit", "barycenter", "least_squares"),
                              level0_anchor_count = 3,
                              level0_local_kk_steps = 3,
                              tinit_factor = 6,
                              seed = 6,
                              trace = c("round", "level"),
                              trace.every = 1,
                              diagnostics = c("none", "light", "full"),
                              target_coords = NULL,
                              diagnostic_sample_size_nonedge = 1000L,
                              diagnostic_sample_size_stress = 500L,
                              diagnostic_nonedge_seed = 1L,
                              diagnostic_stress_seed = 1L) {
  placement_missing <- missing(placement)
  rounds_missing <- missing(rounds)
  final_rounds_missing <- missing(final_rounds)
  num_init_missing <- missing(num_init)
  num_nbrs_missing <- missing(num_nbrs)
  r_missing <- missing(r)
  s_missing <- missing(s)
  repulsion_factor_missing <- missing(repulsion_factor)

  preset <- grip.normalize.preset(preset, fn = "grip.layout.trace")
  resolved <- grip.resolve.preset(
    preset = preset,
    dim = dim,
    placement = placement,
    placement_missing = placement_missing,
    rounds = rounds,
    rounds_missing = rounds_missing,
    final_rounds = final_rounds,
    final_rounds_missing = final_rounds_missing,
    num_init = num_init,
    num_init_missing = num_init_missing,
    num_nbrs = num_nbrs,
    num_nbrs_missing = num_nbrs_missing,
    r = r,
    r_missing = r_missing,
    s = s,
    s_missing = s_missing,
    repulsion_factor = repulsion_factor,
    repulsion_factor_missing = repulsion_factor_missing
  )
  placement <- resolved$placement
  rounds <- resolved$rounds
  final_rounds <- resolved$final_rounds
  num_init <- resolved$num_init
  num_nbrs <- resolved$num_nbrs
  r <- resolved$r
  s <- resolved$s
  repulsion_factor <- resolved$repulsion_factor
  placement <- match.arg(placement)
  final_mode <- match.arg(final_mode)
  level0_insertion_mode <- match.arg(level0_insertion_mode)
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
    adj_list = adj_list,
    weight_list = weight_list,
    edge_weights = edge_weights,
    dim = dim,
    placement = placement,
    seed = seed
  )
  adj_list <- validated$adj_list
  weight_list <- validated$weight_list
  n <- validated$n
  dim <- validated$dim
  seed <- validated$seed
  if (is.null(preset) && final_rounds_missing) {
    final_rounds <- grip.globalrep.default.final_rounds(n)
  }
  tuning <- grip.validate.globalrep.tuning.inputs(
    num_nbrs = num_nbrs,
    r = r,
    s = s,
    repulsion_factor = repulsion_factor,
    coarse_repulsion_factor = coarse_repulsion_factor,
    coarse_repulsion_sample = coarse_repulsion_sample,
    coarse_repulsion_exact_below = coarse_repulsion_exact_below,
    final_anchor_factor = final_anchor_factor,
    final_move_scale_after_first = final_move_scale_after_first,
    level0_insertion_mode = level0_insertion_mode,
    level0_anchor_count = level0_anchor_count,
    level0_local_kk_steps = level0_local_kk_steps
  )
  num_nbrs <- tuning$num_nbrs
  r <- tuning$r
  s <- tuning$s
  repulsion_factor <- tuning$repulsion_factor
  coarse_repulsion_factor <- tuning$coarse_repulsion_factor
  coarse_repulsion_sample <- tuning$coarse_repulsion_sample
  coarse_repulsion_exact_below <- tuning$coarse_repulsion_exact_below
  final_anchor_factor <- tuning$final_anchor_factor
  final_move_scale_after_first <- tuning$final_move_scale_after_first
  level0_insertion_mode <- tuning$level0_insertion_mode
  level0_anchor_count <- tuning$level0_anchor_count
  level0_local_kk_steps <- tuning$level0_local_kk_steps

  comp <- grip.connected.components(adj_list = adj_list, n = n)
  n.comp <- length(unique(comp))
  if (n.comp != 1L) {
    stop(sprintf(
      "grip.layout.trace() currently supports only connected graphs; input graph has %d connected components.",
      n.comp
    ))
  }

  out <- grip_layout_globalrep_trace_adj_cpp(
    adj_list = adj_list,
    weight_list = weight_list,
    n = n,
    dim = dim,
    placement = placement,
    rounds = as.integer(rounds),
    final_rounds = as.integer(final_rounds),
    num_init = as.integer(num_init),
    num_nbrs = num_nbrs,
    r = r,
    s = s,
    repulsion_factor = repulsion_factor,
    coarse_repulsion_factor = coarse_repulsion_factor,
    coarse_repulsion_sample = coarse_repulsion_sample,
    coarse_repulsion_exact_below = coarse_repulsion_exact_below,
    final_anchor_factor = final_anchor_factor,
    final_move_scale_after_first = final_move_scale_after_first,
    level0_insertion_mode = level0_insertion_mode,
    level0_anchor_count = level0_anchor_count,
    level0_local_kk_steps = level0_local_kk_steps,
    final_mode = final_mode,
    tinit_factor = as.integer(tinit_factor),
    seed = seed,
    trace = trace,
    trace_every = trace.every
  )
  out$trace <- trace
  out$trace.every <- trace.every
  out$diagnostics <- grip.trace.compute.diagnostics(
    frames = out$frames,
    meta = out$meta,
    adj.list = adj_list,
    weight.list = weight_list,
    diagnostics = diagnostics,
    target.coords = target_coords,
    sample.size.nonedge = diagnostic_sample_size_nonedge,
    sample.size.stress = diagnostic_sample_size_stress,
    nonedge.seed = diagnostic_nonedge_seed,
    stress.seed = diagnostic_stress_seed
  )
  class(out) <- c("grip_layout_trace", class(out))
  out
}

#' Compute a trace for the legacy GRIP layout
#'
#' @inheritParams grip.layout.legacy
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
#' tr <- grip.layout.trace.legacy(edges, n = 6, dim = 2,
#'                                placement = "barycenter",
#'                                rounds = 3, final_rounds = 2,
#'                                num_init = 3, num_nbrs = 4,
#'                                trace = "level",
#'                                trace.every = 1,
#'                                seed = 1)
#' tr$meta
#' @export
grip.layout.trace.legacy <- function(edges = NULL,
                                     n = NULL,
                                     adj_list = NULL,
                                     weight_list = NULL,
                                     edge_weights = NULL,
                                     dim = 3,
                                     placement = c("barycenter", "circle"),
                                     preset = NULL,
                                     rounds = 20,
                                     final_rounds = 25,
                                     num_init = 36,
                                     num_nbrs = 10,
                                     r = 0.15,
                                     s = 3.0,
                                     repulsion_factor = 1.0,
                                     tinit_factor = 6,
                                     seed = 6,
                                     trace = c("round", "level"),
                                     trace.every = 1) {
  placement_missing <- missing(placement)
  rounds_missing <- missing(rounds)
  final_rounds_missing <- missing(final_rounds)
  num_init_missing <- missing(num_init)
  num_nbrs_missing <- missing(num_nbrs)
  r_missing <- missing(r)
  s_missing <- missing(s)
  repulsion_factor_missing <- missing(repulsion_factor)

  preset <- grip.normalize.preset(preset, fn = "grip.layout.trace.legacy")
  resolved <- grip.resolve.preset(
    preset = preset,
    dim = dim,
    placement = placement,
    placement_missing = placement_missing,
    rounds = rounds,
    rounds_missing = rounds_missing,
    final_rounds = final_rounds,
    final_rounds_missing = final_rounds_missing,
    num_init = num_init,
    num_init_missing = num_init_missing,
    num_nbrs = num_nbrs,
    num_nbrs_missing = num_nbrs_missing,
    r = r,
    r_missing = r_missing,
    s = s,
    s_missing = s_missing,
    repulsion_factor = repulsion_factor,
    repulsion_factor_missing = repulsion_factor_missing
  )
  placement <- resolved$placement
  rounds <- resolved$rounds
  final_rounds <- resolved$final_rounds
  num_init <- resolved$num_init
  num_nbrs <- resolved$num_nbrs
  r <- resolved$r
  s <- resolved$s
  repulsion_factor <- resolved$repulsion_factor
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
    adj_list = adj_list,
    weight_list = weight_list,
    edge_weights = edge_weights,
    dim = dim,
    placement = placement,
    seed = seed
  )
  adj_list <- validated$adj_list
  weight_list <- validated$weight_list
  n <- validated$n
  dim <- validated$dim
  seed <- validated$seed
  tuning <- grip.validate.tuning.inputs(
    num_nbrs = num_nbrs,
    r = r,
    s = s,
    repulsion_factor = repulsion_factor
  )
  num_nbrs <- tuning$num_nbrs
  r <- tuning$r
  s <- tuning$s
  repulsion_factor <- tuning$repulsion_factor

  comp <- grip.connected.components(adj_list = adj_list, n = n)
  n.comp <- length(unique(comp))
  if (n.comp != 1L) {
    stop(sprintf(
      "grip.layout.trace.legacy() currently supports only connected graphs; input graph has %d connected components.",
      n.comp
    ))
  }

  out <- grip_layout_trace_adj_cpp(
    adj_list = adj_list,
    weight_list = weight_list,
    n = n,
    dim = dim,
    placement = placement,
    rounds = as.integer(rounds),
    final_rounds = as.integer(final_rounds),
    num_init = as.integer(num_init),
    num_nbrs = num_nbrs,
    r = r,
    s = s,
    repulsion_factor = repulsion_factor,
    tinit_factor = as.integer(tinit_factor),
    seed = seed,
    trace = trace,
    trace_every = trace.every
  )
  out$trace <- trace
  out$trace.every <- trace.every
  class(out) <- c("grip_layout_trace", class(out))
  out
}
