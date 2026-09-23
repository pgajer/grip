# Validate graph identity before integer coercion, inference, or engine dispatch.
# These checks deliberately do not canonicalize loops, duplicates, or weights.
grip.validate.vertex.count <- function(n) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
      n < 1 || n > .Machine$integer.max || n != trunc(n)) {
    stop("n must be a single finite positive integer within the R integer range",
         call. = FALSE)
  }
  as.integer(n)
}

grip.validate.vertex.ids <- function(x, name, n = NULL) {
  if (!is.numeric(x)) {
    stop(name, " must contain numeric integer vertex ids", call. = FALSE)
  }
  upper <- if (is.null(n)) .Machine$integer.max else min(n, .Machine$integer.max)
  bad <- which(!is.finite(x) | x < 1 | x > upper | x != trunc(x))
  if (length(bad)) {
    i <- bad[[1L]]
    index <- if (is.matrix(x)) paste(arrayInd(i, dim(x)), collapse = ", ") else i
    stop(sprintf("%s[%s] must be a finite integer vertex id within [1, %s]; received %s",
                 name, index, if (is.null(n)) .Machine$integer.max else n,
                 format(x[i], digits = 16)), call. = FALSE)
  }
  invisible(x)
}

grip.validate.graph.arguments <- function(edges = NULL, n = NULL,
                                          adj.list = NULL, weight.list = NULL,
                                          edge.weights = NULL, prepared = NULL) {
  if (!is.null(prepared)) {
    if (any(!vapply(list(edges, adj.list, weight.list, edge.weights), is.null, logical(1)))) {
      stop("prepared cannot be combined with edges, adj.list, weight.list, or edge.weights; supply one graph representation",
           call. = FALSE)
    }
    if (!is.null(n)) {
      n <- grip.validate.vertex.count(n)
      if (!is.null(prepared$n) && !identical(n, grip.validate.vertex.count(prepared$n))) {
        stop("n must match the graph size stored in prepared", call. = FALSE)
      }
    }
    return(invisible(NULL))
  }
  if (!is.null(edges) && !is.null(adj.list)) {
    stop("provide either edges or adj_list, not both", call. = FALSE)
  }
  if (!is.null(weight.list) && is.null(adj.list)) {
    stop("weight.list requires adj.list; use edge.weights with edges", call. = FALSE)
  }
  if (!is.null(edge.weights) && is.null(edges)) {
    stop("edge.weights requires edges; use weight.list with adj.list", call. = FALSE)
  }
  if (!is.null(n)) n <- grip.validate.vertex.count(n)
  if (!is.null(edges)) {
    edges <- as.matrix(edges)
    if (ncol(edges) != 2L) {
      stop("edges must be a two-column integer matrix of 1-based vertex ids", call. = FALSE)
    }
    grip.validate.vertex.ids(edges, "edges", n)
    if (any(edges[, 1L] == edges[, 2L])) {
      stop("self-loops are not supported; remove edges whose endpoints are equal", call. = FALSE)
    }
  }
  if (!is.null(adj.list)) {
    if (!is.list(adj.list)) stop("adj.list must be a list of integer vectors", call. = FALSE)
    if (!is.null(n) && length(adj.list) != n) stop("adj.list length must match n", call. = FALSE)
    for (i in seq_along(adj.list)) {
      # NULL is an established representation of an isolated vertex.
      if (!is.null(adj.list[[i]])) {
        grip.validate.vertex.ids(adj.list[[i]], sprintf("adj_list[[%d]]", i),
                                 if (is.null(n)) length(adj.list) else n)
      }
    }
    grip.validate.undirected.adjacency(adj.list, weight.list)
  }
  invisible(NULL)
}

grip.resolve.graph.n <- function(n, edges, adj.list) {
  if (is.null(n)) {
    n <- if (!is.null(adj.list)) length(adj.list) else if (length(edges)) max(edges) else NULL
  }
  grip.validate.vertex.count(n)
}

# Compare sorted directed entries, retaining parallel-edge multiplicity.
# This checks a contract; it never reorders or canonicalizes the caller's graph.
grip.validate.undirected.adjacency <- function(adj.list, weight.list = NULL) {
  n <- length(adj.list)
  if (!is.null(weight.list)) {
    if (!is.list(weight.list) || length(weight.list) != n) {
      stop("weight.list must be a list parallel to adj.list", call. = FALSE)
    }
    for (i in seq_len(n)) {
      w <- weight.list[[i]]
      if (length(w) != length(adj.list[[i]]) ||
          (length(w) && (!is.numeric(w) || any(!is.finite(w) | w <= 0)))) {
        stop(sprintf("weight.list[[%d]] must contain positive finite lengths parallel to adj.list[[%d]]", i, i), call. = FALSE)
      }
    }
  }
  from <- rep.int(seq_len(n), lengths(adj.list))
  to <- as.integer(unlist(adj.list, use.names = FALSE))
  if (any(from == to)) stop("self-loops are not supported; remove each vertex from its own adjacency", call. = FALSE)
  if (!length(from)) return(invisible(NULL))
  weights <- if (is.null(weight.list)) rep(1, length(from)) else as.double(unlist(weight.list))
  forward <- order(from, to, weights)
  reverse <- order(to, from, weights)
  if (!identical(from[forward], to[reverse]) || !identical(to[forward], from[reverse])) {
    stop("adj.list must be reciprocal: every u -> v needs v -> u with the same multiplicity", call. = FALSE)
  }
  if (!identical(weights[forward], weights[reverse])) {
    stop("weight.list must give matching lengths for reciprocal undirected edges", call. = FALSE)
  }
  invisible(NULL)
}

grip.validate.layout.inputs <- function(edges = NULL,
                                        n = NULL,
                                        adj.list = NULL,
                                        weight.list = NULL,
                                        edge.weights = NULL,
                                        dim = 3,
                                        placement = "barycenter",
                                        seed = 6) {
  grip.validate.graph.arguments(edges, n, adj.list, weight.list, edge.weights)
  if (!is.null(adj.list)) {
    if (!is.list(adj.list)) stop("adj.list must be a list of integer vectors")
    if (is.null(n)) n <- length(adj.list)
    n <- grip.validate.vertex.count(n)
    if (length(adj.list) != n) {
      stop("adj.list length must match n")
    }

    adj.list <- lapply(adj.list, function(v) {
      if (is.null(v)) return(integer())
      if (!is.numeric(v)) stop("adj.list entries must be numeric/integer vectors")
      v <- as.integer(v)
      if (any(!is.finite(v))) stop("adj.list entries must be finite")
      if (any(v <= 0L | v > n)) stop("adj.list must be 1-based and within [1, n]")
      v
    })

    if (!is.null(weight.list)) {
      if (!is.list(weight.list) || length(weight.list) != n) {
        stop("weight.list must be a list parallel to adj.list")
      }
      for (i in seq_len(n)) {
        wi <- weight.list[[i]]
        if (!is.numeric(wi)) {
          stop(sprintf("weight.list[[%d]] must be a numeric vector", i))
        }
        if (length(wi) != length(adj.list[[i]])) {
          stop(sprintf("weight.list[[%d]] must be parallel to adj.list[[%d]]", i, i))
        }
        bad <- which(!is.finite(wi) | wi <= 0)
        if (length(bad) > 0L) {
          j <- bad[[1L]]
          stop(sprintf(
            "weight.list must contain finite values > 0; first invalid at weight.list[[%d]][%d] = %s",
            i,
            j,
            format(wi[j], digits = 16)
          ))
        }
        weight.list[[i]] <- as.double(wi)
      }
    }
  } else {
    if (is.null(edges)) {
      stop("provide either edges or adj_list/weight_list")
    }
    n <- grip.validate.vertex.count(n)
    if (!is.null(edge.weights)) {
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
    }

    converted <- grip.build.adj.from.edges(edges = edges, n = n, edge.weights = edge.weights)
    adj.list <- converted$adj_list
    weight.list <- converted$weight_list
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
    adj_list = adj.list,
    weight_list = weight.list,
    n = n,
    dim = dim,
    seed = seed
  )
}
