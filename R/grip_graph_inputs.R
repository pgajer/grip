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
                                          adj_list = NULL, weight_list = NULL,
                                          edge_weights = NULL, prepared = NULL) {
  if (!is.null(prepared)) {
    if (any(!vapply(list(edges, adj_list, weight_list, edge_weights), is.null, logical(1)))) {
      stop("prepared cannot be combined with edges, adj_list, weight_list, or edge_weights; supply one graph representation",
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
  if (!is.null(edges) && !is.null(adj_list)) {
    stop("provide either edges or adj_list, not both", call. = FALSE)
  }
  if (!is.null(weight_list) && is.null(adj_list)) {
    stop("weight_list requires adj_list; use edge_weights with edges", call. = FALSE)
  }
  if (!is.null(edge_weights) && is.null(edges)) {
    stop("edge_weights requires edges; use weight_list with adj_list", call. = FALSE)
  }
  if (!is.null(n)) n <- grip.validate.vertex.count(n)
  if (!is.null(edges)) {
    edges <- as.matrix(edges)
    if (ncol(edges) != 2L) {
      stop("edges must be a two-column integer matrix of 1-based vertex ids", call. = FALSE)
    }
    grip.validate.vertex.ids(edges, "edges", n)
  }
  if (!is.null(adj_list)) {
    if (!is.list(adj_list)) stop("adj_list must be a list of integer vectors", call. = FALSE)
    for (i in seq_along(adj_list)) {
      # NULL is an established representation of an isolated vertex.
      if (!is.null(adj_list[[i]])) {
        grip.validate.vertex.ids(adj_list[[i]], sprintf("adj_list[[%d]]", i),
                                 if (is.null(n)) length(adj_list) else n)
      }
    }
  }
  invisible(NULL)
}

grip.resolve.graph.n <- function(n, edges, adj_list) {
  if (is.null(n)) {
    n <- if (!is.null(adj_list)) length(adj_list) else if (length(edges)) max(edges) else NULL
  }
  grip.validate.vertex.count(n)
}
