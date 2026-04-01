grip.rectangular.grid.bending.stencils <- function(h, w) {
  h <- as.integer(h)
  w <- as.integer(w)
  grip.validate.scalar(h, "h", lower = 1)
  grip.validate.scalar(w, "w", lower = 1)
  if (h < 3L && w < 3L) {
    return(matrix(integer(), ncol = 3L))
  }

  index <- matrix(seq_len(h * w), nrow = h, ncol = w, byrow = TRUE)
  stencils <- vector("list", max(0L, h * max(0L, w - 2L) + w * max(0L, h - 2L)))
  k <- 1L
  if (w >= 3L) {
    for (r in seq_len(h)) {
      for (c in 2:(w - 1L)) {
        stencils[[k]] <- c(index[r, c - 1L], index[r, c], index[r, c + 1L])
        k <- k + 1L
      }
    }
  }
  if (h >= 3L) {
    for (r in 2:(h - 1L)) {
      for (c in seq_len(w)) {
        stencils[[k]] <- c(index[r - 1L, c], index[r, c], index[r + 1L, c])
        k <- k + 1L
      }
    }
  }
  out <- do.call(rbind, stencils[seq_len(k - 1L)])
  storage.mode(out) <- "integer"
  out
}

grip.validate.bending.stencils <- function(bending_stencils,
                                           n = NULL) {
  if (is.null(bending_stencils)) {
    return(NULL)
  }
  if (!is.matrix(bending_stencils) || ncol(bending_stencils) != 3L) {
    stop("bending_stencils must be NULL or a three-column integer matrix")
  }
  out <- matrix(as.integer(bending_stencils), ncol = 3L)
  if (nrow(out) == 0L) {
    return(out)
  }
  if (any(!is.finite(out))) {
    stop("bending_stencils must contain only finite vertex ids")
  }
  if (!is.null(n) && any(out < 1L | out > as.integer(n))) {
    stop("bending_stencils must contain 1-based vertex ids within [1, nrow(coords)]")
  }
  out
}

grip.flatten.bending.stencils.zero.based <- function(bending_stencils) {
  stencils <- grip.validate.bending.stencils(bending_stencils)
  if (is.null(stencils) || nrow(stencils) == 0L) {
    return(list(
      flat_bend_a = integer(),
      flat_bend_b = integer(),
      flat_bend_c = integer()
    ))
  }
  list(
    flat_bend_a = as.integer(stencils[, 1L] - 1L),
    flat_bend_b = as.integer(stencils[, 2L] - 1L),
    flat_bend_c = as.integer(stencils[, 3L] - 1L)
  )
}
