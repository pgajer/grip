build_test_sierpinski_carpet <- function(level) {
  side <- 3L^level
  grid <- expand.grid(x = 0:(side - 1L), y = 0:(side - 1L))

  keep_cell <- function(x, y) {
    while (x > 0L || y > 0L) {
      if ((x %% 3L) == 1L && (y %% 3L) == 1L) {
        return(FALSE)
      }
      x <- x %/% 3L
      y <- y %/% 3L
    }
    TRUE
  }

  keep <- mapply(keep_cell, grid$x, grid$y)
  cells <- grid[keep, , drop = FALSE]
  coords <- cbind(
    x = cells$x + 0.5,
    y = (side - 1L - cells$y) + 0.5
  )

  id_map <- matrix(0L, nrow = side, ncol = side)
  for (i in seq_len(nrow(cells))) {
    id_map[cells$x[i] + 1L, cells$y[i] + 1L] <- i
  }

  edges <- list()
  for (i in seq_len(nrow(cells))) {
    x <- cells$x[i]
    y <- cells$y[i]
    if (x + 1L < side) {
      nbr <- id_map[x + 2L, y + 1L]
      if (nbr > 0L) edges[[length(edges) + 1L]] <- c(i, nbr)
    }
    if (y + 1L < side) {
      nbr <- id_map[x + 1L, y + 2L]
      if (nbr > 0L) edges[[length(edges) + 1L]] <- c(i, nbr)
    }
  }

  list(edges = do.call(rbind, edges), coords = coords)
}

test_that("grip.geometry.diagnostics returns strong carpet diagnostics on the canonical layout", {
  built <- build_test_sierpinski_carpet(2L)
  diag <- grip.geometry.diagnostics(
    coords = built$coords,
    target.coords = built$coords,
    edges = built$edges,
    family = "sierpinski.carpet",
    sample.size.symmetry = 128L,
    sample.size.wedges = 512L,
    rng.seed = 1L
  )

  expect_s3_class(diag, "data.frame")
  expect_equal(nrow(diag), 1L)
  expect_lt(diag$procrustes.rmse[[1L]], 1e-10)
  expect_gt(diag$global.symmetry.score[[1L]], 0.95)
  expect_lt(diag$local.angle.deviation[[1L]], 1e-10)
  expect_gt(diag$edge.axis.concentration[[1L]], 0.99)
  expect_lt(diag$boundary.waviness[[1L]], 1e-10)
  expect_lt(diag$corridor.waviness[[1L]], 1e-10)
  expect_lt(diag$hole.center.error[[1L]], 1e-10)
})

test_that("grip.geometry.diagnostics detects geometric degradation on a perturbed carpet", {
  built <- build_test_sierpinski_carpet(2L)
  perturbed <- built$coords
  set.seed(11)
  perturbed <- perturbed + matrix(rnorm(length(perturbed), sd = 0.15), ncol = 2L)

  baseline <- grip.geometry.diagnostics(
    coords = built$coords,
    target.coords = built$coords,
    edges = built$edges,
    family = "sierpinski.carpet",
    sample.size.symmetry = 128L,
    sample.size.wedges = 512L,
    rng.seed = 1L
  )
  degraded <- grip.geometry.diagnostics(
    coords = perturbed,
    target.coords = built$coords,
    edges = built$edges,
    family = "sierpinski.carpet",
    sample.size.symmetry = 128L,
    sample.size.wedges = 512L,
    rng.seed = 1L
  )

  expect_lt(baseline$global.symmetry.score[[1L]], 1.000001)
  expect_gt(degraded$procrustes.rmse[[1L]], baseline$procrustes.rmse[[1L]])
  expect_lt(degraded$global.symmetry.score[[1L]], baseline$global.symmetry.score[[1L]])
  expect_gt(degraded$local.angle.deviation[[1L]], baseline$local.angle.deviation[[1L]])
  expect_lt(degraded$edge.axis.concentration[[1L]], baseline$edge.axis.concentration[[1L]])
  expect_gt(degraded$boundary.waviness[[1L]], baseline$boundary.waviness[[1L]])
  expect_gt(degraded$corridor.waviness[[1L]], baseline$corridor.waviness[[1L]])
  expect_gt(degraded$hole.center.error[[1L]], baseline$hole.center.error[[1L]])
})

test_that("grip.geometry.diagnostics leaves carpet-only measures NA for non-carpet graphs", {
  coords <- cbind(
    x = c(-1, 0, 1, -1, 0, 1),
    y = c(0, 0, 0, 1, 1, 1)
  )
  edges <- edges.mesh(2, 3)

  diag <- grip.geometry.diagnostics(
    coords = coords,
    target.coords = coords,
    edges = edges,
    family = "mesh",
    sample.size.symmetry = 64L,
    sample.size.wedges = 128L,
    rng.seed = 1L
  )

  expect_true(is.finite(diag$global.symmetry.score[[1L]]))
  expect_true(is.finite(diag$edge.axis.concentration[[1L]]))
  expect_true(is.na(diag$boundary.waviness[[1L]]))
  expect_true(is.na(diag$corridor.waviness[[1L]]))
  expect_true(is.na(diag$hole.center.error[[1L]]))
})
