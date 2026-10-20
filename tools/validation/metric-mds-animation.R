#!/usr/bin/env Rscript
# Numerical QA for the executable vignette; no browser or ivue required.
# Optional arguments are the original paraboloid and saddle checkpoint RDS files.
pkgload::load_all(".", recompile = FALSE, quiet = TRUE)
source("inst/scripts/metric-mds-animation.R")
reference.files <- commandArgs(trailingOnly = TRUE)
X <- surface.sample()
Y <- X
Y[, 3] <- Y[, 1]^2 - Y[, 2]^2
cases <- list(paraboloid = surface.trace(X), saddle = surface.trace(Y))
stopifnot(identical(X[, 1:2], Y[, 1:2]))
for (name in names(cases)) {
  case <- cases[[name]]
  control <- case$control
  control$trace <- NULL
  plain <- surface.epoch(
    distance.matrix = as.matrix(dist(case$truth)), init = case$initial,
    dim = 3L, backend = "sgd", n.init = 1L, max.iter = case$epochs,
    seed = case$seed, pair.weights = "inverse_squared",
    diagnostics = FALSE, sgd.control = control)
  stopifnot(identical(plain$coords, case$fit$coords),
            identical(plain$continuation$state$native,
                      case$fit$continuation$state$native))
  target <- as.numeric(dist(case$truth))
  endpoint <- which(case$steps$pair_updates %% case$pair.count == 0)
  stopifnot(length(endpoint) == case$epochs + 1L,
            identical(case$steps$pair_updates[endpoint],
                      as.numeric(0:case$epochs * case$pair.count)),
            100*tail(case$errors, 1L) <= 1e-12,
            case$epochs == if (name == "paraboloid") 6L else 5L)
  # Every selected state and displayed error must remain consistent.
  for (i in seq_along(case$raw)) {
    stopifnot(abs(surface.relative.stress(case$raw[[i]], target) -
                    case$errors[i]) < 1e-14,
              abs(surface.relative.stress(case$aligned[[i]], target) -
                    case$errors[i]) < 1e-12)
  }
  if (length(reference.files)) {
    stopifnot(length(reference.files) == 2L)
    old <- readRDS(reference.files[[match(name, names(cases))]])
    stopifnot(is.list(old$three), !is.null(old$three$fit),
              identical(unname(case$truth), unname(old$truth)),
              identical(case$fit$coords, old$three$fit$coords),
              identical(case$fit$continuation$state$native,
                        old$three$fit$continuation$state$native))
    stopifnot(max(abs(case$raw[[1L]] - old$three$initial)) < 1e-14)
    for (i in seq_along(endpoint)[-1L])
      stopifnot(identical(case$raw[[endpoint[i]]], old$three$frames[[i]]))
  }
  cat(name, ":", case$epochs, "epochs;", length(case$raw),
      "saved states; 100 S =", 100*tail(case$errors, 1L),
      "; traced/resumed and untraced/uninterrupted states and RNG identical.\n")
}
cat("Matched x,y, all saved errors, rigid alignment and epoch endpoints passed.\n")
