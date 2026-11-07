## ----setup, include=FALSE-----------------------------------------------------
library(grip)
options(rgl.useNULL = TRUE)
knitr::opts_chunk$set(fig.retina=1, dpi=72, collapse = TRUE, comment = "#>", fig.width = 11,
                    fig.height = 6, out.width = "100%")
source(system.file("scripts", "metric-mds-animation.R", package = "grip"),
       local = TRUE)

## ----sample-------------------------------------------------------------------
set.seed(20261020)
n <- 250L
u <- runif(n)
theta <- runif(n, 0, 2*pi)
r <- sqrt(((1 + u*(5^(3/2) - 1))^(2/3) - 1)/4)
paraboloid <- cbind(x = r*cos(theta), y = r*sin(theta), z = r^2)
saddle <- paraboloid
saddle[, "z"] <- saddle[, "x"]^2 - saddle[, "y"]^2
stopifnot(identical(paraboloid[, 1:2], saddle[, 1:2]))

## ----initial------------------------------------------------------------------
truth <- paraboloid
target <- as.numeric(dist(truth))
relative.stress <- function(X) mean((as.numeric(dist(X))/target - 1)^2)

set.seed(20261021)
initial <- matrix(rnorm(n*3L), n, 3L)
initial <- sweep(initial, 2, colMeans(initial))
ratio <- as.numeric(dist(initial))/target
initial <- initial * sum(ratio)/sum(ratio^2)

## ----one-epoch----------------------------------------------------------------
control <- list(
  scheduler = "hybrid", learning.rate = 0.5, final.rate = 0.01,
  switch.ratio = 0.4, schedule.epochs = 500L,
  checkpoint.every = 1L, keep.state = TRUE,
  trace = list(every = 32L, max.bytes = 64*1024^2)
)
fit <- surface.epoch(
  distance.matrix = as.matrix(dist(truth)), init = initial, dim = 3L,
  backend = "sgd", n.init = 1L, max.iter = 1L, seed = 20261022L,
  pair.weights = "inverse_squared", diagnostics = FALSE,
  sgd.control = control
)
head(fit$update_trace$steps)

## ----inspect------------------------------------------------------------------
recorded <- fit$update_trace
first <- recorded$coordinates[[1L]]
last <- tail(recorded$coordinates, 1L)[[1L]]
c(initial = relative.stress(first), after_one_epoch = relative.stress(last))

## ----continue-----------------------------------------------------------------
next.fit <- surface.epoch(resume = fit, max.iter = 1L)
stopifnot(identical(next.fit$update_trace$coordinates[[1L]], last))
tail(next.fit$update_trace$steps, 1L)

## ----fit-surfaces-------------------------------------------------------------
cases <- list(
  paraboloid = surface.trace(paraboloid),
  saddle = surface.trace(saddle)
)

## ----align--------------------------------------------------------------------
align <- function(X, truth) {
  center <- colMeans(truth)
  centered <- sweep(X, 2, colMeans(X))
  decomposition <- svd(crossprod(centered, sweep(truth, 2, center)))
  sweep(centered %*% decomposition$u %*% t(decomposition$v), 2, center, "+")
}
aligned <- align(last, truth)
stopifnot(max(abs(as.numeric(dist(aligned)) - as.numeric(dist(last)))) < 1e-12)

## ----animation, echo=FALSE----------------------------------------------------
surface.animation(cases)

## ----results, echo=FALSE------------------------------------------------------
knitr::kable(do.call(rbind, lapply(names(cases), function(name) {
  case <- cases[[name]]
  data.frame(Surface = name, Epochs = case$epochs,
             `Saved frames` = length(case$aligned),
             `Final 100 S` = sprintf("%.6g", 100*tail(case$errors, 1L)),
             `RMS relative error (%)` =
               sprintf("%.6g", 100*sqrt(tail(case$errors, 1L))),
             Stop = case$stop.reason, check.names = FALSE)
})), row.names = FALSE)

## ----every-update, eval=FALSE-------------------------------------------------
# # Run sequentially: allow room for R objects as well as native recording.
# cases <- list(
#   paraboloid = surface.trace(paraboloid, every = 1L, max.bytes = 256*1024^2),
#   saddle = surface.trace(saddle, every = 1L, max.bytes = 256*1024^2)
# )
# surface.animation(cases)

## ----helper-files-------------------------------------------------------------
basename(system.file("scripts", c("metric-mds-animation.R",
                                 "metric-mds-animation.js"), package = "grip"))

## ----save, eval=FALSE---------------------------------------------------------
# source(system.file("scripts", "metric-mds-animation.R", package = "grip"))
# X <- surface.sample()
# Y <- X
# Y[, "z"] <- Y[, "x"]^2 - Y[, "y"]^2
# cases <- list(paraboloid = surface.trace(X), saddle = surface.trace(Y))
# options(grip.vignette.interactive = TRUE)
# widget <- surface.animation(cases)
# htmlwidgets::saveWidget(widget, "metric-mds-surfaces.html", selfcontained = TRUE)

