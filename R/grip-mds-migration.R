#' Migrating classical and metric stress MDS workflows
#'
#' In grip through version 0.2.0, `metric.mds()` used `stats::cmdscale()`.
#' From version 0.2.0.9000, that algorithm is named [classical.mds()], and
#' [metric.mds()] minimizes raw distance stress, using native SGD by default.
#' This is an intentional behavioral change rather than a deprecated alias.
#'
#' @section Preserve an existing analysis:
#' Replace `metric.mds(...)` with `classical.mds(...)`, and explicit
#' `init = "metric_mds"` with `init = "classical_mds"`. Both [edge.kk()]
#' and [kernel.gram.gkk()] default to classical initialization, preserving
#' their previous default algorithm without requiring smacof.
#' `add` and `eig` are classical-scaling arguments and are accepted only by
#' `classical.mds()`. Existing saved coordinates, figure labels, and cache
#' fields named `metric_mds` from older versions remain classical results.
#' Renaming a call does not require recomputing those fits.
#'
#' @section Request stress minimization:
#' Use `metric.mds(...)` or `edge.kk(init = "metric_mds", ...)` to use the
#' default native SGD backend. For control over multiple starts and budgets,
#' call `metric.mds()` first and pass `coords = fit$coords` to `edge.kk()`.
#' Check `metadata$starts`, `metadata$termination`, and the independently
#' calculated stress diagnostics. Both backends are local optimizers; an
#' iteration limit or a finite result does not certify an optimum.
#'
#' @section Preserve a previous SMACOF analysis:
#' Install the optional smacof package and explicitly request
#' `metric.mds(backend = "smacof", ...)`. The `eps` tolerance belongs to that
#' backend and must not be supplied to SGD. To retain a SMACOF initializer in
#' a refinement workflow, compute it explicitly and pass its coordinates to
#' `edge.kk()` or `kernel.gram.gkk()`. No automatic backend fallback is made.
#'
#' @section Native SGD backend:
#' `metric.mds(backend = "sgd", max.iter = 30, ...)` explicitly selects the
#' default backend with a short trial budget, without requiring smacof.
#' SGD uses the same uniform pair weights and shortest-path targets. Its controls
#' remain provisional calibration choices, and `max.iter` counts complete pair
#' passes. Schedule completion is reported as `iteration_limit`, not convergence.
#' Inspect per-start records and `metadata$sgd` checkpoint histories. Selecting
#' `init = "metric_mds"` in refinement methods follows the default SGD backend;
#' their default `init = "classical_mds"` is unchanged. A different optimizer
#' cannot repair inaccurate graph distances.
#'
#' @section Reproducibility:
#' Raw stress is measured against the original graph-distance units after
#' restoring scale from the backend. Target-normalized raw stress and
#' scale-profiled Stress-1 have the same optimal shapes when scale is free
#' and pair weights match, but literal Stress-1 at the returned raw-stress
#' scale is a different diagnostic value. Record the algorithm, objective,
#' backend/package versions, graph identity, scale policy, and start settings
#' in new result manifests. Recompute downstream refinement when changing its
#' initializer, and use new cache identities rather than replacing historical
#' classical results in place.
#'
#' @md
#' @name grip-mds-migration
NULL
