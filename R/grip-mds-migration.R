#' Migrating classical and metric stress MDS workflows
#'
#' In grip through version 0.2.0, `metric.mds()` used `stats::cmdscale()`.
#' From version 0.2.0.9000, that algorithm is named [classical.mds()], and
#' [metric.mds()] is a ratio-SMACOF wrapper that minimizes raw distance stress.
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
#' Install the optional smacof package and use `metric.mds(...)` or
#' `edge.kk(init = "metric_mds", ...)`. For control over multiple starts and
#' tolerances, call `metric.mds()` first and pass its coordinates to `edge.kk()`.
#' Check `metadata$starts`, `metadata$termination`, and the independently
#' calculated stress diagnostics. SMACOF is a local optimizer; an iteration
#' limit or a finite result does not certify an optimum. No fallback to
#' classical MDS is made when smacof is unavailable.
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
