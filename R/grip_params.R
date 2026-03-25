#' Extract layout parameters from a comparison summary row
#'
#' Converts one row of the summary data frame returned by
#' \code{\link{grip.compare.layouts}} into a named list of parameters suitable
#' for passing to \code{\link{grip.layout}} or back into a
#' \code{grip.compare.layouts(candidates = ...)} call.
#'
#' @param summary.row A one-row data frame (or the first row is used) from the
#'   \code{$summary} element of \code{grip.compare.layouts()}.
#' @return A named list with elements \code{placement}, \code{rounds},
#'   \code{final_rounds}, \code{num_init}, \code{num_nbrs}, \code{r},
#'   \code{s}, \code{repulsion_factor}, \code{tinit_factor}, and optionally
#'   \code{preset} (included only when the summary row records a non-empty
#'   preset).
#' @examples
#' \dontrun{
#' edges <- edges.mesh(6, 6)
#' cmp <- grip.compare.layouts(edges, n = 36, dim = 2,
#'                             candidates = c("default", "mesh"),
#'                             seeds = 1)
#' params <- grip.params.from.summary(cmp$summary[1, ])
#' coords <- grip.layout(edges, n = 36, dim = 2,
#'                        placement      = params$placement,
#'                        rounds         = params$rounds,
#'                        final_rounds   = params$final_rounds,
#'                        num_init       = params$num_init,
#'                        num_nbrs       = params$num_nbrs,
#'                        r              = params$r,
#'                        s              = params$s,
#'                        repulsion_factor = params$repulsion_factor,
#'                        seed = 42)
#' }
#' @export
grip.params.from.summary <- function(summary.row) {
  if (!is.data.frame(summary.row) || nrow(summary.row) == 0L) {
    stop("summary.row must be a non-empty data frame")
  }
  row <- summary.row[1L, , drop = FALSE]
  out <- list(
    placement        = row$placement[[1L]],
    rounds           = row$rounds[[1L]],
    final_rounds     = row$final.rounds[[1L]],
    num_init         = row$num.init[[1L]],
    num_nbrs         = row$num.nbrs[[1L]],
    r                = row$r[[1L]],
    s                = row$s[[1L]],
    repulsion_factor = row$repulsion.factor[[1L]],
    tinit_factor     = row$tinit.factor[[1L]]
  )
  if (!is.na(row$preset[[1L]]) && nzchar(row$preset[[1L]])) {
    out$preset <- row$preset[[1L]]
  }
  out
}
