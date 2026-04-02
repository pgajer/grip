\name{grip.build.misf}
\alias{grip.build.misf}
\title{Build a maximal independent set filtration (MISF) for a graph}
\usage{
grip.build.misf(
  edges = NULL,
  n = NULL,
  adj_list = NULL,
  weight_list = NULL,
  edge_weights = NULL,
  num_init = 24L,
  num_nbrs = 20L,
  seed = 6L
)
}
\arguments{
  \item{edges}{Two-column integer matrix of undirected edges (1-based vertex ids).}
  \item{n}{Number of vertices.}
  \item{adj_list}{Optional adjacency list (1-based integer vectors).}
  \item{weight_list}{Optional positive edge-weight list parallel to \code{adj_list}.}
  \item{edge_weights}{Optional positive vector parallel to \code{edges}.}
  \item{num_init}{Target top-level active-set size used by the current GRIP MISF builder.}
  \item{num_nbrs}{Retained local neighborhood budget per MISF level, matching the current GRIP refinement schedule metadata.}
  \item{seed}{Optional integer seed passed to the current GRIP graph RNG.}
}
\value{
An object of class \code{"grip_misf"} containing the nested MISF levels,
vertex depths, the internal MISF order, per-level sizes, and the retained
per-level neighborhood schedule.
}
\description{
\code{grip.build.misf()} exposes the maximal independent set filtration already
constructed internally by the GRIP layout engine. The result is graph-first:
you can supply either an edge list plus \code{n}, or an adjacency/weight-list
pair.
}
\details{
The current implementation returns the same MISF structure used by GRIP's
multiscale layout core. Edge weights are accepted and passed through the
underlying graph object, but the current MISF construction itself is driven by
the graph topology and hop-distance thresholds rather than weighted
shortest-path distances.
}
\examples{
edges <- edges.mesh(4, 4)
misf <- grip.build.misf(edges = edges, n = 16, num_init = 6, seed = 1)
misf$misf_size
misf$levels[[1L]]
}
