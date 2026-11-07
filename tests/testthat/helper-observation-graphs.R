# Observation-based sparse preparation is optional and needs an unreleased
# dgraphs API. CRAN installations can still test the remaining grip features.
skip_if_no_observation_graphs <- function() {
  skip_if_not_installed("dgraphs")
  if (!"search.data" %in% names(formals(dgraphs::create.sknn.graph))) {
    skip("dgraphs with search.data and streaming component repair is required")
  }
}
