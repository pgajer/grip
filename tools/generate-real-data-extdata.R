#!/usr/bin/env Rscript

if (!requireNamespace("igraph", quietly = TRUE)) {
  stop("Install 'igraph' to regenerate the real-data extdata files.")
}

out.dir <- file.path("inst", "extdata")
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

write_edges_csv <- function(graph.name, file.name) {
  g <- igraph::make_graph(graph.name)
  edges <- igraph::as_edgelist(g, names = FALSE)
  edges <- as.data.frame(edges)
  names(edges) <- c("source", "target")
  utils::write.csv(edges, file.path(out.dir, file.name), row.names = FALSE)
}

write_edges_csv("Zachary", "karate-club-edges.csv")
write_edges_csv("Krackhardt_Kite", "krackhardt-kite-edges.csv")

cat("Wrote extdata files to", out.dir, "\n")
