#!/usr/bin/env Rscript
lib <- Sys.getenv("GRIP_RJOURNAL_PACKAGE_LIBRARY")
if(nzchar(lib)) .libPaths(c(lib,.libPaths()))
args <- commandArgs(trailingOnly=TRUE); stopifnot(length(args)>=1)
out <- normalizePath(args[1]); reps <- if(length(args)>1) as.integer(strsplit(args[2],",")[[1]]) else 1:5
kmax <- if(length(args)>2) as.integer(args[3]) else 20L
for(r in reps) {
  dest <- file.path(out,sprintf("graphs-%02d.rds",r))
  previous<-if(file.exists(dest)) readRDS(dest) else NULL
  x <- readRDS(file.path(out,sprintf("cloud-%02d.rds",r)))$coords
  graphs <- if(is.null(previous)) list() else previous$graphs
  started <- proc.time()[[3]]
  for(k in 3:kmax) {
    if(as.character(k) %in% names(graphs)) next
    tick <- proc.time()[[3]]
    g <- dgraphs::create.sknn.graph(x,k=k,neighbor.method="exact",
      connect.components=TRUE,connect.method="component.mst",prune.method="none",
      prune.edges=FALSE,edge.weight="distance")
    ig <- igraph::graph_from_edgelist(g$edge_matrix,directed=FALSE)
    D <- igraph::distances(ig,weights=g$edge_weight)
    stopifnot(igraph::is_connected(ig),all(is.finite(D)),
      g$n_mst_edges_added==g$n_components_before-1L,
      max(abs(g$edge_weight-sqrt(rowSums((x[g$edge_matrix[,1],]-x[g$edge_matrix[,2],])^2))))<1e-12)
    graphs[[as.character(k)]] <- list(k=k,edges=g$edge_matrix,weights=g$edge_weight,
      distances=D,mst_edges=g$mst_edge_matrix,mst_weights=g$mst_edge_weight,
      components_before=g$n_components_before,bridges=g$n_mst_edges_added,
      elapsed=proc.time()[[3]]-tick)
    cat(sprintf("cloud=%d k=%d edges=%d bridges=%d seconds=%.2f\n",r,k,nrow(g$edge_matrix),g$n_mst_edges_added,proc.time()[[3]]-tick));flush.console()
  }
  saveRDS(list(replicate=r,graphs=graphs,elapsed=proc.time()[[3]]-started+
    if(is.null(previous)) 0 else previous$elapsed,session=sessionInfo()),dest)
}
