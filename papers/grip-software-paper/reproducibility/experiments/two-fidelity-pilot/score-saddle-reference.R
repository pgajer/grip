# Reference scoring for an existing single cloud and its fitted layouts.
# Does not refit layouts or change the sampled-saddle experiment's scoring convention.
saddle.reference.mesh <- function(X, triangles, C, subdivisions = 2L) {
  uv <- X[, 1:2, drop = FALSE]
  for (level in seq_len(subdivisions)) {
    edges <- rbind(triangles[, 1:2], triangles[, 2:3], triangles[, c(3, 1)])
    keys <- paste(pmin(edges[,1],edges[,2]),pmax(edges[,1],edges[,2]),sep=":")
    unique.edges <- edges[!duplicated(keys), , drop=FALSE]
    midpoint <- (uv[unique.edges[,1],,drop=FALSE]+uv[unique.edges[,2],,drop=FALSE])/2
    ids <- nrow(uv) + match(keys, unique(keys))
    nf <- nrow(triangles); ab <- ids[seq_len(nf)]; bc <- ids[nf+seq_len(nf)]; ca <- ids[2*nf+seq_len(nf)]
    triangles <- rbind(cbind(triangles[,1],ab,ca),cbind(ab,triangles[,2],bc),
                      cbind(ca,bc,triangles[,3]),cbind(ab,bc,ca))
    uv <- rbind(uv,midpoint)
  }
  list(coords=cbind(x=uv[,1],y=uv[,2],z=C*(uv[,1]^2-uv[,2]^2)),triangles=triangles)
}

score.saddle.reference <- function(X, graphs, mds, mds.edge.kk, C,
                                    sample_sizes=c(2000L,8000L), seed=1901L,
                                    subdivisions=2L, progress=TRUE) {
  stopifnot(ncol(X)==3L, all(names(graphs) %in% names(mds)),
            all(names(graphs) %in% names(mds.edge.kk)),
            max(abs(X[,3]-C*(X[,1]^2-X[,2]^2)))<1e-10)
  triangles <- geometry::delaunayn(X[,1:2])
  reference <- saddle.reference.mesh(X, triangles, C, subdivisions)
  records <- list(); counter <- 0L
  for (key in names(graphs)) {
    if (progress) message("Scoring k = ",key)
    g <- graphs[[key]]
    prepared <- grip::prepare.geodesic.kk(g$edge_matrix,n=nrow(X),
        edge_weights=g$edge_weight,tie_mode="single")
    candidates <- list("Original saddle"=X,"Metric MDS"=mds[[key]],
                       "MDS + edge-KK"=mds.edge.kk[[key]])
    for (method in names(candidates)) {
      Z <- candidates[[method]]
      score <- grip::score.gmds(Z,prepared=prepared,scale_mode="profiled",edge_length_epsilon=0)
      pairs <- prepared$pair_matrix
      chord <- sqrt(rowSums((Z[pairs[,1],,drop=FALSE]-Z[pairs[,2],,drop=FALSE])^2))
      dg <- prepared$pair_graph_distance
      a <- sum(chord*dg)/sum(dg^2)
      stress1 <- sqrt(sum((chord-a*dg)^2)/sum(chord^2))
      # Independent conversion of target-normalized chord residual to Stress-1.
      stopifnot(abs(stress1-score$metric.chord.stress*a*sqrt(sum(dg^2)/sum(chord^2)))<1e-10)
      for (alignment in c("rigid","similarity")) {
        fit <- grip::score.coordinates(Z,X,alignment=alignment)
        for (ns in sample_sizes) {
          surf <- grip::score.surface(fit$coords,triangles,reference$coords,reference$triangles,
                                      sample_size=ns,seed=seed)
          counter <- counter+1L
          records[[counter]] <- data.frame(n=nrow(X),k=as.integer(key),method=method,
            alignment=alignment,path_rel=score$gmds.stress,edge_rel=score$edge.rel.rmse,
            stress1=stress1,coordinate_rmse=fit$rmse,coordinate_relative_rmse=fit$relative_rmse,
            alignment_scale=fit$scale,surface_rms=surf$rms,surface_mc_se=surf$rms_mc_se,
            surface_forward_rms=surf$forward_rms,surface_reverse_rms=surf$reverse_rms,
            surface_area=surf$area,reference_area=surf$reference_area,
            zero_area_faces=surf$zero_area_faces,sample_size=ns,surface_seed=seed,
            subdivisions=subdivisions)
        }
      }
    }
  }
  list(scores=do.call(rbind,records),triangles=triangles,reference=reference)
}
