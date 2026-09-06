# All-route check adapted from the independent September 6, 2026 audit.
# Explicit optional validation against original saved caches; not a build dependency.
args <- commandArgs(TRUE)
stopifnot(length(args)==2)
root <- normalizePath(args[1]); dest <- args[2]
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
ans <- list(); pos <- 0L
sources <- character()
for (surface in c("saddle", "paraboloid")) for (rep in 1:3) {
  key <- paste0(surface, "-disk-rep", rep, "-n240-r64")
  inputs <- file.path(root, "output/mds-edge-kk-radius/inputs")
  if (!file.exists(file.path(inputs, paste0(key,"-truth.csv")))) {
    key <- paste0(surface, "-base_disk-rep", rep, "-n240-r64")
  }
  truth <- as.matrix(read.csv(file.path(inputs,paste0(key,"-truth.csv")),header=FALSE))
  ref <- as.matrix(read.csv(file.path(inputs,paste0(key,"-distance.csv")),header=FALSE))
  xc <- scale(truth, scale=FALSE)
  for (regime in c("ambient", "geodesic")) {
    fitdir <- file.path(root,"output/mds-edge-kk-radius/fits-v2",key)
    result <- readRDS(file.path(fitdir,paste0(regime,"-k32-result.rds")))$value
    p <- readRDS(file.path(fitdir,paste0(regime,"-k32-prepared.rds")))$value$p
    sources <- c(sources,
      file.path(inputs,paste0(key,"-truth.csv")),file.path(inputs,paste0(key,"-distance.csv")),
      file.path(fitdir,paste0(regime,"-k32-result.rds")),file.path(fitdir,paste0(regime,"-k32-prepared.rds")),
      file.path(root,"tools/experiments/mds-edge-kk-radius/summary/coordinates",paste0(key,".rds")))
    for (method in c("stress", "stress_fixed_primary")) {
      z <- result$candidates[[method]]
      portable <- readRDS(file.path(root,"tools/experiments/mds-edge-kk-radius/summary/coordinates",paste0(key,".rds")))
      matches <- Filter(function(v) v$scores$regime[1]==regime && v$scores$k[1]==32 && method %in% names(v$candidates),portable)
      stopifnot(length(matches)==1, max(abs(z-matches[[1]]$candidates[[method]]))<1e-12)
      # Independently sum every recorded route in ordinary R, in physical units.
      edge_lengths <- sqrt(rowSums((z[p$edges[,1],]-z[p$edges[,2],])^2))
      input_edges <- p$edge_targets*result$scale_unit
      edge_scale <- sum(edge_lengths*input_edges)/sum(input_edges^2)
      flat <- sqrt(rowSums((z[p$flat_edge_u+1L,]-z[p$flat_edge_v+1L,])^2))
      offsets <- p$flat_pair_edge_offsets
      stopifnot(length(offsets)-1L==28680L)
      paths <- vapply(seq_len(length(offsets)-1L),function(i) sum(flat[seq.int(offsets[i]+1L,offsets[i+1L])]),numeric(1))
      targets <- p$pair_graph_distance*result$scale_unit
      path_scale <- sum(paths*targets)/sum(targets^2)
      path_rel <- sqrt(sum((paths-path_scale*targets)^2)/sum((path_scale*targets)^2))
      path_reference <- sqrt(sum((paths/edge_scale-ref[p$pair_matrix])^2)/sum(ref[p$pair_matrix]^2))
      recorded <- result$scores[result$scores$method==method,]
      zc <- scale(z,scale=FALSE)
      rot <- svd(crossprod(zc,xc))
      aligned <- sum(rot$d)/sum(zc^2)*zc %*% (rot$u %*% t(rot$v))
      coord <- sqrt(sum((aligned-xc)^2)/sum(xc^2))
      xy <- sqrt(sum((aligned[,1:2]-xc[,1:2])^2)/sum(xc[,1:2]^2))
      height <- sqrt(sum((aligned[,3]-xc[,3])^2)/sum(xc[,3]^2))
      height_only <- sqrt(sum(xc[,1:2]^2)/sum(xc^2))
      pos <- pos+1L
      ans[[pos]] <- data.frame(case=key,regime=regime,method=method,path_rel=path_rel,path_reference=path_reference,
        path_rel_discrepancy=abs(path_rel-recorded$path_rel),path_reference_discrepancy=abs(path_reference-recorded$path_reference),
        coordinate_error=coord,coordinate_discrepancy=abs(coord-recorded$procrustes),horizontal_error=xy,
        vertical_error=height,height_only_error_upper_bound=height_only,min_normalized_edge=min(p$edge_targets))
    }
  }
}
out <- do.call(rbind,ans)
write.csv(out,file.path(dest,"independent-check.csv"),row.names=FALSE)
print(out[,c("case","regime","method","path_rel_discrepancy","path_reference_discrepancy")])
stopifnot(max(out$path_rel_discrepancy,out$path_reference_discrepancy,out$coordinate_discrepancy)<1e-9)
sha <- function(path) substr(system2("shasum",c("-a","256",shQuote(path)),stdout=TRUE),1,64)
manifest <- list(scope="All recorded routes in 24 selected configurations, no refitting; caches required to rerun this optional check.",
  coordinate_source_commit="a43edfef1bbee70e825e790677abb35f12dfff53",
  sources=lapply(sort(unique(sources)),function(path) list(path=substring(path,nchar(root)+2L),sha256=sha(path))),
  output_sha256=sha(file.path(dest,"independent-check.csv")),configurations=nrow(out),routes_per_configuration=28680L)
jsonlite::write_json(manifest,file.path(dest,"source-manifest.json"),pretty=TRUE,auto_unbox=TRUE)
cat("PASS:",nrow(out),"configurations, all 28,680 recorded routes per configuration.\n")
