# Run from the grip repository root; evaluates saved fits, never refits them.
pkgload::load_all(".", quiet=TRUE)
base <- "papers/grip-software-paper"
source(file.path(base,"reproducibility/experiments/two-fidelity-pilot/score-saddle-reference.R"))
args <- commandArgs(TRUE)
input <- if(length(args)) args[1] else file.path(base,"build/two-fidelity-pilot")
files <- list.files(input, "^fit-[0-9]+-k[0-9]+[.]rds$", full.names=TRUE)
stopifnot(length(files)==5L)
rows <- list(); clouds <- list()
for (file in files) {
  f <- readRDS(file); message("Reference diagnostics: cloud ",f$replicate)
  X <- f$coords; tri <- geometry::delaunayn(X[,1:2])
  ref <- saddle.reference.mesh(X,tri,.8,2L)
  aligned <- list()
  for (method in names(f$candidates)) for (alignment in c("rigid","similarity")) {
    fit <- grip::score.coordinates(f$candidates[[method]],X,alignment=alignment)
    if(alignment=="similarity") {
      old <- f$scores$procrustes[match(method,f$scores$method)]
      stopifnot(abs(fit$relative_rmse-old)<1e-10)
      aligned[[method]] <- fit$coords
    }
    for(ns in c(2000L,8000L)) {
      s <- grip::score.surface(fit$coords,tri,ref$coords,ref$triangles,
                             sample_size=ns,seed=1901L+f$replicate)
      rows[[length(rows)+1L]] <- data.frame(replicate=f$replicate,k=f$k,
        method=method,alignment=alignment,sample_size=ns,
        coordinate_relative_rmse=fit$relative_rmse,coordinate_rmse=fit$rmse,
        alignment_scale=fit$scale,surface_rms=s$rms,surface_mc_se=s$rms_mc_se)
    }
  }
  clouds[[as.character(f$replicate)]] <- list(coords=X,triangles=tri,
    candidates=f$candidates,aligned=aligned,k=f$k)
}
scores <- do.call(rbind,rows)
stopifnot(nrow(scores)==60L,all(is.finite(scores$surface_rms)))
out <- file.path(base,"reproducibility/precomputed")
bundle <- list(scores=scores,clouds=clouds,provenance=list(
  fit_version="grip 0.2.0",diagnostics_commit="b72f61d",
  ivue_commit="872f9d4",reference_subdivisions=2L,
  surface_seed="1901 + replicate",input_md5=unname(tools::md5sum(files)),
  R=R.version.string,geometry=as.character(packageVersion("geometry"))))
saveRDS(bundle,file.path(out,"saddle-reference-diagnostics.rds"),compress="xz")
write.csv(scores,file.path(out,"saddle-reference-diagnostics.csv"),row.names=FALSE)
# Preserve the separately run n=500 sensitivity experiment, not as pilot data.
single_input <- if(length(args)>1L) args[2] else "output/single-saddle-reference-2026-09-04"
single <- read.csv(file.path(single_input,"scores.csv"))
write.csv(single,file.path(out,"single-saddle-reference.csv"),row.names=FALSE)
print(aggregate(cbind(coordinate_relative_rmse,surface_rms)~method+alignment,
                subset(scores,sample_size==8000),median))
