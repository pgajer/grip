# Fresh weighted fits; the historical vignette_results.rds is never overwritten.
# Run from the source checkout with an optimized native build already available.
source('inst/scripts/portable-session-info.R')
pkgload::load_all('.',recompile=FALSE,quiet=TRUE)
data(hmp.u01.gc.coarse)
g <- hmp.u01.gc.coarse
fit_one <- function(expr) {
  warnings <- character()
  elapsed <- system.time(value <- withCallingHandlers(expr,warning=function(w) {
    warnings <<- c(warnings,conditionMessage(w));invokeRestart('muffleWarning')
  }))
  list(value=value,elapsed_seconds=unname(elapsed['elapsed']),warnings=warnings)
}
mds <- fit_one(metric.mds(adj.list=g$adj_list,weight.list=g$weight_list,
  n=length(g$adj_list),dim=3,backend='sgd',init='random',seed=1,
  max.iter=200,pair.weights='uniform',diagnostics=FALSE))
weighted <- fit_one(grip(adj.list=g$adj_list,weight.list=g$weight_list,
  n=length(g$adj_list),dim=3,metric='edge_length',seed=1))
mds$value$prepared <- NULL # Do not bundle the rebuildable dense distance matrix.
record <- list(sgd=mds,weighted_grip=weighted,
  settings=list(seed=1,dimension=3,backend='sgd',approximation='full',
    pair_weights='uniform',initialization='random',max_iter=200,
    weighted_grip_metric='edge_length'),
  data_md5=tools::md5sum('data/hmp.u01.gc.coarse.rda'),
  script_md5=tools::md5sum('inst/scripts/precompute-hmp-layout-selector.R'),
  source_md5=tools::md5sum(c(list.files('R',full.names=TRUE,pattern='[.]R$'),
    list.files('src',full.names=TRUE,pattern='[.](cpp|h)$'))),
  session=portable_session_info())
stopifnot(identical(dim(mds$value$coords),c(1828L,3L)),
  identical(dim(weighted$value),c(1828L,3L)),
  all(is.finite(mds$value$coords)),all(is.finite(weighted$value)))
saveRDS(record,'inst/extdata/hmp_u01_gc_coarse/weighted_layouts.rds')
cat('Saved fresh weighted fits.\n')
