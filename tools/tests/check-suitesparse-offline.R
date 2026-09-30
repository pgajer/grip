# Rscript --vanilla tools/tests/check-suitesparse-offline.R /path/to/library
# Validate the installed package, not pkgload's source environment.
args <- commandArgs(TRUE)
if (length(args)) .libPaths(c(args[1], .Library))
library(grip)
for (fn in c('url', 'socketConnection', 'system', 'system2')) {
  trace(fn, tracer=quote(stop('External access forbidden in offline validation')),
        print=FALSE, where=baseenv())
}
trace('download.file', tracer=quote(stop('Network access forbidden')),
      print=FALSE, where=asNamespace('utils'))
catalogue <- suitesparse.graphs()
stopifnot(nrow(catalogue)==85L, sum(catalogue$gallery)==80L)
for (id in catalogue$id) {
  stopifnot(identical(suitesparse.graph(id), suitesparse.graph(sub('^[^/]+/', '', id))))
}
stopifnot(identical(suitesparse.graph('lesmis')$vertex_data$label[1], 'Myriel'))
stopifnot(identical(suitesparse.graph('dwt_307')$n, 307L))
cat('All 85 installed graphs load with network and subprocess calls forbidden.\n')
