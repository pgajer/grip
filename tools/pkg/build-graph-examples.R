#!/usr/bin/env Rscript
# Run from the package root: Rscript tools/pkg/build-graph-examples.R <fresh-output-dir> <candidate-library>
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L, requireNamespace('callr', quietly = TRUE))
if (dir.exists(args[1])) stop('Choose a fresh results directory to retain previous runs.')
dir.create(args[1], recursive = TRUE)
out <- normalizePath(args[1]); lib <- normalizePath(args[2], mustWork = TRUE)
.libPaths(c(lib, .libPaths()))
library(grip, lib.loc = lib)
source('inst/scripts/graph-examples.R')
# The installed candidate must contain zheng.graphs; install before running this script.
cases <- graph_example_cases()
recipe <- file.path(out, 'executed-recipe.R')
file.copy('inst/scripts/graph-examples.R', recipe)
bundle <- list(cases = cases, fits = list(), protocol = list(
  seed = 2026L, dimension = 3L, n_init = 3L, mds_max_iter = 80L,
  edge_kk_max_iter = 100L, stiffness_method = 'uniform', scale_mode = 'identity',
  per_graph_timeout_seconds = 180, settings = 'Illustrative choices, not performance tuning',
  R = R.version.string, platform = R.version$platform, grip_version = as.character(packageVersion('grip')),
  recorded_utc = format(Sys.time(), tz = 'UTC', usetz = TRUE),
  source_commit = system('git rev-parse HEAD', intern = TRUE),
  fitting_recipe_source = lapply(c('graph_example_cases','graph_example_fit','graph_example_scores'),
    function(name) list(name = name, source = deparse(get(name)))),
  source_md5 = as.list(tools::md5sum(c('inst/scripts/graph-examples.R', 'R/metric_mds.R',
    'R/metric_mds_sgd.R', 'src/metric_mds_sgd.cpp', 'R/gmds_layout_interface.R', 'data/zheng.graphs.rda')))))
saveRDS(bundle, file.path(out, 'progress.rds'))
for (id in names(cases)) {
  cat('Fitting', id, 'with', cases[[id]]$n, 'vertices\n'); flush.console()
  bundle$fits[[id]] <- tryCatch(callr::r(function(graph, recipe, lib) {
    .libPaths(c(lib, .libPaths())); library(grip, lib.loc = lib); source(recipe)
    graph_example_fit(graph)
  }, args = list(cases[[id]], recipe, lib), libpath = .libPaths(), timeout = 180),
  error = function(e) list(error = conditionMessage(e)))
  saveRDS(bundle, file.path(out, 'progress.rds'))
}
saveRDS(bundle, file.path(out, 'layouts.rds'), compress = 'xz', version = 2)
failed <- names(Filter(function(f) !is.null(f$error) || !is.null(f$refined$error), bundle$fits))
cat('Saved all', length(cases), 'outcomes. Unsuccessful:', paste(failed, collapse = ', '), '\n')
