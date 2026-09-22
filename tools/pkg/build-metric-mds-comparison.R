#!/usr/bin/env Rscript
# Explicit developer benchmark; never executed by a vignette or package install.
# Usage: Rscript tools/pkg/build-metric-mds-comparison.R <fresh-output-dir> <candidate-library> <install-log>
args <- commandArgs(TRUE)
if (length(args) != 3L) stop('Supply a fresh output directory, installed candidate library, and clean-install compilation log.')
build.log <- readLines(args[[3]], warn = FALSE)
compile.line <- grep('-c metric_mds_sgd.cpp', build.log, value = TRUE, fixed = TRUE)
if (length(compile.line) != 1L) stop('The log must show one fresh compilation of metric_mds_sgd.cpp; install with --preclean.')
flags <- regmatches(compile.line, gregexpr('-O[0-3sg]', compile.line))[[1]]
if (!length(flags) || tail(flags, 1L) == '-O0') stop('Benchmark requires an optimized build, not reused development objects or -O0.')
if (!requireNamespace('callr', quietly = TRUE)) stop('Install callr before running this benchmark.')
if (!requireNamespace('smacof', quietly = TRUE)) stop('Install smacof before running this benchmark.')
out <- args[[1]]
if (dir.exists(out)) stop('Output already exists; preserve prior attempts and choose a fresh directory.')
dir.create(out, recursive = TRUE)
out <- normalizePath(out)
lib <- normalizePath(args[[2]], mustWork = TRUE)
.libPaths(c(lib, .libPaths()))
library(grip, lib.loc = lib)
stopifnot(identical(eval(formals(grip::metric.mds)$backend), c('sgd','smacof')))
recipe <- file.path(out, 'executed-recipe.R')
file.copy('inst/scripts/metric-mds-comparison.R', recipe)
source(recipe)
started <- proc.time()[['elapsed']]
protocol <- list(sides = c(8L,12L), seeds = 1:5, budgets = c(10L,30L,100L),
  per_fit_timeout_seconds = 30, campaign_seconds = 1800,
  representative_seed = 1L, representative_budget = 100L,
  order = 'Alternate backend order for each dataset/seed/budget block',
  tuning = 'Provisional calibration choices; no outcome-driven tuning',
  timing = 'Wall-clock fit includes wrapper, normalization and checkpoint scoring; excludes shared preparation, process startup, namespace loading and independent final scoring',
  selection = 'Show each independent schedule; no best-of-budgets selection or equal-iteration cost claim',
  source_commit = system('git rev-parse HEAD', intern = TRUE),
  source_md5 = as.list(tools::md5sum(c('R/metric_mds.R','R/metric_mds_sgd.R',
    'src/metric_mds_sgd.cpp','inst/scripts/metric-mds-comparison.R',
    'tools/pkg/build-metric-mds-comparison.R'))),
  fitting_recipe_source = lapply(c('comparison_cases','comparison_fit','comparison_start','comparison_score'),
    function(name) list(name = name, source = deparse(get(name)))),
  RNG = RNGkind(),
  compilation = compile.line,
  installed_native_md5 = unname(tools::md5sum(list.files(system.file('libs', package = 'grip'),
    pattern = '[.](so|dll|dylib)$', full.names = TRUE, recursive = TRUE))),
  R = R.version.string, platform = R.version$platform,
  system = as.list(Sys.info()[c('sysname','release','machine')]),
  package_versions = vapply(c('grip','smacof','ivue','rgl'), function(p)
    if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else NA_character_, ''),
  recorded_utc = format(Sys.time(), tz = 'UTC', usetz = TRUE))
saveRDS(protocol, file.path(out, 'protocol.rds'))
cases <- comparison_cases(protocol$sides)
bundle <- list(protocol = protocol, cases = cases, baseline = list(), fits = list())
block <- 0L
for (id in names(cases)) {
  case <- cases[[id]]
  cat('Preparing', id, '\n'); flush.console()
  begin <- proc.time()[['elapsed']]
  prepared <- grip::prepare.graph.geodesic.mds(edges = case$edges, n = case$n, edge_weights = case$weights)
  bundle$cases[[id]]$preparation_seconds <- proc.time()[['elapsed']] - begin
  D <- prepared$distance_matrix
  bundle$cases[[id]]$distances <- D
  if (case$target == 'Ambient Euclidean') stopifnot(isTRUE(all.equal(D, as.matrix(dist(case$X)), check.attributes = FALSE, tolerance = 1e-10)))
  begin <- proc.time()[['elapsed']]
  classical <- grip::classical.mds(prepared = prepared, dim = case$dimension, diagnostics = FALSE)$coords
  v <- as.double(dist(classical)); target <- as.double(as.dist(D))
  classical <- classical * sum(v * target) / sum(v^2)
  elapsed <- proc.time()[['elapsed']] - begin
  bundle$baseline[[id]] <- list(coords = classical, score = comparison_score(classical, D), seconds = elapsed)
  for (seed in protocol$seeds) for (budget in protocol$budgets) {
    block <- block + 1L
    backends <- if (block %% 2L) c('sgd','smacof') else c('smacof','sgd')
    for (backend in backends) {
      worker.begin <- proc.time()[['elapsed']]
      remaining <- protocol$campaign_seconds - (worker.begin - started)
      result <- if (remaining <= 0) simpleError('Campaign allowance exhausted') else tryCatch(
        callr::r(function(case, prepared, seed, budget, backend, recipe, lib) {
          .libPaths(c(lib, .libPaths()))
          library(grip, lib.loc = lib)
          loadNamespace('smacof')
          source(recipe)
          comparison_fit(case, prepared, seed, budget, backend)
        }, args = list(case, prepared, seed, budget, backend, recipe, lib),
        libpath = .libPaths(), timeout = min(remaining, protocol$per_fit_timeout_seconds)),
        error = function(e) e)
      if (inherits(result, 'error')) {
        result <- list(coords = NULL, row = data.frame(case = id, seed = seed, budget = budget,
          backend = backend, status = if (remaining <= 0) 'not_run_campaign_limit' else 'worker_error_or_timeout',
          seconds = NA_real_, raw_stress = NA_real_, error = NA_real_, iterations = NA_integer_,
          termination = NA_character_, converged = FALSE, warnings = '', failure = conditionMessage(result)))
      }
      result$row$worker_seconds <- proc.time()[['elapsed']] - worker.begin
      key <- paste(id, seed, budget, backend, sep = '_')
      bundle$fits[[key]] <- result
      saveRDS(bundle, file.path(out, 'progress.rds'), compress = FALSE)
      if (length(bundle$fits) %% 20L == 0L) { cat('Recorded', length(bundle$fits), '/ 420 fits\n'); flush.console() }
    }
  }
}
bundle$protocol$total_wall_seconds <- proc.time()[['elapsed']] - started
saveRDS(bundle, file.path(out, 'benchmark.rds'), compress = 'xz')
utils::write.csv(comparison_rows(bundle), file.path(out, 'results.csv'), row.names = FALSE)
utils::write.csv(comparison_summary(bundle), file.path(out, 'summary.csv'), row.names = FALSE)
cat('Completed', length(bundle$fits), 'fits. Review benchmark.rds before publishing.\n')
