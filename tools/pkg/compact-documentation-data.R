#!/usr/bin/env Rscript
# Derive small installed documentation bundles without rerunning any fits.
args <- commandArgs(TRUE)
stopifnot(length(args) <= 1L, all(args %in% "--check"))
check <- identical(args, "--check")
for (kind in c("metric-mds-comparison", "sparse-mds-comparison")) {
  source <- file.path("tools", "pkg", "fixtures", "documentation-data", kind, "benchmark.rds")
  target <- file.path("inst", "extdata", kind, "benchmark.rds")
  full <- readRDS(source)
  compact <- full
  for (i in seq_along(compact$fits)) {
    row <- compact$fits[[i]]$row
    keep <- row$seed == 1L &&
      (kind == "sparse-mds-comparison" || row$budget == 100L)
    if (!keep) compact$fits[[i]]$coords <- NULL
  }
  if (kind == "sparse-mds-comparison")
    for (i in seq_along(compact$evaluation)) compact$evaluation[[i]]$targets <- NULL
  compact$storage <- list(
    kind = "documentation",
    summary.records = "All attempted fits, failures, scores and timings retained unchanged",
    coordinates = if (kind == "sparse-mds-comparison") "Seed 1" else "Seed 1, budget 100",
    full.bundle = source
  )
  stopifnot(identical(lapply(full$fits, `[[`, "row"),
                      lapply(compact$fits, `[[`, "row")),
            identical(full$cases, compact$cases),
            identical(full$protocol, compact$protocol))
  for (i in seq_along(compact$fits))
    if (!is.null(compact$fits[[i]]$coords))
      stopifnot(identical(full$fits[[i]]$coords, compact$fits[[i]]$coords))
  if (check) stopifnot(identical(readRDS(target), compact)) else {
    saveRDS(compact, target, compress = "xz", version = 2)
    stopifnot(identical(readRDS(target), compact))
  }
  cat(kind, ":", length(full$fits), "unchanged summary records;",
      sum(vapply(compact$fits, function(f) !is.null(f$coords), TRUE)),
      "unchanged display layouts;", file.info(target)$size, "bytes\n")
}
