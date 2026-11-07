## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(dev="jpeg", dev.args=list(quality=90), fig.retina=1, dpi=96, collapse = TRUE, comment = '#>', fig.width = 9,
                     fig.height = 4, out.width = '100%')
library(grip)
recipe <- system.file('scripts', 'metric-mds-comparison.R', package = 'grip')
source(recipe, local = TRUE)
source(system.file('scripts','vignette-views.R',package='grip'),local=TRUE)
source(system.file('scripts', 'metric-mds-surface-alignment.R', package = 'grip'), local = TRUE)
bundle <- readRDS(system.file('extdata', 'metric-mds-comparison', 'benchmark.rds', package = 'grip'))
rows <- comparison_rows(bundle)
last <- rows[rows$budget == 100L & rows$status == 'ok', ]
paired <- merge(last[last$backend == 'sgd', c('case','seed','error')],
                last[last$backend == 'smacof', c('case','seed','error')],
                by = c('case','seed'), suffixes = c('.sgd','.smacof'))
have.ivue <- isTRUE(getOption('grip.vignette.interactive', FALSE)) &&
  requireNamespace('ivue', quietly = TRUE) && requireNamespace('rgl', quietly = TRUE) &&
  all(c('plot3D.groups','color.scale.groups','layer3D.callback','layer3D.axes','camera.zup') %in% getNamespaceExports('ivue')) &&
  'limits' %in% names(formals(ivue::plot3D.plain))

## ----run-status---------------------------------------------------------------
knitr::kable(as.data.frame(table(Backend = rows$backend, Status = rows$status)),
             col.names = c('Backend','Outcome','Fits'))
knitr::kable(do.call(rbind, lapply(split(rows, rows$backend), function(z) {
  data.frame(Backend = z$backend[1], Fits = nrow(z),
    `Successful fits` = sum(z$status == 'ok'),
    `Total fit seconds` = sum(z$seconds, na.rm = TRUE),
    `Total worker seconds` = sum(z$worker_seconds), check.names = FALSE)
})), row.names = FALSE, digits = 3)

## ----individual-view, eval=FALSE----------------------------------------------
# bundle <- readRDS(system.file('extdata', 'metric-mds-comparison',
#                               'benchmark.rds', package = 'grip'))
# case <- bundle$cases[['saddle_graph-144']]
# fits <- comparison_representative(bundle, case$id)
# surface <- comparison.surface.saved(case, fits)
# comparison_view(case, fits, show = 'sgd', surface.alignments = surface$registrations)
# comparison_view(case, fits, show = 'smacof', surface.alignments = surface$registrations)
# comparison_view(case, fits, show = 'overlay', surface.alignments = surface$registrations)

## ----tradeoff-gaussian, echo=FALSE, fig.cap="Separate schedules: median normalized error versus median fitting time; bars show the interquartile ranges over five starts. Labels are pass/iteration allowances."----
comparison_tradeoff(bundle, c("gaussian-64", "gaussian-144"))

## ----table-gaussian-64, echo=FALSE--------------------------------------------
case <- bundle$cases[["gaussian-64"]]
fits <- comparison_representative(bundle, case$id)
surface <- if (grepl("^saddle_graph-", case$id) ) comparison.surface.saved(case, fits) else NULL
reference <- comparison_score(case$X[,seq_len(case$dimension),drop=FALSE], case$distances)
classical <- bundle$baseline[[case$id]]
tab <- data.frame(Method = c("Generating reference", "Classical MDS", "SGD", "SMACOF"),
  Raw_stress = c(reference[["raw_stress"]], classical$score[["raw_stress"]], fits$sgd$row$raw_stress, fits$smacof$row$raw_stress),
  Normalized_error = c(reference[["error"]], classical$score[["error"]], fits$sgd$row$error, fits$smacof$row$error),
  Fit_seconds = c(NA, classical$seconds, fits$sgd$row$seconds, fits$smacof$row$seconds),
  Iterations = c(NA, NA, fits$sgd$row$iterations, fits$smacof$row$iterations),
  Termination = c("Reference only", "Spectral reference", if (fits$sgd$row$status == "ok") fits$sgd$row$termination else fits$sgd$row$status, if (fits$smacof$row$status == "ok") fits$smacof$row$termination else fits$smacof$row$status))
tab$Raw_stress <- format(tab$Raw_stress, scientific = TRUE, digits = 4)
tab$Normalized_error <- format(tab$Normalized_error, scientific = TRUE, digits = 4)
knitr::kable(tab, digits = 3, col.names = c("Method", "Raw stress", "Normalized error", "Fit seconds", "Passes/iterations", "Termination"))

## ----overlay-gaussian-64, echo=FALSE, fig.height=3.5--------------------------
if (case$dimension == 2L) {
  comparison_planar(case, fits)
} else if (have.ivue) {
  comparison_view(case, fits, limits=if (!is.null(surface)) surface$limits else NULL)
} else {
  documentation.full.view(case,fits)
}

## ----table-gaussian-144, echo=FALSE-------------------------------------------
case <- bundle$cases[["gaussian-144"]]
fits <- comparison_representative(bundle, case$id)
surface <- if (grepl("^saddle_graph-", case$id) ) comparison.surface.saved(case, fits) else NULL
reference <- comparison_score(case$X[,seq_len(case$dimension),drop=FALSE], case$distances)
classical <- bundle$baseline[[case$id]]
tab <- data.frame(Method = c("Generating reference", "Classical MDS", "SGD", "SMACOF"),
  Raw_stress = c(reference[["raw_stress"]], classical$score[["raw_stress"]], fits$sgd$row$raw_stress, fits$smacof$row$raw_stress),
  Normalized_error = c(reference[["error"]], classical$score[["error"]], fits$sgd$row$error, fits$smacof$row$error),
  Fit_seconds = c(NA, classical$seconds, fits$sgd$row$seconds, fits$smacof$row$seconds),
  Iterations = c(NA, NA, fits$sgd$row$iterations, fits$smacof$row$iterations),
  Termination = c("Reference only", "Spectral reference", if (fits$sgd$row$status == "ok") fits$sgd$row$termination else fits$sgd$row$status, if (fits$smacof$row$status == "ok") fits$smacof$row$termination else fits$smacof$row$status))
tab$Raw_stress <- format(tab$Raw_stress, scientific = TRUE, digits = 4)
tab$Normalized_error <- format(tab$Normalized_error, scientific = TRUE, digits = 4)
knitr::kable(tab, digits = 3, col.names = c("Method", "Raw stress", "Normalized error", "Fit seconds", "Passes/iterations", "Termination"))

## ----overlay-gaussian-144, echo=FALSE, fig.height=3.5-------------------------
if (case$dimension == 2L) {
  comparison_planar(case, fits)
} else if (have.ivue) {
  comparison_view(case, fits, limits=if (!is.null(surface)) surface$limits else NULL)
} else {
  documentation.full.view(case,fits)
}

## ----tradeoff-plane, echo=FALSE, fig.cap="Separate schedules: median normalized error versus median fitting time; bars show the interquartile ranges over five starts. Labels are pass/iteration allowances."----
comparison_tradeoff(bundle, c("plane-64", "plane-144"))

## ----table-plane-64, echo=FALSE-----------------------------------------------
case <- bundle$cases[["plane-64"]]
fits <- comparison_representative(bundle, case$id)
surface <- if (grepl("^saddle_graph-", case$id) ) comparison.surface.saved(case, fits) else NULL
reference <- comparison_score(case$X[,seq_len(case$dimension),drop=FALSE], case$distances)
classical <- bundle$baseline[[case$id]]
tab <- data.frame(Method = c("Generating reference", "Classical MDS", "SGD", "SMACOF"),
  Raw_stress = c(reference[["raw_stress"]], classical$score[["raw_stress"]], fits$sgd$row$raw_stress, fits$smacof$row$raw_stress),
  Normalized_error = c(reference[["error"]], classical$score[["error"]], fits$sgd$row$error, fits$smacof$row$error),
  Fit_seconds = c(NA, classical$seconds, fits$sgd$row$seconds, fits$smacof$row$seconds),
  Iterations = c(NA, NA, fits$sgd$row$iterations, fits$smacof$row$iterations),
  Termination = c("Reference only", "Spectral reference", if (fits$sgd$row$status == "ok") fits$sgd$row$termination else fits$sgd$row$status, if (fits$smacof$row$status == "ok") fits$smacof$row$termination else fits$smacof$row$status))
tab$Raw_stress <- format(tab$Raw_stress, scientific = TRUE, digits = 4)
tab$Normalized_error <- format(tab$Normalized_error, scientific = TRUE, digits = 4)
knitr::kable(tab, digits = 3, col.names = c("Method", "Raw stress", "Normalized error", "Fit seconds", "Passes/iterations", "Termination"))

## ----overlay-plane-64, echo=FALSE, fig.height=3.5-----------------------------
if (case$dimension == 2L) {
  comparison_planar(case, fits)
} else if (have.ivue) {
  comparison_view(case, fits, limits=if (!is.null(surface)) surface$limits else NULL)
} else {
  documentation.full.view(case,fits)
}

## ----table-plane-144, echo=FALSE----------------------------------------------
case <- bundle$cases[["plane-144"]]
fits <- comparison_representative(bundle, case$id)
surface <- if (grepl("^saddle_graph-", case$id) ) comparison.surface.saved(case, fits) else NULL
reference <- comparison_score(case$X[,seq_len(case$dimension),drop=FALSE], case$distances)
classical <- bundle$baseline[[case$id]]
tab <- data.frame(Method = c("Generating reference", "Classical MDS", "SGD", "SMACOF"),
  Raw_stress = c(reference[["raw_stress"]], classical$score[["raw_stress"]], fits$sgd$row$raw_stress, fits$smacof$row$raw_stress),
  Normalized_error = c(reference[["error"]], classical$score[["error"]], fits$sgd$row$error, fits$smacof$row$error),
  Fit_seconds = c(NA, classical$seconds, fits$sgd$row$seconds, fits$smacof$row$seconds),
  Iterations = c(NA, NA, fits$sgd$row$iterations, fits$smacof$row$iterations),
  Termination = c("Reference only", "Spectral reference", if (fits$sgd$row$status == "ok") fits$sgd$row$termination else fits$sgd$row$status, if (fits$smacof$row$status == "ok") fits$smacof$row$termination else fits$smacof$row$status))
tab$Raw_stress <- format(tab$Raw_stress, scientific = TRUE, digits = 4)
tab$Normalized_error <- format(tab$Normalized_error, scientific = TRUE, digits = 4)
knitr::kable(tab, digits = 3, col.names = c("Method", "Raw stress", "Normalized error", "Fit seconds", "Passes/iterations", "Termination"))

## ----overlay-plane-144, echo=FALSE, fig.height=3.5----------------------------
if (case$dimension == 2L) {
  comparison_planar(case, fits)
} else if (have.ivue) {
  comparison_view(case, fits, limits=if (!is.null(surface)) surface$limits else NULL)
} else {
  documentation.full.view(case,fits)
}

## ----tradeoff-paraboloid, echo=FALSE, fig.cap="Separate schedules: median normalized error versus median fitting time; bars show the interquartile ranges over five starts. Labels are pass/iteration allowances."----
comparison_tradeoff(bundle, c("paraboloid-64", "paraboloid-144"))

## ----table-paraboloid-64, echo=FALSE------------------------------------------
case <- bundle$cases[["paraboloid-64"]]
fits <- comparison_representative(bundle, case$id)
surface <- if (grepl("^saddle_graph-", case$id) ) comparison.surface.saved(case, fits) else NULL
reference <- comparison_score(case$X[,seq_len(case$dimension),drop=FALSE], case$distances)
classical <- bundle$baseline[[case$id]]
tab <- data.frame(Method = c("Generating reference", "Classical MDS", "SGD", "SMACOF"),
  Raw_stress = c(reference[["raw_stress"]], classical$score[["raw_stress"]], fits$sgd$row$raw_stress, fits$smacof$row$raw_stress),
  Normalized_error = c(reference[["error"]], classical$score[["error"]], fits$sgd$row$error, fits$smacof$row$error),
  Fit_seconds = c(NA, classical$seconds, fits$sgd$row$seconds, fits$smacof$row$seconds),
  Iterations = c(NA, NA, fits$sgd$row$iterations, fits$smacof$row$iterations),
  Termination = c("Reference only", "Spectral reference", if (fits$sgd$row$status == "ok") fits$sgd$row$termination else fits$sgd$row$status, if (fits$smacof$row$status == "ok") fits$smacof$row$termination else fits$smacof$row$status))
tab$Raw_stress <- format(tab$Raw_stress, scientific = TRUE, digits = 4)
tab$Normalized_error <- format(tab$Normalized_error, scientific = TRUE, digits = 4)
knitr::kable(tab, digits = 3, col.names = c("Method", "Raw stress", "Normalized error", "Fit seconds", "Passes/iterations", "Termination"))

## ----overlay-paraboloid-64, echo=FALSE, fig.height=3.5------------------------
if (case$dimension == 2L) {
  comparison_planar(case, fits)
} else if (have.ivue) {
  comparison_view(case, fits, limits=if (!is.null(surface)) surface$limits else NULL)
} else {
  documentation.full.view(case,fits)
}

## ----table-paraboloid-144, echo=FALSE-----------------------------------------
case <- bundle$cases[["paraboloid-144"]]
fits <- comparison_representative(bundle, case$id)
surface <- if (grepl("^saddle_graph-", case$id) ) comparison.surface.saved(case, fits) else NULL
reference <- comparison_score(case$X[,seq_len(case$dimension),drop=FALSE], case$distances)
classical <- bundle$baseline[[case$id]]
tab <- data.frame(Method = c("Generating reference", "Classical MDS", "SGD", "SMACOF"),
  Raw_stress = c(reference[["raw_stress"]], classical$score[["raw_stress"]], fits$sgd$row$raw_stress, fits$smacof$row$raw_stress),
  Normalized_error = c(reference[["error"]], classical$score[["error"]], fits$sgd$row$error, fits$smacof$row$error),
  Fit_seconds = c(NA, classical$seconds, fits$sgd$row$seconds, fits$smacof$row$seconds),
  Iterations = c(NA, NA, fits$sgd$row$iterations, fits$smacof$row$iterations),
  Termination = c("Reference only", "Spectral reference", if (fits$sgd$row$status == "ok") fits$sgd$row$termination else fits$sgd$row$status, if (fits$smacof$row$status == "ok") fits$smacof$row$termination else fits$smacof$row$status))
tab$Raw_stress <- format(tab$Raw_stress, scientific = TRUE, digits = 4)
tab$Normalized_error <- format(tab$Normalized_error, scientific = TRUE, digits = 4)
knitr::kable(tab, digits = 3, col.names = c("Method", "Raw stress", "Normalized error", "Fit seconds", "Passes/iterations", "Termination"))

## ----overlay-paraboloid-144, echo=FALSE, fig.height=3.5-----------------------
if (case$dimension == 2L) {
  comparison_planar(case, fits)
} else if (have.ivue) {
  comparison_view(case, fits, limits=if (!is.null(surface)) surface$limits else NULL)
} else {
  documentation.full.view(case,fits)
}

## ----tradeoff-saddle, echo=FALSE, fig.cap="Separate schedules: median normalized error versus median fitting time; bars show the interquartile ranges over five starts. Labels are pass/iteration allowances."----
comparison_tradeoff(bundle, c("saddle-64", "saddle-144"))

## ----table-saddle-64, echo=FALSE----------------------------------------------
case <- bundle$cases[["saddle-64"]]
fits <- comparison_representative(bundle, case$id)
surface <- if (grepl("^saddle_graph-", case$id) ) comparison.surface.saved(case, fits) else NULL
reference <- comparison_score(case$X[,seq_len(case$dimension),drop=FALSE], case$distances)
classical <- bundle$baseline[[case$id]]
tab <- data.frame(Method = c("Generating reference", "Classical MDS", "SGD", "SMACOF"),
  Raw_stress = c(reference[["raw_stress"]], classical$score[["raw_stress"]], fits$sgd$row$raw_stress, fits$smacof$row$raw_stress),
  Normalized_error = c(reference[["error"]], classical$score[["error"]], fits$sgd$row$error, fits$smacof$row$error),
  Fit_seconds = c(NA, classical$seconds, fits$sgd$row$seconds, fits$smacof$row$seconds),
  Iterations = c(NA, NA, fits$sgd$row$iterations, fits$smacof$row$iterations),
  Termination = c("Reference only", "Spectral reference", if (fits$sgd$row$status == "ok") fits$sgd$row$termination else fits$sgd$row$status, if (fits$smacof$row$status == "ok") fits$smacof$row$termination else fits$smacof$row$status))
tab$Raw_stress <- format(tab$Raw_stress, scientific = TRUE, digits = 4)
tab$Normalized_error <- format(tab$Normalized_error, scientific = TRUE, digits = 4)
knitr::kable(tab, digits = 3, col.names = c("Method", "Raw stress", "Normalized error", "Fit seconds", "Passes/iterations", "Termination"))

## ----overlay-saddle-64, echo=FALSE, fig.height=3.5----------------------------
if (case$dimension == 2L) {
  comparison_planar(case, fits)
} else if (have.ivue) {
  comparison_view(case, fits, limits=if (!is.null(surface)) surface$limits else NULL)
} else {
  documentation.full.view(case,fits)
}

## ----table-saddle-144, echo=FALSE---------------------------------------------
case <- bundle$cases[["saddle-144"]]
fits <- comparison_representative(bundle, case$id)
surface <- if (grepl("^saddle_graph-", case$id) ) comparison.surface.saved(case, fits) else NULL
reference <- comparison_score(case$X[,seq_len(case$dimension),drop=FALSE], case$distances)
classical <- bundle$baseline[[case$id]]
tab <- data.frame(Method = c("Generating reference", "Classical MDS", "SGD", "SMACOF"),
  Raw_stress = c(reference[["raw_stress"]], classical$score[["raw_stress"]], fits$sgd$row$raw_stress, fits$smacof$row$raw_stress),
  Normalized_error = c(reference[["error"]], classical$score[["error"]], fits$sgd$row$error, fits$smacof$row$error),
  Fit_seconds = c(NA, classical$seconds, fits$sgd$row$seconds, fits$smacof$row$seconds),
  Iterations = c(NA, NA, fits$sgd$row$iterations, fits$smacof$row$iterations),
  Termination = c("Reference only", "Spectral reference", if (fits$sgd$row$status == "ok") fits$sgd$row$termination else fits$sgd$row$status, if (fits$smacof$row$status == "ok") fits$smacof$row$termination else fits$smacof$row$status))
tab$Raw_stress <- format(tab$Raw_stress, scientific = TRUE, digits = 4)
tab$Normalized_error <- format(tab$Normalized_error, scientific = TRUE, digits = 4)
knitr::kable(tab, digits = 3, col.names = c("Method", "Raw stress", "Normalized error", "Fit seconds", "Passes/iterations", "Termination"))

## ----overlay-saddle-144, echo=FALSE, fig.height=3.5---------------------------
if (case$dimension == 2L) {
  comparison_planar(case, fits)
} else if (have.ivue) {
  comparison_view(case, fits, limits=if (!is.null(surface)) surface$limits else NULL)
} else {
  documentation.full.view(case,fits)
}

## ----tradeoff-paraboloid_graph, echo=FALSE, fig.cap="Separate schedules: median normalized error versus median fitting time; bars show the interquartile ranges over five starts. Labels are pass/iteration allowances."----
comparison_tradeoff(bundle, c("paraboloid_graph-64", "paraboloid_graph-144"))

## ----table-paraboloid_graph-64, echo=FALSE------------------------------------
case <- bundle$cases[["paraboloid_graph-64"]]
fits <- comparison_representative(bundle, case$id)
surface <- if (grepl("^saddle_graph-", case$id) ) comparison.surface.saved(case, fits) else NULL
reference <- comparison_score(case$X[,seq_len(case$dimension),drop=FALSE], case$distances)
classical <- bundle$baseline[[case$id]]
tab <- data.frame(Method = c("Generating reference", "Classical MDS", "SGD", "SMACOF"),
  Raw_stress = c(reference[["raw_stress"]], classical$score[["raw_stress"]], fits$sgd$row$raw_stress, fits$smacof$row$raw_stress),
  Normalized_error = c(reference[["error"]], classical$score[["error"]], fits$sgd$row$error, fits$smacof$row$error),
  Fit_seconds = c(NA, classical$seconds, fits$sgd$row$seconds, fits$smacof$row$seconds),
  Iterations = c(NA, NA, fits$sgd$row$iterations, fits$smacof$row$iterations),
  Termination = c("Reference only", "Spectral reference", if (fits$sgd$row$status == "ok") fits$sgd$row$termination else fits$sgd$row$status, if (fits$smacof$row$status == "ok") fits$smacof$row$termination else fits$smacof$row$status))
tab$Raw_stress <- format(tab$Raw_stress, scientific = TRUE, digits = 4)
tab$Normalized_error <- format(tab$Normalized_error, scientific = TRUE, digits = 4)
knitr::kable(tab, digits = 3, col.names = c("Method", "Raw stress", "Normalized error", "Fit seconds", "Passes/iterations", "Termination"))

## ----overlay-paraboloid_graph-64, echo=FALSE, fig.height=3.5------------------
if (case$dimension == 2L) {
  comparison_planar(case, fits)
} else if (have.ivue) {
  comparison_view(case, fits, limits=if (!is.null(surface)) surface$limits else NULL)
} else {
  documentation.full.view(case,fits)
}

## ----table-paraboloid_graph-144, echo=FALSE-----------------------------------
case <- bundle$cases[["paraboloid_graph-144"]]
fits <- comparison_representative(bundle, case$id)
surface <- if (grepl("^saddle_graph-", case$id) ) comparison.surface.saved(case, fits) else NULL
reference <- comparison_score(case$X[,seq_len(case$dimension),drop=FALSE], case$distances)
classical <- bundle$baseline[[case$id]]
tab <- data.frame(Method = c("Generating reference", "Classical MDS", "SGD", "SMACOF"),
  Raw_stress = c(reference[["raw_stress"]], classical$score[["raw_stress"]], fits$sgd$row$raw_stress, fits$smacof$row$raw_stress),
  Normalized_error = c(reference[["error"]], classical$score[["error"]], fits$sgd$row$error, fits$smacof$row$error),
  Fit_seconds = c(NA, classical$seconds, fits$sgd$row$seconds, fits$smacof$row$seconds),
  Iterations = c(NA, NA, fits$sgd$row$iterations, fits$smacof$row$iterations),
  Termination = c("Reference only", "Spectral reference", if (fits$sgd$row$status == "ok") fits$sgd$row$termination else fits$sgd$row$status, if (fits$smacof$row$status == "ok") fits$smacof$row$termination else fits$smacof$row$status))
tab$Raw_stress <- format(tab$Raw_stress, scientific = TRUE, digits = 4)
tab$Normalized_error <- format(tab$Normalized_error, scientific = TRUE, digits = 4)
knitr::kable(tab, digits = 3, col.names = c("Method", "Raw stress", "Normalized error", "Fit seconds", "Passes/iterations", "Termination"))

## ----overlay-paraboloid_graph-144, echo=FALSE, fig.height=3.5-----------------
if (case$dimension == 2L) {
  comparison_planar(case, fits)
} else if (have.ivue) {
  comparison_view(case, fits, limits=if (!is.null(surface)) surface$limits else NULL)
} else {
  documentation.full.view(case,fits)
}

## ----tradeoff-saddle_graph, echo=FALSE, fig.cap="Separate schedules: median normalized error versus median fitting time; bars show the interquartile ranges over five starts. Labels are pass/iteration allowances."----
comparison_tradeoff(bundle, c("saddle_graph-64", "saddle_graph-144"))

## ----table-saddle_graph-64, echo=FALSE----------------------------------------
case <- bundle$cases[["saddle_graph-64"]]
fits <- comparison_representative(bundle, case$id)
surface <- if (grepl("^saddle_graph-", case$id) ) comparison.surface.saved(case, fits) else NULL
reference <- comparison_score(case$X[,seq_len(case$dimension),drop=FALSE], case$distances)
classical <- bundle$baseline[[case$id]]
tab <- data.frame(Method = c("Generating reference", "Classical MDS", "SGD", "SMACOF"),
  Raw_stress = c(reference[["raw_stress"]], classical$score[["raw_stress"]], fits$sgd$row$raw_stress, fits$smacof$row$raw_stress),
  Normalized_error = c(reference[["error"]], classical$score[["error"]], fits$sgd$row$error, fits$smacof$row$error),
  Fit_seconds = c(NA, classical$seconds, fits$sgd$row$seconds, fits$smacof$row$seconds),
  Iterations = c(NA, NA, fits$sgd$row$iterations, fits$smacof$row$iterations),
  Termination = c("Reference only", "Spectral reference", if (fits$sgd$row$status == "ok") fits$sgd$row$termination else fits$sgd$row$status, if (fits$smacof$row$status == "ok") fits$smacof$row$termination else fits$smacof$row$status))
tab$Raw_stress <- format(tab$Raw_stress, scientific = TRUE, digits = 4)
tab$Normalized_error <- format(tab$Normalized_error, scientific = TRUE, digits = 4)
knitr::kable(tab, digits = 3, col.names = c("Method", "Raw stress", "Normalized error", "Fit seconds", "Passes/iterations", "Termination"))

## ----surface-shape-saddle_graph-64, echo=FALSE--------------------------------
if (!is.null(surface)) knitr::kable(surface$audit, digits=4, row.names=FALSE, col.names=c("Backend", "Subdivisions per edge", "Surface RMS before", "Surface RMS after")) else cat("Install FNN to calculate surface alignment.")

## ----surface-overlay-saddle_graph-64, echo=FALSE------------------------------
documentation.full.view(case,fits)

## ----table-saddle_graph-144, echo=FALSE---------------------------------------
case <- bundle$cases[["saddle_graph-144"]]
fits <- comparison_representative(bundle, case$id)
surface <- if (grepl("^saddle_graph-", case$id) ) comparison.surface.saved(case, fits) else NULL
reference <- comparison_score(case$X[,seq_len(case$dimension),drop=FALSE], case$distances)
classical <- bundle$baseline[[case$id]]
tab <- data.frame(Method = c("Generating reference", "Classical MDS", "SGD", "SMACOF"),
  Raw_stress = c(reference[["raw_stress"]], classical$score[["raw_stress"]], fits$sgd$row$raw_stress, fits$smacof$row$raw_stress),
  Normalized_error = c(reference[["error"]], classical$score[["error"]], fits$sgd$row$error, fits$smacof$row$error),
  Fit_seconds = c(NA, classical$seconds, fits$sgd$row$seconds, fits$smacof$row$seconds),
  Iterations = c(NA, NA, fits$sgd$row$iterations, fits$smacof$row$iterations),
  Termination = c("Reference only", "Spectral reference", if (fits$sgd$row$status == "ok") fits$sgd$row$termination else fits$sgd$row$status, if (fits$smacof$row$status == "ok") fits$smacof$row$termination else fits$smacof$row$status))
tab$Raw_stress <- format(tab$Raw_stress, scientific = TRUE, digits = 4)
tab$Normalized_error <- format(tab$Normalized_error, scientific = TRUE, digits = 4)
knitr::kable(tab, digits = 3, col.names = c("Method", "Raw stress", "Normalized error", "Fit seconds", "Passes/iterations", "Termination"))

## ----surface-shape-saddle_graph-144, echo=FALSE-------------------------------
if (!is.null(surface)) knitr::kable(surface$audit, digits=4, row.names=FALSE, col.names=c("Backend", "Subdivisions per edge", "Surface RMS before", "Surface RMS after")) else cat("Install FNN to calculate surface alignment.")

## ----surface-overlay-saddle_graph-144, echo=FALSE-----------------------------
documentation.full.view(case,fits)

## ----tradeoff-helix_graph, echo=FALSE, fig.cap="Separate schedules: median normalized error versus median fitting time; bars show the interquartile ranges over five starts. Labels are pass/iteration allowances."----
comparison_tradeoff(bundle, c("helix_graph-64", "helix_graph-144"))

## ----table-helix_graph-64, echo=FALSE-----------------------------------------
case <- bundle$cases[["helix_graph-64"]]
fits <- comparison_representative(bundle, case$id)
surface <- if (grepl("^saddle_graph-", case$id) ) comparison.surface.saved(case, fits) else NULL
reference <- comparison_score(case$X[,seq_len(case$dimension),drop=FALSE], case$distances)
classical <- bundle$baseline[[case$id]]
tab <- data.frame(Method = c("Generating reference", "Classical MDS", "SGD", "SMACOF"),
  Raw_stress = c(reference[["raw_stress"]], classical$score[["raw_stress"]], fits$sgd$row$raw_stress, fits$smacof$row$raw_stress),
  Normalized_error = c(reference[["error"]], classical$score[["error"]], fits$sgd$row$error, fits$smacof$row$error),
  Fit_seconds = c(NA, classical$seconds, fits$sgd$row$seconds, fits$smacof$row$seconds),
  Iterations = c(NA, NA, fits$sgd$row$iterations, fits$smacof$row$iterations),
  Termination = c("Reference only", "Spectral reference", if (fits$sgd$row$status == "ok") fits$sgd$row$termination else fits$sgd$row$status, if (fits$smacof$row$status == "ok") fits$smacof$row$termination else fits$smacof$row$status))
tab$Raw_stress <- format(tab$Raw_stress, scientific = TRUE, digits = 4)
tab$Normalized_error <- format(tab$Normalized_error, scientific = TRUE, digits = 4)
knitr::kable(tab, digits = 3, col.names = c("Method", "Raw stress", "Normalized error", "Fit seconds", "Passes/iterations", "Termination"))

## ----overlay-helix_graph-64, echo=FALSE, fig.height=3.5-----------------------
if (case$dimension == 2L) {
  comparison_planar(case, fits)
} else if (have.ivue) {
  comparison_view(case, fits, limits=if (!is.null(surface)) surface$limits else NULL)
} else {
  documentation.full.view(case,fits)
}

## ----table-helix_graph-144, echo=FALSE----------------------------------------
case <- bundle$cases[["helix_graph-144"]]
fits <- comparison_representative(bundle, case$id)
surface <- if (grepl("^saddle_graph-", case$id) ) comparison.surface.saved(case, fits) else NULL
reference <- comparison_score(case$X[,seq_len(case$dimension),drop=FALSE], case$distances)
classical <- bundle$baseline[[case$id]]
tab <- data.frame(Method = c("Generating reference", "Classical MDS", "SGD", "SMACOF"),
  Raw_stress = c(reference[["raw_stress"]], classical$score[["raw_stress"]], fits$sgd$row$raw_stress, fits$smacof$row$raw_stress),
  Normalized_error = c(reference[["error"]], classical$score[["error"]], fits$sgd$row$error, fits$smacof$row$error),
  Fit_seconds = c(NA, classical$seconds, fits$sgd$row$seconds, fits$smacof$row$seconds),
  Iterations = c(NA, NA, fits$sgd$row$iterations, fits$smacof$row$iterations),
  Termination = c("Reference only", "Spectral reference", if (fits$sgd$row$status == "ok") fits$sgd$row$termination else fits$sgd$row$status, if (fits$smacof$row$status == "ok") fits$smacof$row$termination else fits$smacof$row$status))
tab$Raw_stress <- format(tab$Raw_stress, scientific = TRUE, digits = 4)
tab$Normalized_error <- format(tab$Normalized_error, scientific = TRUE, digits = 4)
knitr::kable(tab, digits = 3, col.names = c("Method", "Raw stress", "Normalized error", "Fit seconds", "Passes/iterations", "Termination"))

## ----overlay-helix_graph-144, echo=FALSE, fig.height=3.5----------------------
if (case$dimension == 2L) {
  comparison_planar(case, fits)
} else if (have.ivue) {
  comparison_view(case, fits, limits=if (!is.null(surface)) surface$limits else NULL)
} else {
  documentation.full.view(case,fits)
}

## ----experiment-figures, results='asis', echo=FALSE---------------------------
kinds <- c('gaussian','plane','paraboloid','saddle','paraboloid_graph','saddle_graph','helix_graph')
for (kind in kinds) {
  ids <- paste(kind, c(64,144), sep = '-')
  case <- bundle$cases[[ids[1]]]
  cat('\n\n## ', case$label, ': ', case$target, '\n\n', sep = '')
  child <- c(
    sprintf('```{r tradeoff-%s, echo=FALSE, fig.cap="Separate schedules: median normalized error versus median fitting time; bars show the interquartile ranges over five starts. Labels are pass/iteration allowances."}', kind),
    sprintf('comparison_tradeoff(bundle, c("%s", "%s"))', ids[1], ids[2]), '```')
  cat(knitr::knit_child(text = child, envir = environment(), quiet = TRUE), sep = '\n')
  for (id in ids) {
    child <- c(sprintf('### %s points', bundle$cases[[id]]$n),
      sprintf('```{r table-%s, echo=FALSE}', id),
      sprintf('case <- bundle$cases[["%s"]]', id),
      'fits <- comparison_representative(bundle, case$id)',
      'surface <- if (grepl("^saddle_graph-", case$id) ) comparison.surface.saved(case, fits) else NULL',
      'reference <- comparison_score(case$X[,seq_len(case$dimension),drop=FALSE], case$distances)',
      'classical <- bundle$baseline[[case$id]]',
      'tab <- data.frame(Method = c("Generating reference", "Classical MDS", "SGD", "SMACOF"),',
      '  Raw_stress = c(reference[["raw_stress"]], classical$score[["raw_stress"]], fits$sgd$row$raw_stress, fits$smacof$row$raw_stress),',
      '  Normalized_error = c(reference[["error"]], classical$score[["error"]], fits$sgd$row$error, fits$smacof$row$error),',
      '  Fit_seconds = c(NA, classical$seconds, fits$sgd$row$seconds, fits$smacof$row$seconds),',
      '  Iterations = c(NA, NA, fits$sgd$row$iterations, fits$smacof$row$iterations),',
      '  Termination = c("Reference only", "Spectral reference", if (fits$sgd$row$status == "ok") fits$sgd$row$termination else fits$sgd$row$status, if (fits$smacof$row$status == "ok") fits$smacof$row$termination else fits$smacof$row$status))',
      'tab$Raw_stress <- format(tab$Raw_stress, scientific = TRUE, digits = 4)',
      'tab$Normalized_error <- format(tab$Normalized_error, scientific = TRUE, digits = 4)',
      'knitr::kable(tab, digits = 3, col.names = c("Method", "Raw stress", "Normalized error", "Fit seconds", "Passes/iterations", "Termination"))', '```',
      'Fixed seed 1, allowance 100. Gray: generating reference; blue: SGD; orange: SMACOF. The table above reports the exact displayed fits. Alignment preserves scale.', '')
    if (kind != 'saddle_graph') child <- c(child,
      sprintf('```{r overlay-%s, echo=FALSE, fig.height=3.5}', id),
      'if (case$dimension == 2L) {',
      '  comparison_planar(case, fits)',
      '} else if (have.ivue) {',
      '  comparison_view(case, fits, limits=if (!is.null(surface)) surface$limits else NULL)',
      '} else {',
      '  documentation.full.view(case,fits)',
      '}', '```')
    cat(knitr::knit_child(text = child, envir = environment(), quiet = TRUE), sep = '\n')
    if (kind == 'saddle_graph') {
      child <- c('#### Surface-shape alignment', '',
        'The gray reference stays fixed. The blue and orange layouts are moved to match surface shape rather than vertex labels. Choose SGD, SMACOF, or both.', '',
        sprintf('```{r surface-shape-%s, echo=FALSE}', id),
        'if (!is.null(surface)) knitr::kable(surface$audit, digits=4, row.names=FALSE, col.names=c("Backend", "Subdivisions per edge", "Surface RMS before", "Surface RMS after")) else cat("Install FNN to calculate surface alignment.")', '```', '',
        sprintf('```{r surface-overlay-%s, echo=FALSE}', id),
        'documentation.full.view(case,fits)',
        '```', '')
      cat(knitr::knit_child(text=child, envir=environment(), quiet=TRUE), sep='\n')
    }
  }
}

## ----weighted-api, eval=FALSE-------------------------------------------------
# g <- suitesparse.graph('dwt_307')
# prepared <- prepare.graph.geodesic.mds(edges=g$edges,n=g$n,edge.weights=g$edge_weights)
# start <- withr::with_seed(1,matrix(rnorm(g$n*3),g$n,3))
# sgd <- metric.mds(prepared = prepared, dim = 3, init = start,
#                  backend = "sgd", pair.weights = "inverse_squared",
#                  max.iter = 100, diagnostics = FALSE)
# smacof <- metric.mds(prepared = prepared, dim = 3, init = start,
#                     backend = "smacof", pair.weights = "inverse_squared",
#                     max.iter = 100, eps = 1e-8, diagnostics = FALSE)

## ----weighting-results, include=FALSE-----------------------------------------
weighting.bundle <- readRDS(system.file('extdata','metric-mds-weighting',
                                      'benchmark.rds',package='grip'))
weighting.rows <- comparison_rows(weighting.bundle)
weighting.summary <- comparison_summary(weighting.bundle)
weighting.pairs <- merge(
  weighting.summary[weighting.summary$pair_weights=='uniform',c('case','backend','uniform_error','relative_error')],
  weighting.summary[weighting.summary$pair_weights=='inverse_squared',c('case','backend','uniform_error','relative_error')],
  by=c('case','backend'),suffixes=c('.uniform','.inverse'))

## ----weighting-overlay-saddle_graph-64, echo=FALSE----------------------------
if (surface.required) {
  if (!is.null(surface) && have.ivue) comparison_view(case,fits,weighting=TRUE, surface.alignments=surface$registrations)
  else documentation.full.view(case,fits,weighting=TRUE)
} else if (have.ivue) { comparison_view(case,fits,weighting=TRUE) } else {
documentation.full.view(case,fits,weighting=TRUE)
}

## ----weighting-overlay-paraboloid_graph-64, echo=FALSE------------------------
if (surface.required) {
  if (!is.null(surface) && have.ivue) comparison_view(case,fits,weighting=TRUE, surface.alignments=surface$registrations)
  else documentation.full.view(case,fits,weighting=TRUE)
} else if (have.ivue) { comparison_view(case,fits,weighting=TRUE) } else {
documentation.full.view(case,fits,weighting=TRUE)
}

## ----weighting-overlay-helix_graph-64, echo=FALSE-----------------------------
if (surface.required) {
  if (!is.null(surface) && have.ivue) comparison_view(case,fits,weighting=TRUE, surface.alignments=surface$registrations)
  else documentation.full.view(case,fits,weighting=TRUE)
} else if (have.ivue) { comparison_view(case,fits,weighting=TRUE) } else {
documentation.full.view(case,fits,weighting=TRUE)
}

## ----weighting-figures, results='asis', echo=FALSE----------------------------
for (id in names(weighting.bundle$cases)) {
  case <- weighting.bundle$cases[[id]]
  fits <- comparison_weighting_fits(weighting.bundle,id)
  surface.required <- grepl('^saddle', id)
  surface <- if (surface.required)
    comparison.surface.saved(case, fits) else NULL
  cat('\n\n### ',case$label,': uniform versus inverse-squared weights\n\n',sep='')
  med <- weighting.summary[weighting.summary$case==id,]
  med$backend <- toupper(med$backend)
  med$pair_weights <- ifelse(med$pair_weights=='uniform','Uniform','Inverse-squared')
  cat('Medians across five starting shapes:\n\n')
  cat(knitr::kable(med[,c('backend','pair_weights','successful','attempted',
                         'uniform_error','relative_error','seconds')],
    row.names=FALSE,digits=4,col.names=c('Backend','Weighting','Returned','Attempted',
      'Uniform error','Relative error','Fit seconds')),sep='\n')
  shown <- do.call(rbind,lapply(fits,`[[`,'row'))
  shown$backend <- toupper(shown$backend)
  shown$pair_weights <- ifelse(shown$pair_weights=='uniform','Uniform','Inverse-squared')
  cat('\n\nDisplayed runs (seed 1):\n\n')
  cat(knitr::kable(shown[,c('backend','pair_weights','uniform_error','relative_error',
                          'seconds','termination','status')], row.names=FALSE,digits=4,
    col.names=c('Backend','Weighting','Uniform error','Relative error','Fit seconds','Termination','Status')),sep='\n')
  child <- c(
    sprintf('```{r weighting-overlay-%s, echo=FALSE}',id),
    'if (surface.required) {',
    '  if (!is.null(surface) && have.ivue) comparison_view(case,fits,weighting=TRUE, surface.alignments=surface$registrations)',
    '  else documentation.full.view(case,fits,weighting=TRUE)',
    '} else if (have.ivue) { comparison_view(case,fits,weighting=TRUE) } else {',
    'documentation.full.view(case,fits,weighting=TRUE)','}','```')
  cat(knitr::knit_child(text=child,envir=environment(),quiet=TRUE),sep='\n')
}

