## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(fig.retina=1, dpi=72, 
  collapse = TRUE,
  comment = "#>",
  fig.width = 8.8,
  fig.height = 5.2
)

## -----------------------------------------------------------------------------
library(grip)
source(system.file('scripts','vignette-views.R',package='grip'),local=TRUE)

## -----------------------------------------------------------------------------
mesh.edges <- edges.sierpinski.triangle(4)

mesh.trace <- trace.grip(
  mesh.edges,
  n = max(mesh.edges),
  dim = 2,
  preset = "carpet",
  rounds = 48,
  final.rounds = 96,
  trace = "level",
  diagnostics = "light",
  seed = 1
)

knitr::kable(mesh.trace$meta)
knitr::kable(mesh.trace$diagnostics[, c(
  "frame",
  "edge.length.cv",
  "median.edge.length",
  "sampled.nonedge.sep.ratio"
)], digits = 3)

## ----trace-view, echo=FALSE---------------------------------------------------
if (!documentation.interactive()) {
  documentation.trace(mesh.trace, mesh.edges)
} else htmltools::tagList(lapply(seq_along(mesh.trace$frames), function(i) {
  z <- mesh.trace$frames[[i]]
  active <- stats::complete.cases(z)
  e <- mesh.edges[active[mesh.edges[,1]] & active[mesh.edges[,2]],,drop=FALSE]
  ids <- integer(nrow(z)); ids[active] <- seq_len(sum(active))
  e <- matrix(ids[e],ncol=2)
  z <- z[active,,drop=FALSE]
  if (ncol(z)==2L) z <- cbind(z,0)
  htmltools::tagList(
    htmltools::tags$h3(paste('Recorded stage',i,'of',length(mesh.trace$frames))),
    ivue::plot3D.plain(z,layers=list(ivue::layer3D.edges(e)),
      axes=FALSE,aspect='equal',width=900,height=420,
      description=paste('Stage',i,':',sum(active),'active vertices. Original planar coordinates.'), controls = FALSE)
  )
}))

## -----------------------------------------------------------------------------
mesh.trace.round <- trace.grip(
  mesh.edges,
  n = max(mesh.edges),
  dim = 2,
  preset = "carpet",
  rounds = 6,
  final.rounds = 8,
  trace = "round",
  trace.every = 4,
  diagnostics = "light",
  seed = 1
)

data.frame(
  trace_mode = c("level", "round"),
  n.frames = c(length(mesh.trace$frames), length(mesh.trace.round$frames))
)

## ----diagnostic-curve, echo=FALSE---------------------------------------------
plot(mesh.trace$diagnostics$frame,mesh.trace$diagnostics$edge.length.cv,
  type='b',pch=16,col='#1769AA',xlab='Recorded frame',ylab='Edge-length variation')

