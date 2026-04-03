## ---- Figure: Multiscale trace of a mesh layout (Section 6) ------------------
##
## Four-panel trace showing coarsest init -> early -> mid -> final layout
## of a 6x6 mesh. Adapted from grip-examples.Rmd trace section.
##
## Output: figures/trace-mesh-6x6.pdf

source("00-common.R")

## Generate trace
mesh_edges <- edges.mesh(6, 6)
mesh_n     <- max(mesh_edges)

## Build the trace call; include diagnostics only if the installed version
## supports it (added after the initial CRAN-prep source tree).
trace_args <- list(
  edges        = mesh_edges,
  n            = mesh_n,
  dim          = 2,
  preset       = "mesh",
  rounds       = 12,
  final_rounds = 16,
  trace        = "level",
  seed         = 31
)
if ("diagnostics" %in% names(formals(grip.layout.trace))) {
  trace_args$diagnostics <- "light"
}
mesh_trace <- do.call(grip.layout.trace, trace_args)

## Select 4 representative frames: first, 1/3, 2/3, last
nf  <- length(mesh_trace$frames)
sel <- unique(c(1L,
                max(2L, round(nf / 3)),
                max(3L, round(2 * nf / 3)),
                nf))

## Build labels from trace metadata
labels <- paste0(
  mesh_trace$meta$phase[sel],
  " (frame ", sel, "/", nf, ")"
)

save_pdf("trace-mesh-6x6.pdf", width = 12, height = 3.2, {
  op <- par(mfrow = c(1, length(sel)), mar = c(1, 1, 2.8, 1), bg = "white")

  for (i in seq_along(sel)) {
    plot_trace_frame(
      mesh_trace$frames[[ sel[i] ]],
      mesh_edges,
      main = labels[i]
    )
  }

  par(op)
})

cat("Wrote:", file.path(fig_dir, "trace-mesh-6x6.pdf"), "\n")
