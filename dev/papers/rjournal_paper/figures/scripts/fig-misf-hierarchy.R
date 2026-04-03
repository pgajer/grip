## ---- Figure: MIS filtration hierarchy illustration (Section 2) ---------------
##
## Shows the coarse-to-fine MIS hierarchy on a small mesh:
##   Panel 1: Coarsest level (initial seed vertices)
##   Panel 2: Intermediate level
##   Panel 3: Full graph with all vertices placed
##
## Uses grip.layout.trace to extract the hierarchy from the trace frames:
## at each frame, vertices with non-NA coordinates are the "active" set.
## Falls back to grip.build.misf if available.
##
## Output: figures/misf-hierarchy.pdf

source("00-common.R")

## Use a small mesh so the hierarchy is visually clear
mesh_edges <- edges.mesh(6, 6)
mesh_n     <- max(mesh_edges)

## Compute a reference layout for consistent vertex positions across panels
coords <- grip.layout(mesh_edges, n = mesh_n, dim = 2, preset = "mesh", seed = 1)

if (exists("grip.build.misf", mode = "function")) {
  ## ---- Primary path: use grip.build.misf directly ----
  misf <- grip.build.misf(mesh_edges, n = mesh_n)

  vertex_depth <- misf$vertex_depth
  max_level    <- max(vertex_depth)

  show_levels <- unique(c(0L,
                          as.integer(floor(max_level / 2)),
                          max_level))

  panel_labels <- paste0("Level ", show_levels,
                         " (", sapply(show_levels, function(l) sum(vertex_depth >= l)),
                         " vertices)")

  active_sets <- lapply(show_levels, function(l) vertex_depth >= l)

} else {
  ## ---- Fallback: infer hierarchy from trace frames ----
  cat("grip.build.misf not found in installed package; using trace fallback.\n")

  mesh_trace <- grip.layout.trace(
    mesh_edges,
    n            = mesh_n,
    dim          = 2,
    preset       = "mesh",
    rounds       = 12,
    final_rounds = 16,
    trace        = "level",
    seed         = 31
  )

  nf <- length(mesh_trace$frames)

  ## Pick coarsest, mid, and final frames
  sel <- unique(c(1L,
                  max(2L, round(nf / 2)),
                  nf))

  active_sets <- lapply(sel, function(idx) {
    complete.cases(mesh_trace$frames[[idx]][, 1:2, drop = FALSE])
  })

  panel_labels <- paste0(
    c("Coarsest", "Intermediate", "Final")[seq_along(sel)],
    " (", sapply(active_sets, sum), " vertices)"
  )
}

## ---- Draw figure ----
save_pdf("misf-hierarchy.pdf", width = 12, height = 4, {
  op <- par(mfrow = c(1, length(active_sets)),
            mar = c(1, 1, 3, 1), bg = "white")

  for (i in seq_along(active_sets)) {
    active <- active_sets[[i]]

    plot(coords[, 1], coords[, 2],
         asp = 1, type = "n", axes = FALSE,
         xlab = "", ylab = "", main = panel_labels[i])

    ## All edges very faint
    segments(coords[mesh_edges[, 1], 1], coords[mesh_edges[, 1], 2],
             coords[mesh_edges[, 2], 1], coords[mesh_edges[, 2], 2],
             col = "gray92")

    ## Edges between active vertices
    active_edge_mask <- active[mesh_edges[, 1]] & active[mesh_edges[, 2]]
    if (any(active_edge_mask)) {
      ae <- mesh_edges[active_edge_mask, , drop = FALSE]
      segments(coords[ae[, 1], 1], coords[ae[, 1], 2],
               coords[ae[, 2], 1], coords[ae[, 2], 2],
               col = "gray70")
    }

    ## Inactive vertices as faint dots
    points(coords[!active, 1], coords[!active, 2],
           pch = 16, cex = 0.5, col = "gray85")

    ## Active vertices as solid dots
    points(coords[active, 1], coords[active, 2],
           pch = 16, cex = 1.0, col = "#1F3B73")
  }

  par(op)
})

cat("Wrote:", file.path(fig_dir, "misf-hierarchy.pdf"), "\n")
