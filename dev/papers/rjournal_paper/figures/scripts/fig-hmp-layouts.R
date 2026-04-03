## ---- Figure: HMP microbial network layouts (Section 7.3) --------------------
##
## Three-panel 3D layout of the coarsened HMP graph, colored by CST.
## Shows: default preset | best preset | best local-search candidate.
## Adapted from grip-hmp-u01-coarse.Rmd.
##
## Output: figures/hmp-layouts.pdf

source("00-common.R")

## Load precomputed results (full pipeline is too expensive for a script)
res_path <- system.file(
  "extdata", "hmp_u01_gc_coarse", "vignette_results.rds",
  package = "grip"
)
if (!nzchar(res_path)) {
  stop("Precomputed HMP vignette results not found. ",
       "Run the HMP vignette precomputation first.")
}
hmp_results <- readRDS(res_path)

## Load the bundled graph for edge matrix + labels
data(hmp.u01.gc.coarse)

## Build edge matrix from adjacency list
edge_matrix_from_adj <- function(adj_list) {
  edges <- list()
  idx <- 0L
  for (u in seq_along(adj_list)) {
    nbrs <- adj_list[[u]]
    nbrs <- nbrs[nbrs > u]
    if (!length(nbrs)) next
    for (v in nbrs) {
      idx <- idx + 1L
      edges[[idx]] <- c(u, v)
    }
  }
  do.call(rbind, edges)
}

hmp_edges <- edge_matrix_from_adj(hmp.u01.gc.coarse$adj_list)
hmp_cst   <- hmp.u01.gc.coarse$vertex_data$cst

## CST colour palette
cst_levels <- sort(unique(hmp_cst))
cst_cols   <- setNames(grDevices::hcl.colors(length(cst_levels), "Dark 3"),
                        cst_levels)
hmp_cst_col <- cst_cols[hmp_cst]

## Extract layouts
lay_default <- hmp_results$layouts$preset$default
lay_tree    <- hmp_results$layouts$preset$tree
lay_torus   <- hmp_results$layouts$preset$torus

## Also get a top local candidate if available
top_local_names <- hmp_results$top_local_candidates
lay_local <- if (length(top_local_names) > 0) {
  hmp_results$layouts$local[[ top_local_names[1] ]]
} else {
  lay_torus
}
local_label <- if (length(top_local_names) > 0) {
  paste0("Local search (", top_local_names[1], ")")
} else {
  "torus preset"
}

## Draw figure
save_pdf("hmp-layouts.pdf", width = 14, height = 4.8, {
  op <- par(mfrow = c(1, 3), mar = c(1.2, 1.2, 3, 1.2), bg = "white")

  grip.plot(lay_default, hmp_edges,
            projection = "ortho",
            main       = "Default preset",
            vertex.col = hmp_cst_col,
            edge.col   = "gray90")

  grip.plot(lay_tree, hmp_edges,
            projection = "ortho",
            main       = "Tree preset",
            vertex.col = hmp_cst_col,
            edge.col   = "gray90")

  grip.plot(lay_local, hmp_edges,
            projection = "ortho",
            main       = local_label,
            vertex.col = hmp_cst_col,
            edge.col   = "gray90")

  ## Add legend to the last panel
  legend("bottomleft",
         legend = cst_levels,
         col    = cst_cols,
         pch    = 16,
         bty    = "n",
         cex    = 0.65,
         title  = "CST")

  par(op)
})

## Also print the preset comparison table for reference
cat("\nPreset comparison summary:\n")
print(hmp_results$preset_summary[, c("candidate",
                                      "sampled.stress.mean",
                                      "cluster.separation.mean",
                                      "score.composite")])

cat("\nWrote:", file.path(fig_dir, "hmp-layouts.pdf"), "\n")
