# Minimal executable example retained in Supplement S3.
# Run with library(grip); source("edge-kk-workflow.R"), or load the named
# chunks with knitr::read_chunk() when rendering the supplement.

## ---- geodesic-workflow ----
workflow.edges <- edges.mesh(4, 4)
workflow.grid <- as.matrix(expand.grid(
  x = seq(-1, 1, length.out = 4),
  y = seq(-1, 1, length.out = 4)
))
initial <- cbind(
  1.45 * workflow.grid[, 1] + 0.22 * workflow.grid[, 2]^2,
  0.55 * workflow.grid[, 2] + 0.10 * workflow.grid[, 1]
)

prepared <- prepare.edge.kk(workflow.edges, n = 16)
before <- score.gmds(initial, prepared = prepared)
fit <- edge.kk(
  coords = initial,
  prepared = prepared,
  density_mix_schedule = 1,
  max_iter = 8,
  return_trace = FALSE
)
after <- score.gmds(fit$coords, prepared = prepared)

## ---- edge-kk-workflow ----
display.coords <- function(coords) {
  centered <- sweep(coords, 2, colMeans(coords), "-", check.margin = FALSE)
  centered / sqrt(mean(rowSums(centered^2)))
}

workflow.layouts <- list(
  "Initial mesh coordinates" = display.coords(initial),
  "After edge-KK refinement" = display.coords(fit$coords)
)
workflow.scores <- c(before$edge.rel.rmse, after$edge.rel.rmse)
workflow.range <- range(do.call(rbind, workflow.layouts))
workflow.pad <- 0.12 * diff(workflow.range)

op <- par(mfrow = c(1, 2), mar = c(1, 1, 3.1, 1))
for (i in seq_along(workflow.layouts)) {
  coords <- workflow.layouts[[i]]
  plot(
    coords[, 1], coords[, 2], type = "n", asp = 1, axes = FALSE,
    xlab = "", ylab = "",
    xlim = workflow.range + c(-workflow.pad, workflow.pad),
    ylim = workflow.range + c(-workflow.pad, workflow.pad),
    main = sprintf("%s\nEdge relative RMSE = %.3f",
                   names(workflow.layouts)[i], workflow.scores[i])
  )
  segments(
    coords[workflow.edges[, 1], 1], coords[workflow.edges[, 1], 2],
    coords[workflow.edges[, 2], 1], coords[workflow.edges[, 2], 2],
    col = "gray82", lwd = 1
  )
  points(coords[, 1], coords[, 2], pch = 16, cex = 0.55,
         col = if (i == 1L) "#666666" else "#1F3B73")
}
par(op)
