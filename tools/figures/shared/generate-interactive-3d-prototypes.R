#!/usr/bin/env Rscript

options(rgl.useNULL = TRUE)

suppressPackageStartupMessages({
  library(pkgload)
  library(rgl)
  library(htmlwidgets)
  library(plotly)
})

repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

pkgload::load_all(repo_root, quiet = TRUE, export_all = FALSE, helpers = FALSE)

out_dir <- file.path(repo_root, "dev", "manual", "interactive-prototypes")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

edges_to_segments <- function(edges, coords) {
  x <- y <- z <- numeric(0)
  for (i in seq_len(nrow(edges))) {
    e <- edges[i, ]
    x <- c(x, coords[e[1], 1], coords[e[2], 1], NA_real_)
    y <- c(y, coords[e[1], 2], coords[e[2], 2], NA_real_)
    z <- c(z, coords[e[1], 3], coords[e[2], 3], NA_real_)
  }
  list(x = x, y = y, z = z)
}

edges <- edges.torus(8, 8)
n <- max(edges)
coords <- grip.layout(edges, n = n, dim = 3, preset = "torus", seed = 21)
segs <- edges_to_segments(edges, coords)

rgl::open3d(useNULL = TRUE)
rgl::bg3d(color = "white")
rgl::segments3d(
  segs$x, segs$y, segs$z,
  color = "#B6BCC6",
  alpha = 0.45,
  lwd = 1
)
rgl::points3d(
  coords[, 1], coords[, 2], coords[, 3],
  color = "#1F3B73",
  size = 5
)
rgl::view3d(theta = 35, phi = 22, zoom = 0.8)
rgl_widget <- rgl::rglwidget(width = 900, height = 650)
rgl_path <- file.path(out_dir, "torus-rglwidget.html")
htmlwidgets::saveWidget(rgl_widget, file = rgl_path, selfcontained = TRUE)
rgl::close3d()

plotly_obj <- plot_ly()
plotly_obj <- add_trace(
  plotly_obj,
  x = segs$x,
  y = segs$y,
  z = segs$z,
  type = "scatter3d",
  mode = "lines",
  line = list(color = "rgba(182,188,198,0.55)", width = 2),
  hoverinfo = "none",
  showlegend = FALSE
)
plotly_obj <- add_trace(
  plotly_obj,
  x = coords[, 1],
  y = coords[, 2],
  z = coords[, 3],
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 3.5, color = "#1F3B73"),
  hoverinfo = "none",
  showlegend = FALSE
)
plotly_obj <- layout(
  plotly_obj,
  scene = list(
    xaxis = list(visible = FALSE),
    yaxis = list(visible = FALSE),
    zaxis = list(visible = FALSE),
    camera = list(eye = list(x = 1.45, y = 1.1, z = 0.8))
  ),
  paper_bgcolor = "white",
  plot_bgcolor = "white",
  margin = list(l = 0, r = 0, b = 0, t = 0)
)
plotly_path <- file.path(out_dir, "torus-plotly.html")
htmlwidgets::saveWidget(plotly_obj, file = plotly_path, selfcontained = TRUE)

sizes <- data.frame(
  artifact = c("torus-rglwidget.html", "torus-plotly.html"),
  bytes = c(file.info(rgl_path)$size, file.info(plotly_path)$size)
)
write.csv(sizes, file.path(out_dir, "artifact-sizes.csv"), row.names = FALSE)

cat("Wrote:\n")
cat(" - ", rgl_path, "\n", sep = "")
cat(" - ", plotly_path, "\n", sep = "")
cat(" - ", file.path(out_dir, "artifact-sizes.csv"), "\n", sep = "")
