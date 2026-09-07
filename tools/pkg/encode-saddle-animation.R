#!/usr/bin/env Rscript
# Capture the filled WebGL meshes rather than use ivue's point/edge-only GIF
# renderer. The captures are frames from the shared ivue animation timeline.
overlay <- "--overlay-reference" %in% commandArgs(trailingOnly = TRUE)
out <- if (overlay) "output/s4-3c-rotation" else "output/readme-saddle"
paths <- file.path(out, "frames", sprintf("frame-%03d.png", 0:239))
stopifnot(all(file.exists(paths)))
frames <- magick::image_scale(magick::image_read(paths), if (overlay) "720x" else "1350x")
gif <- magick::image_animate(frames, fps = 20, loop = 0, optimize = TRUE)
magick::image_write(gif, file.path(out, "saddle-rotation.gif"))
dir.create("man/figures", recursive = TRUE, showWarnings = FALSE)
readme.file <- if (overlay) "man/figures/readme-saddle-overlay-rotation.gif" else
  "man/figures/readme-saddle-rotation.gif"
stopifnot(file.copy(file.path(out, "saddle-rotation.gif"),
  readme.file, overwrite = TRUE))
message("Wrote ", file.path(out, "saddle-rotation.gif"))
