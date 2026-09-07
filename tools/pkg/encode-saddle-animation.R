#!/usr/bin/env Rscript
# Capture the filled WebGL meshes rather than use ivue's point/edge-only GIF
# renderer. The captures are frames from the shared ivue animation timeline.
out <- "output/readme-saddle"
paths <- file.path(out, "frames", sprintf("frame-%03d.png", 0:239))
stopifnot(all(file.exists(paths)))
frames <- magick::image_scale(magick::image_read(paths), "1350x")
gif <- magick::image_animate(frames, fps = 20, loop = 0, optimize = TRUE)
magick::image_write(gif, file.path(out, "saddle-rotation.gif"))
message("Wrote ", file.path(out, "saddle-rotation.gif"))
