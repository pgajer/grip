#!/usr/bin/env Rscript
# Run from the repository root. Requires ivue with layer3D.axes(), camera.zup()
# and animate.frames(), plus rgl, htmlwidgets and Pandoc. No fits are recomputed.
options(rgl.useNULL = TRUE)
base <- "papers/grip-software-paper/reproducibility/precomputed"
pilot <- readRDS(file.path(base, "two-fidelity-saddle.rds"))
reference <- readRDS(file.path(base, "saddle-reference-diagnostics.rds"))
cloud <- reference$clouds[[as.character(pilot$representative)]]
Z <- cloud$aligned
titles <- c("Original saddle", "metric-MDS", "metric-MDS + edge-KK")
dir.create("output/readme-saddle", recursive = TRUE, showWarnings = FALSE)
out <- normalizePath("output/readme-saddle")

# Preserve the saved alignment, vertex colors and parameter-plane connectivity.
colors <- colorRampPalette(c("#173D65", "#86AFC4", "#D9B18B", "#8E4921"))(100)[
  pmin(100, pmax(1, 1 + floor((cloud$coords[, 1] + 1) * 49.5)))]
extent <- max(abs(do.call(rbind, Z))) * 1.20
limits <- matrix(rep(c(-extent, extent), each = 3), 3)
corners <- as.matrix(expand.grid(limits[1, ], limits[2, ], limits[3, ]))
common <- ivue::layer3D.callback(function(ctx) {
  rgl::points3d(corners, alpha = 0, size = 1)
})
mesh <- ivue::layer3D.mesh(cloud$triangles, col = "gray75", alpha = .23,
  edge.col = "gray45", edge.alpha = .23, edge.width = .65)
axes <- ivue::layer3D.axes(limits = limits, head.length = .04,
  head.angle = pi/8, width = 2, cex = 1.2)
camera <- ivue::camera.zup(elevation = 20, turn = -135, fov = 0, zoom = .5)
views <- lapply(seq_along(Z), function(j) {
  w <- ivue::plot3D.plain(Z[[j]], col = colors, point.type = "sphere",
    sphere.radius = .009, axes = FALSE, aspect = "equal",
    layers = list(mesh, axes, common), camera = camera,
    width = 600L, height = 510L)
  w$elementId <- paste0("saddle-view-", j)
  w
})

# ivue currently animates point/edge coordinates, not filled faces. Its player
# supplies the single timeline; the adapter below applies one camera matrix to
# all three intact mesh scenes on each tick (including scrubbing and reset).
fps <- 20
count <- 240L
angles <- (seq_len(count) - 1) * 360 / count
frames <- lapply(angles * pi/180, function(a) matrix(c(cos(a), sin(a), 0), 1))
timeline <- ivue::animate.frames(frames, fps = fps, max.frames = NULL,
  labels = sprintf("Rotation: %05.1f°", angles), col = "transparent",
  width = 64L, height = 64L)
timeline.id <- timeline$elementId
player.id <- timeline$x$players[[1]]
matrices <- lapply(angles, function(a) as.vector(
  ivue::camera.zup(elevation = 20, turn = -135 + a, fov = 0, zoom = .5)$userMatrix))
config <- jsonlite::toJSON(list(player = player.id,
  scenes = vapply(views, `[[`, "", "elementId"), matrices = matrices),
  auto_unbox = TRUE, digits = 15)
adapter <- paste0("(function() { const config = ", config, ";
  function connect() {
    const player = document.getElementById(config.player);
    const scenes = config.scenes.map(id => document.getElementById(id)?.rglinstance);
    if (!player?.rgltimer || scenes.some(s => !s)) {
      setTimeout(connect, 50); return;
    }
    const timer = player.rgltimer, tick = timer.Tick;
    function setFrame(value) {
      const frame = ((Math.round(value) % config.matrices.length) + config.matrices.length)
        % config.matrices.length;
      scenes.forEach(scene => {
        scene.getObj(scene.scene.rootSubscene).par3d.userMatrix.load(config.matrices[frame]);
        scene.drawScene();
      });
      document.body.dataset.frame = frame;
    }
    timer.Tick = function() { tick.call(this); setFrame(this.value); };
    window.saddlePreview = {timer, scenes, setFrame: function(i) {
      timer.value = i; timer.Tick();
    }};
    timer.Tick();
    document.body.dataset.ready = 'true';
    if (!new URLSearchParams(location.search).has('paused')) timer.PlayButton.click();
  }
  connect();
})();")
css <- paste0("
  body {margin:0; padding:30px 24px; background:white; color:#26323c;
    font:16px/1.5 system-ui,sans-serif;}
  main {max-width:1800px; margin:auto;}
  h1 {font-size:25px; font-weight:600; margin:0 0 8px;}
  p {margin:0 0 20px; color:#56616b;}
  #figure {display:grid; grid-template-columns:repeat(3,minmax(0,1fr)); background:white;}
  h2 {text-align:center; font-size:21px; font-weight:500; margin:15px 0 0;}
  .rglWebGL {max-width:100%;}
  #", timeline.id, " {position:absolute!important; width:64px!important;
    height:64px!important; left:-10000px;}
  #", player.id, " {display:flex; align-items:center; justify-content:center;
    flex-wrap:wrap; gap:8px; margin:20px 0;}
  input[type=button] {padding:7px 12px; font:inherit; cursor:pointer;}
  input[type=range] {width:300px;} output {min-width:145px;}
  .note {text-align:center; font-size:14px;}
")
page <- htmltools::tags$html(
  htmltools::tags$head(htmltools::tags$title("Saddle comparison — synchronized rotation"),
    htmltools::tags$meta(charset = "utf-8"),
    htmltools::tags$style(htmltools::HTML(css))),
  htmltools::tags$body(htmltools::tags$main(
    htmltools::tags$h1("Saddle configurations"),
    htmltools::tags$p("A shared z-axis rotation · 12 seconds per revolution · common scale and triangulation"),
    htmltools::tags$div(id = "figure", lapply(seq_along(views), function(j)
      htmltools::tags$section(htmltools::tags$h2(titles[j]), views[[j]]))),
    timeline,
    htmltools::tags$p(class = "note", "Pause or drag the slider to inspect the same angle in all three views."),
    htmltools::tags$script(htmltools::HTML(adapter)))))
raw <- file.path(out, "saddle-rotation-unbundled.html")
htmltools::save_html(page, raw, libdir = "saddle-libs")
stopifnot(rmarkdown::pandoc_available())
wrapped <- file.path(out, "saddle-rotation-raw.md")
template <- file.path(out, "saddle-template.html")
writeLines(c("```{=html}", readLines(raw, warn = FALSE), "```"), wrapped)
writeLines("$body$", template)
rmarkdown::pandoc_convert(wrapped, from = "markdown", to = "html",
  output = file.path(out, "saddle-rotation.html"),
  options = c("--standalone", "--embed-resources", "--template", template))
writeLines(c("Generated by: Rscript tools/pkg/render-saddle-preview.R",
  paste("ivue:", as.character(packageVersion("ivue"))),
  paste("rgl:", as.character(packageVersion("rgl"))),
  paste("Input MD5:", paste(tools::md5sum(file.path(base, c(
    "two-fidelity-saddle.rds", "saddle-reference-diagnostics.rds"))), collapse = " ")),
  "240 frames; 20 fps; 12 seconds per revolution; shared z-up camera matrices.",
  "Saved aligned coordinates and fixed triangulation; no refitting or deformation."),
  file.path(out, "provenance.txt"))
message("Wrote ", file.path(out, "saddle-rotation.html"))
