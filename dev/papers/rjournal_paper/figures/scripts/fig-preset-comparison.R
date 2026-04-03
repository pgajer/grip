## ---- Figure: Default vs preset-tuned layout (Section 5) ---------------------
##
## Side-by-side comparison of grip.layout() with no preset vs with the
## "carpet" preset on a Sierpinski carpet level-3. Shows the value of
## validated parameter presets.
##
## Output: figures/preset-comparison.pdf

source("00-common.R")

## Generate a Sierpinski carpet (level 3 for visual clarity)
carpet_edges <- edges.sierpinski.carpet(3)
carpet_n     <- max(carpet_edges)

## Default layout (no preset)
coords_default <- grip.layout(
  carpet_edges,
  n    = carpet_n,
  dim  = 2,
  seed = 1
)

## Preset-tuned layout
coords_preset <- grip.layout(
  carpet_edges,
  n      = carpet_n,
  dim    = 2,
  preset = "carpet",
  seed   = 1
)

## Score both
score_default <- grip.score.layout(coords_default, edges = carpet_edges, n = carpet_n)
score_preset  <- grip.score.layout(coords_preset,  edges = carpet_edges, n = carpet_n)

cat("\n--- Preset comparison scores ---\n")
cat("Default:  stress =", round(score_default$sampled.stress, 4),
    "  edge-length CV =", round(score_default$edge.length.cv, 4), "\n")
cat("Carpet preset: stress =", round(score_preset$sampled.stress, 4),
    "  edge-length CV =", round(score_preset$edge.length.cv, 4), "\n")

## Draw figure
save_pdf("preset-comparison.pdf", width = 10, height = 4.5, {
  op <- par(mfrow = c(1, 2), mar = c(1.2, 1.2, 3, 1.2), bg = "white")

  grip.plot(
    coords_default, carpet_edges,
    main       = paste0("Default (stress = ",
                        round(score_default$sampled.stress, 3), ")"),
    vertex.col = "black",
    edge.col   = "gray82",
    pch        = 16,
    cex        = 0.45
  )

  grip.plot(
    coords_preset, carpet_edges,
    main       = paste0("Carpet preset (stress = ",
                        round(score_preset$sampled.stress, 3), ")"),
    vertex.col = "#1F3B73",
    edge.col   = "gray82",
    pch        = 16,
    cex        = 0.45
  )

  par(op)
})

cat("Wrote:", file.path(fig_dir, "preset-comparison.pdf"), "\n")
