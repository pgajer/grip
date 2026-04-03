## ---- Figure: Zachary karate club 3-step workflow (Section 7.2) --------------
##
## Shows the final tuned 3D layout of the karate club, colored by faction.
## The three-step workflow (preset shortlist -> local search -> metadata-aware
## comparison) is demonstrated in the paper text; this script produces the
## visual result.
##
## Output: figures/karate-workflow.pdf

source("00-common.R")

## Load karate club data
karate_edges <- as.matrix(
  read.csv(system.file("extdata", "karate-club-edges.csv", package = "grip"))
)
karate_n <- max(karate_edges)
labels_df <- read.csv(
  system.file("extdata", "karate-club-membership.csv", package = "grip")
)
karate_club <- labels_df$club[order(labels_df$vertex)]

## Step 1: Preset shortlisting
karate_presets <- grip.compare.layouts(
  edges      = karate_edges,
  n          = karate_n,
  dim        = 3,
  candidates = c("default", "tree", "mesh"),
  seeds      = 1:3,
  sample.size.stress = 1000L,
  edge.crossings     = "never",
  return.layouts     = TRUE
)

cat("Preset summary:\n")
print(karate_presets$summary[, c("candidate", "score.composite")])

## Step 2: Local parameter search around the best preset
best_preset <- karate_presets$summary[1, , drop = FALSE]
karate_local <- grip.compare.layouts(
  edges = karate_edges,
  n     = karate_n,
  dim   = 3,
  search = list(
    candidate.prefix = "karate.local",
    placement        = best_preset$placement[[1L]],
    rounds           = best_preset$rounds[[1L]],
    final_rounds     = c(128L, 160L, 192L),
    num_nbrs         = c(4L, 8L, 12L),
    repulsion_factor = c(0, 0.5)
  ),
  seeds = 1:2,
  sample.size.stress = 1000L,
  edge.crossings     = "never",
  return.layouts     = TRUE
)

cat("\nLocal search top 3:\n")
print(head(karate_local$summary[, c("candidate", "score.composite")], 3))

## Step 3: Metadata-aware comparison (best local vs default)
karate_meta <- grip.compare.layouts(
  edges      = karate_edges,
  n          = karate_n,
  dim        = 3,
  candidates = list(
    default = NULL,
    tuned   = grip.params.from.summary(
      karate_local$summary[1, , drop = FALSE]
    )
  ),
  clusters           = karate_club,
  seeds              = 1:3,
  sample.size.stress = 1000L,
  edge.crossings     = "never",
  return.layouts     = TRUE
)

cat("\nMetadata-aware comparison:\n")
print(karate_meta$summary[, c("candidate", "score.composite",
                               "cluster.separation.mean")])

## Extract the best layout for plotting
best_name   <- karate_meta$summary$candidate[1]
best_coords <- karate_meta$layouts[[best_name]][[1]]

## Faction colors
faction_cols <- ifelse(karate_club == 1, "#1F3B73", "#D35400")

## Draw the figure: two panels (default vs tuned)
default_coords <- karate_meta$layouts[["default"]][[1]]
tuned_coords   <- karate_meta$layouts[["tuned"]][[1]]

save_pdf("karate-workflow.pdf", width = 10, height = 4.5, {
  op <- par(mfrow = c(1, 2), mar = c(1.2, 1.2, 3, 1.2), bg = "white")

  grip.plot(
    default_coords, karate_edges,
    projection = "ortho",
    main       = "Default GRIP",
    vertex.col = faction_cols,
    edge.col   = "gray82"
  )

  grip.plot(
    tuned_coords, karate_edges,
    projection = "ortho",
    main       = "Tuned GRIP (3-step workflow)",
    vertex.col = faction_cols,
    edge.col   = "gray82"
  )

  par(op)
})

cat("Wrote:", file.path(fig_dir, "karate-workflow.pdf"), "\n")
