## ---- Master script: generate all Paper 1 figures ----------------------------
##
## Usage:
##   cd dev/papers/rjournal_paper/figures/scripts
##   Rscript run-all.R
##
## Or from R:
##   setwd("dev/papers/rjournal_paper/figures/scripts")
##   source("run-all.R")
##
## All figures are written to ../figures/ as PDF files.

scripts <- c(
  "fig-architecture.R",
  "fig-misf-hierarchy.R",
  "fig-trace-mesh.R",
  "fig-preset-comparison.R",
  "fig-triptych-saddle.R",
  "fig-2d-vs-3d.R",
  "fig-karate-workflow.R",
  "fig-hmp-layouts.R",
  "fig-comparison-panel.R"
)

cat("=== Generating all Paper 1 figures ===\n\n")

for (s in scripts) {
  cat("--- Running:", s, "---\n")
  tryCatch(
    source(s, local = new.env(parent = globalenv())),
    error = function(e) {
      cat("  ERROR in", s, ":", conditionMessage(e), "\n")
    }
  )
  cat("\n")
}

cat("=== Done. Check ../figures/ for output PDFs. ===\n")
