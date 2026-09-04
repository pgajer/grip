#!/usr/bin/env Rscript
# Run from the grip repository root. IVUE_SOURCE can select a source checkout.
ivue_source <- Sys.getenv("IVUE_SOURCE")
if(nzchar(ivue_source)) pkgload::load_all(ivue_source,quiet=TRUE)
stopifnot(all(c("layer3D.mesh","layer3D.surface") %in% getNamespaceExports("ivue")))
paper <- normalizePath("papers/grip-software-paper")
rmarkdown::render(file.path(paper,"supplement/S4-interactive-saddle.Rmd"),
  output_dir=file.path(paper,"supplement"),
  intermediates_dir=file.path(paper,"build/supplement"),
  knit_root_dir=paper,quiet=TRUE)
