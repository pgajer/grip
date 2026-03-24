# HMP/U01 coarsened graph provenance

This directory contains raw bundled artifacts for the coarsened HMP/U01
giant-component graph shipped with the `grip` package.

Source workflow:

- upstream project: `ZB/chm_paper`
- upstream dataset: publicly available HMP+U01 16S amplicon data
- upstream branch: `>=1%` feature screen with `relative_abundance_pca`
- selected graph parameter: `k = 3`
- upstream usable sample count on this branch: `6521`
- original graph giant component size: `6474`
- coarsening trajectory: `6474 -> 3449 -> 1828`

Bundled files:

- `coarse_graph_edges.tsv.gz`: weighted undirected edge list on the coarsened graph
- `coarse_vertex_labels.tsv.gz`: coarse-vertex labels and aggregated metadata
- `coarse_membership.tsv.gz`: exact mapping from original giant-component samples to supernodes
- `coarsening_rounds.tsv`: vertex counts after each coarsening round

The packaged data object `hmp.u01.gc.coarse` is generated from these raw files
by `data-raw/hmp_u01_gc_coarse.R`.
