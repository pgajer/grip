# HMP-only graph provenance

This directory contains the exact graph and sample-level inputs used for the
Human Microbiome Project example in the `grip` R Journal paper. The graph was
rebuilt on 2026-08-22 by `data-raw/hmp_gc.R`.

## Cohort boundary

The paper graph retains only samples that satisfy both of these positive
criteria before any feature screening:

- `Project` begins with `HMP`;
- `16S_Platform` is exactly `Illumina`.

These criteria identify 4,411 eligible samples. Every vertex distributed here
is therefore an explicitly identified HMP Illumina 16S rRNA amplicon sample.

## Representation and graph construction

Features detected in at least 1% of eligible samples are retained; 169 features
pass this screen. Samples must retain at least 500 reads and at least 25% of
their original reads after screening. Counts are converted to relative
abundance, represented by the first 10 principal components, and used to build
a symmetric 3-nearest-neighbor graph with `dgraphs::create.sknn.graph()`.
Symmetrization uses the union rule: an undirected edge is present when either
sample selects the other as a neighbor. No component-connecting or graph-
coarsening procedure is applied.

The largest connected component contains 4,391 vertices and 9,067 edges.
Twenty eligible samples lie in three smaller components and are not used in
the paper example.

## Distributed files

- `graph_edges.tsv.gz`: undirected edges and PCA-space edge lengths;
- `vertex_metadata.tsv.gz`: sample identifiers, HMP/platform fields, and
  retained-read diagnostics for every graph vertex;
- `feature_manifest.tsv.gz`: the complete feature screen and retention flag;
- `graph_summary.tsv`: cohort, graph, and preprocessing counts.

The raw upstream abundance and metadata tables are not redistributed, and no
public acquisition procedure for those exact tables is currently supplied.
Consequently, this directory reproduces the paper from the exact final graph
but does not independently regenerate that graph from primary abundance data.
When the upstream tables are available separately, their locations must be
supplied with `GRIP_HMP_METADATA_TSV` and `GRIP_HMP_FEATURE_MATRIX_TSV`;
`data-raw/hmp_gc.R` then records and executes the complete filtering and graph-
construction procedure. The distributed edge, vertex, and feature-manifest
files remain sufficient to inspect the cohort boundary and reconstruct the
exact final graph used by the paper.
