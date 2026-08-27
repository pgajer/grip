# UMB-HMP-only graph provenance

This directory contains the exact graph and sample-level inputs used for the
University of Maryland Baltimore Human Microbiome Project (UMB-HMP) example in
the `grip` R Journal paper. UMB-HMP was the Ravel-led longitudinal vaginal
microbiome demonstration project, not the NIH HMP healthy-reference cohort.
The parent study enrolled 135 nonpregnant women of reproductive age for daily
vaginal self-sampling over ten weeks (Ravel et al., 2013, DOI:
10.1186/2049-2618-1-29). Later work identifies the cohort as UMB-HMP and its
profiles as V3--V4 16S rRNA data (Lee et al., 2023, DOI:
10.1371/journal.pcbi.1011295). The graph was rebuilt on 2026-08-22 by
`data-raw/hmp_gc.R`; the independently runnable supplemental builder is
`../scripts/build-hmp-only-graph.R`.

## Cohort boundary

The paper graph retains only samples that satisfy both of these positive
criteria before any feature screening:

- `Project` begins with `HMP`;
- `16S_Platform` is exactly `Illumina`.

These criteria identify 4,411 eligible samples. Every vertex distributed here
is therefore an explicitly identified UMB-HMP Illumina 16S rRNA amplicon
sample. The cited publications establish the cohort and assay lineage; the
upstream tables in this supplement define the exact 4,411-profile extraction
used for the graph.

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

The `upstream/` subdirectory supplies the minimal UMB-HMP-only feature-count and
technical metadata tables required by the supplemental builder. These tables
contain 4,411 eligible samples and 231 input features; they exclude U01
records, clinical fields, participant identifiers, and host sequence. Their
manifest, checksums, provenance and licensing information are included in the
same directory. The historical combined working tables are neither required
nor distributed.

Running `scripts/build-hmp-only-graph.R` with its default paths reconstructs
the graph from these minimal tables. Validation against the graph supplied with
the manuscript confirmed identical vertex identifiers, edge endpoints, edge
weights, and feature-manifest rows. The rebuilt graph therefore preserves the
exact graph-defining outputs used in the paper; its `source_rows` summary is
4,411 rather than the 6,536 rows in the historical combined HMP/U01 working
table because the archive begins after the explicit UMB-HMP-only cohort extraction.
