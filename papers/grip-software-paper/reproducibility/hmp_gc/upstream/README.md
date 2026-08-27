# UMB-HMP-only upstream inputs

Version: 1.0.0

Prepared: 2026-08-24

This directory contains the minimal inputs from the Ravel-led University of
Maryland Baltimore Human Microbiome Project (UMB-HMP) needed to regenerate the
graph used in the `grip` R Journal paper. UMB-HMP was a longitudinal vaginal
microbiome HMP demonstration project, not the NIH HMP healthy-reference cohort.
The source working tables also contained a related U01 cohort and extensive
clinical metadata. Those rows and fields are not included here.

## Files

- `hmp_illumina_metadata.tsv.gz`: 4,411 UMB-HMP Illumina sample identifiers and
  the three technical cohort fields used by the graph builder;
- `hmp_illumina_feature_counts.tsv.gz`: microbial feature counts for the same
  4,411 samples, in identical row order;
- `MANIFEST.tsv`: descriptions, dimensions, byte sizes, and SHA-256 checksums;
- `SHA256SUMS`: checksums in standard `shasum -c` format;
- `DATA_USE_AND_LICENSE.md`: provenance, disclosure boundary, and reuse terms;
- `zenodo-metadata.json`: metadata prepared for a versioned Zenodo deposit.

The metadata table contains no clinical variables, participant identifiers,
host sequence, or U01 records. The feature table contains microbial abundance
counts only.

## Rebuilding the graph

From the parent `reproducibility` directory, run:

```sh
Rscript scripts/build-hmp-only-graph.R
```

The builder reads these two files by default and writes regenerated artifacts
under `generated/hmp_gc/`. Alternative paths and output locations can be
supplied with `GRIP_HMP_METADATA_TSV`, `GRIP_HMP_FEATURE_MATRIX_TSV`, and
`GRIP_HMP_OUTPUT_DIR`.

The regenerated graph has 4,391 vertices and 9,067 edges. Validation on
2026-08-24 confirmed exact agreement with the paper artifacts for vertex IDs,
edge endpoints, edge weights, and the feature-screening manifest. The
regenerated `source_rows` summary is 4,411 rather than 6,536 because this
deposit begins after the documented UMB-HMP/Illumina cohort extraction; all
graph-defining quantities are unchanged.

## Extraction record

`scripts/prepare-hmp-upstream-deposit.R` creates this directory from the
historical combined working tables. It applies the positive UMB-HMP and Illumina
filters, selects only the four technical metadata columns, aligns the feature
rows by `canonical_sample_id`, validates the counts, and writes the checksums
and manifest. The combined tables are not needed to rebuild or audit the graph.
