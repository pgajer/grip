# Complete documentation benchmarks

These are the original frozen metric-MDS and sparse-MDS benchmark bundles,
including every saved layout and independent evaluation target. No fits were
rerun when the installed documentation bundles were reduced for CRAN.

The installed bundles retain every attempted-fit record, failure, score, timing
and metadata entry, plus seed-1 layouts used by the guides (budget 100 for the
full-MDS comparison). Sparse evaluation-source identities and modes remain;
the full target matrices are retained here. The separate higher-dimensional
bundle is unchanged.

Run `Rscript tools/pkg/compact-documentation-data.R` from the repository root
to regenerate the installed bundles, or add `--check` to verify them. Rendering
uses saved scores and coordinates; changing the storage does not change the
scientific comparisons. Load these complete bundles to recompute scores or
inspect other starts and budgets.
