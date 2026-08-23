# Reproducing the `grip` software-paper results

This directory makes the precomputed components of the paper inspectable and
reproduces all reported results from the supplied final graph. The manuscript
evaluates the smaller examples directly.

## Requirements

- R and `grip` 0.2.0 or later;
- `dgraphs` to rebuild the HMP graph from the upstream count and metadata
  tables;
- `igraph` and `graphlayouts` for the cross-package benchmark;
- the packages listed in the manuscript-level `_Rpackages.txt`.

## Included artifacts

- `precomputed/vs_alternatives/benchmark_results.rds`: cross-package layouts,
  scores, timings, and generation-session information;
- `BENCHMARK_PROVENANCE.md`: timing boundary, repeat policy, hardware,
  software versions, and the benchmark medians reported in the paper;
- `hmp_gc/`: the exact HMP-only Illumina 16S edge list, vertex metadata,
  feature-screening manifest, graph summary, and provenance record;
- `SHA256SUMS`: checksums for every supplied binary or tabular input and
  precomputed result.

The manuscript uses the full 4,391-vertex giant component of the HMP-only
graph.

## HMP-only graph

`scripts/build-hmp-only-graph.R` documents and executes the cohort filter,
feature screening, and graph construction when the two combined upstream
tables are supplied separately. It retains only rows explicitly labeled as HMP
and Illumina, then constructs the symmetric 3-nearest-neighbor graph described
in `hmp_gc/PROVENANCE.md`:

```sh
GRIP_HMP_METADATA_TSV=/path/to/hmp_analysis_metadata.tsv \
GRIP_HMP_FEATURE_MATRIX_TSV=/path/to/hmp_feature_matrix.tsv \
  Rscript scripts/build-hmp-only-graph.R
```

The raw upstream tables are not redistributed, and no public acquisition
procedure for those exact tables is currently supplied. The supplement
therefore does not independently regenerate the HMP graph from primary
abundance data. The supplied HMP-only edge, sample, and feature-manifest files
do reproduce the exact final graph used by the paper and allow the inclusion
boundary and every downstream reported result to be audited directly.

## Cross-package benchmark

From this directory, run:

```sh
GRIP_VS_ALTERNATIVES_OUTPUT=generated/benchmark_results.rds \
  Rscript scripts/precompute-vs-alternatives.R
```

The benchmark resets the layout seed before each of five calls and reports the
median elapsed time and interquartile range. Timing covers layout generation
only; scoring and garbage collection are excluded. Every method receives the
same unweighted HMP topology, and the common score is sampled hop-distance
stress. Elapsed times remain machine-dependent; `BENCHMARK_PROVENANCE.md`
records the machine used for the supplied artifact.

## Fast article build

The manuscript reads the supplied file under `precomputed/` so reviewers do
not need to rerun the heavier benchmark. The scripts document the graph-
construction boundary and execute the derivation of the benchmark artifact
from the supplied final HMP graph.
