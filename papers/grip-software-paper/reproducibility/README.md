# Reproducing the `grip` software-paper results

This directory makes the precomputed components of the paper inspectable and
reproduces all reported results, including reconstruction of the UMB-HMP-only graph
from the supplied minimal upstream tables. The manuscript evaluates the smaller
examples directly.

## Requirements

- R and `grip` 0.2.0 or later;
- `dgraphs` to rebuild the UMB-HMP graph from the upstream count and metadata
  tables;
- `igraph` and `graphlayouts` for the cross-package benchmark;
- the packages listed in the manuscript-level `_Rpackages.txt`.

## Included artifacts

- `precomputed/vs_alternatives/benchmark_results.rds`: cross-package layouts,
  scores, timings, and generation-session information;
- `BENCHMARK_PROVENANCE.md`: timing boundary, repeat policy, hardware,
  software versions, and the benchmark medians reported in the paper;
- `hmp_gc/`: the exact UMB-HMP-only Illumina 16S edge list, vertex metadata,
  feature-screening manifest, graph summary, and provenance record;
- `hmp_gc/upstream/`: the minimal UMB-HMP-only feature-count and technical
  metadata tables needed to reconstruct that graph, together with an archive
  manifest, checksums, provenance and licensing information;
- `SHA256SUMS`: checksums for every supplied binary or tabular input and
  precomputed result.

The manuscript uses the full 4,391-vertex giant component of the UMB-HMP-only
graph. UMB-HMP denotes the Ravel-led longitudinal vaginal microbiome HMP
demonstration project, not the NIH HMP healthy-reference cohort.

## Paper-to-artifact map

| Manuscript component | Supplied input | Regeneration or validation path |
|---|---|---|
| Fixed weighted graph comparison | `precomputed/vs_alternatives/benchmark_results.rds`, component `weighted_saddle` | `scripts/precompute-vs-alternatives.R` |
| UMB-HMP-only graph construction summary | `hmp_gc/` graph files and `hmp_gc/upstream/` input tables | `scripts/build-hmp-only-graph.R` |
| Repeated HMP runtime and shared-scoring benchmark | `precomputed/vs_alternatives/benchmark_results.rds`, component `hmp` | `scripts/precompute-vs-alternatives.R` |
| Benchmark hardware, software, timing boundary, and repeat policy | `BENCHMARK_PROVENANCE.md` and the RDS `benchmark_metadata` component | Recorded by `scripts/precompute-vs-alternatives.R` |
| Small-graph figures and tables | Evaluated directly from `grip-software-paper.Rmd` | Render the manuscript; no precomputed supplement input is used |

The table identifies the manuscript result that consumes each supplied
artifact. Figure and table numbering may change during editing, so the mapping
uses stable section descriptions and RDS component names.

## Software environment and validation

`BENCHMARK_PROVENANCE.md` records the released R version, operating system,
hardware, BLAS, package versions, timing boundary, and random-number policy used
for the supplied benchmark artifact. The RDS file also retains its generation
session metadata. The complete package list needed to render the article is in
the manuscript-level `_Rpackages.txt` file.

From this directory, verify the distributed inputs and artifacts with:

```sh
shasum -a 256 -c SHA256SUMS
```

The graph builder validates dimensions, identifiers, edge endpoints, edge
weights, and the feature manifest against the supplied HMP graph. The benchmark
script records the software session with its regenerated artifact.

## UMB-HMP-only graph

`scripts/build-hmp-only-graph.R` documents and executes the cohort filter,
feature screening, and graph construction from the supplied minimal UMB-HMP-only
tables. It retains only rows explicitly labeled as HMP and Illumina, then
constructs the symmetric 3-nearest-neighbor graph described in
`hmp_gc/PROVENANCE.md`:

```sh
Rscript scripts/build-hmp-only-graph.R
```

The minimal inputs exclude U01 records, clinical fields, participant
identifiers, and host sequence. They were extracted from the historical
combined working tables by `scripts/prepare-hmp-upstream-deposit.R`; the
combined tables are neither required nor distributed. A clean rebuild from the
minimal inputs reproduces the exact vertex identifiers, edge endpoints, edge
weights, and feature manifest used in the paper. Paths can still be overridden
with `GRIP_HMP_METADATA_TSV`, `GRIP_HMP_FEATURE_MATRIX_TSV`, and
`GRIP_HMP_OUTPUT_DIR`.

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
