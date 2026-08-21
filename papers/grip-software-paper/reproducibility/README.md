# Reproducing the grip Software-Paper Results

This directory makes the two precomputed components of the paper inspectable
without relying on paths elsewhere in the development repository. The paper
itself evaluates the smaller examples directly.

## Requirements

- R and grip 0.2.0 or later
- igraph and graphlayouts for the cross-package benchmark
- the packages listed in the manuscript-level `_Rpackages.txt`

## Included artifacts

- `precomputed/vs_alternatives/benchmark_results.rds`: cross-package layouts,
  scores, timings, and generation-session information.
- `precomputed/hmp_u01_gc_coarse/vignette_results.rds`: HMP candidate layouts,
  scores, and timings used in the paper.
- `hmp_u01_gc_coarse/`: the exact weighted edge list, vertex labels,
  original-to-coarse membership, coarsening history, provenance, and checksums
  for the HMP graph used in the paper.
- `SHA256SUMS`: checksums for every supplied binary or tabular input and
  precomputed result.

## Cross-package benchmark

From this directory, run:

```sh
GRIP_VS_ALTERNATIVES_OUTPUT=generated/benchmark_results.rds \
  Rscript scripts/precompute-vs-alternatives.R
```

The script uses fixed layout and scoring seeds. Elapsed times remain
machine-dependent.

## HMP graph and candidate search

First rebuild the graph object from the supplied tables:

```sh
Rscript scripts/build-hmp-data-object.R
```

Then recreate the candidate-search artifact:

```sh
GRIP_HMP_DATA_RDS=generated/hmp_u01_gc_coarse.rds \
GRIP_HMP_RESULTS_OUTPUT=generated/hmp_vignette_results.rds \
  Rscript scripts/precompute-hmp-results.R
```

The candidate search uses fixed seeds. Timings depend on the machine. The
provenance record states the boundary of reproducibility: the supplied tables
reproduce the exact coarsened graph used by the paper, but they do not rerun the
earlier amplicon-processing and iKNN graph-selection workflow.

## Fast article build

The manuscript reads the supplied files under `precomputed/` so that reviewers
do not need to rerun either heavier workflow. The scripts above document and
execute the complete derivation from the supplied graph inputs to the saved
paper artifacts.
