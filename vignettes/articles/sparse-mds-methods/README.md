# Sparse metric MDS: methods and synthetic validation

[LaTeX source](methods.tex) · [Scores](results/fits.csv) · [Provenance](results/provenance.json)

This note describes uniform landmarks, Euclidean regional sparse stress and
graph-geodesic regional sparse stress as implemented on October 7, 2026 in
`grip` and `dgraphs`. It covers internally selected pivots, regional weights and
component-MST repair. Repair streams distances without retaining a full distance
matrix, but its worst-case distance work remains quadratic. Dense calculations
occur only in the deliberately small references.

The experiment completed 288/288 fits on four 120-point quadratic forms:
latent dimensions 2 and 3, index 0 and 1, outputs 2 and 3, three seeds,
16 pivots, q=5/10/20 and 100 SGD epochs with the published weight-based schedule.
The note distinguishes Euclidean from geodesic targets. These results establish
neither convergence nor large-data performance or V4 parameter adequacy.

A [retinal relative-stress article](../retinal-relative-stress/README.md)
examines how the asymmetric regional update rule and its reported proxy differ
from all-pairs relative stress. Its small-data comparison uses common starting
coordinates, exact all-pairs evaluation and 500-epoch SGD trajectories.

## Saved scientific results

- [All 288 scores](results/fits.csv), [graph diagnostics](results/graphs.csv),
  [PCA baselines](results/pca-baselines.csv), [paired comparisons](results/paired_comparisons.csv),
  [summary statistics](results/summary.csv), and [neighbor agreement](results/neighbor-agreement.csv).
- [Figures 1–4: primary comparisons](results/error-comparison.pdf) and
  [Figures 5–8: paired neighborhood sensitivity](results/neighborhood-sensitivity.pdf).
  Seed ranges are descriptive, not confidence intervals.
- [Benchmark-era source fingerprints](results/source_manifest.csv) retain the
  original MD5 identities with paths relative to a common parent of the package
  checkouts. [Result provenance](results/provenance.json) records the generation
  time and SHA-256 checksums of the bundled tables and figures.

At that benchmark revision, 6,665 grip checks and 112 focused dgraphs checks
passed without failures, warnings or skips. The final run included corrected
handling of coincident PCA scores for distinct original observations. Source
and compiled-library hashes were stable throughout the run. These historical
checks do not establish the status of subsequent package revisions.

The final R process took 41.12 seconds of serial wall time and 443,318,272 bytes
(422.8 MiB) of maximum resident memory. This includes loading, preparation,
small dense references, fitting, scoring and result saving; it excludes native
compilation and rendering. It is not a sparse-only memory benchmark.

## Render the frozen results

From this directory, with `NOTE_OUTPUT` set to a directory outside the source tree:

```sh
python3 build_report.py --pdf --output-dir "$NOTE_OUTPUT"
```

The builder checks the saved result checksums, regenerates the numerical tables,
and compiles the note twice. It copies the supplements beside the exported PDF
and writes build records to the requested output directory. Rendering does not
run fits or test the current implementation against the historical results.
Python 3.9 or later and a LaTeX installation with the packages named in
`methods.tex` are required. Omit `--pdf` to refresh only the source and supporting
files. A citation check can be run with an independently maintained evidence file:

```sh
make citation-check CITATION_EVIDENCE="$CITATION_EVIDENCE" BUILD_DIR="$NOTE_OUTPUT"
```

## Run a new scientific experiment

Use the `grip` and `dgraphs` source checkouts. Compile optimized native libraries
with `pkgbuild::compile_dll(debug=FALSE, force=TRUE)`, then load them with
`pkgload::load_all(..., recompile=FALSE)`. From the grip root, with `RESULT_DIR`
set to a new output directory:

```sh
Rscript --vanilla inst/scripts/sparse-mds-methods-validation.R "$PWD" "$RESULT_DIR"
Rscript --vanilla tools/pkg/summarize-sparse-methods.R "$RESULT_DIR" "$RESULT_DIR"
```

The runner refuses to overwrite an existing completed score table. Data seeds
are `20261020 + 10*d + index`, preparation seeds `20261030 + 1:3`, and fitting
seeds `1:3`; initialization is matched across methods within each case,
dimension and seed. Settings, input coordinates, preparations, layouts, source
identities and session information are saved by the runner. Compare a new run
with the frozen results separately: current source versions can differ from
the benchmark-era implementation, and the frozen tables must not be silently
replaced with new measurements.
