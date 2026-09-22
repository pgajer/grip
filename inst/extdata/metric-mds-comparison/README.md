# Synthetic metric-MDS comparison

This public reproducibility bundle belongs to the installed
`metric-mds-backends` vignette. It contains synthetic data only.

`benchmark.rds` records 14 datasets, public graph constructions and target
matrices, the full protocol and source fingerprints, classical-MDS references,
and every attempted SGD/SMACOF fit, including coordinates, timing, warnings and
failures. Five paired starts and independent allowances of 10, 30 and 100 give
420 planned fits. Displayed fits are always seed 1 and allowance 100.

The PNG files are images of the ivue overlays used as a fallback when the
optional ivue package is unavailable. Nonplanar examples are fitted in 3D;
the planar controls are fitted and drawn in 2D by the vignette itself.

Reproduction scripts:

- `inst/scripts/metric-mds-comparison.R`: public recipes; also installed.
- `tools/pkg/build-metric-mds-comparison.R`: bounded, isolated timing experiment.
- `tools/pkg/render-metric-mds-comparison.R`: ivue views of the saved fits.
- `tools/pkg/capture-metric-mds-comparison.cjs`: offline view checks and PNG capture.

The local graph-path targets are not continuous surface geodesics. The ambient
quadform targets are Euclidean distances represented by complete weighted graphs.
Figure alignment permits translation and orthogonal rotation/reflection, not scale.
Timing is platform-specific; these small examples do not establish large-data scaling.

The published measurements use a clean optimized (`-O2`) installation. The
previous unoptimized attempt remains a separate private diagnostic record and
is disclosed in the vignette and protocol; it is not mixed into these results.
Execution-time source fingerprints are retained. The protocol also stores the
exact fitting/fixture function source; subsequent script changes improve plots
and dependency preflight checks without changing those numerical functions.
