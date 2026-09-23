# Display registrations for the vignettes

These five RDS files contain surface-shape display registrations derived from
the existing metric-MDS, weighting, and sparse-MDS experiment bundles. They do
not replace or modify fitted coordinates, stress, timings, seeds, or outcomes.

`make documentation-assets` rebuilds them explicitly. Ordinary vignette builds
verify the reference, triangles, fitted display coordinates, and alignment
algorithm against the saved fingerprint, then reuse the registration. A stale
record stops the build rather than silently showing a different alignment.

The full-MDS records include the denser surface audits. Sparse records use the
independently fitted scale multipliers already reported in that experiment,
then translation, rotation and reflection only. The large saddle retains all
4,096 vertices and 7,938 triangles, with one area-weighted sample per triangle
for registration. See the complete experiment vignettes for the algorithm.

The source experiment bundles remain under metric-mds-comparison,
metric-mds-weighting, and sparse-mds-comparison. Each document render writes
build-manifest.tsv with checksums of sources and all saved RDS inputs.
