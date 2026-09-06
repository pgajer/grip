# Data-to-geodesic embedding research workspace

Status: **planned/deferred**, separate from the focused methods paper in
[`../gmds_manuscript/`](../gmds_manuscript/). This is a scientific planning and
literature workspace; no principal `simods_data_geodesic_mds.tex` exists here.
A targeted September 6 location audit found no successor master in the related
geodesic-data-geometry, geosmooth, legacy geodesic-MDS, grip or shared-paper trees.
The historical SIMODS target is a proposal, not a submission commitment.

The broader agenda combines graph construction and intrinsic-distance estimation
from `geodesic_data_geometry`, embedding methods from grip, and non-oracle graph
selection. Its scope is not merged into the focused methods manuscript.

- `notes/simods_data_geodesic_mds_paper_scope.md`: historical scientific scope.
- `notes/`: scientific benchmark and graph-selection plans, plus generated dashboards.
- `evidence/`: claim, resource and figure/table registries.
- `literature/` and `paper_reviews_latex/`: scholarly reviews, source evidence and figures.
- `scripts/`: dashboard builders using local inputs.
- `provenance/migration.json`: original paths and hashes.

Regenerate the four dashboards with `make dashboards` (Python and PyYAML).
Their contents describe the deferred research agenda; linked legacy experiments
remain in geodesicMDS for reproducibility. Generated scratch crops and compilation
intermediates were omitted from migration; original files are backed up.

This directory is the canonical workspace for that agenda. Future writing should
establish one root manuscript source here rather than creating drafts in `build/`.
