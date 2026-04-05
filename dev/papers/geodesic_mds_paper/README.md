# Geodesic MDS Paper Workspace

This is the canonical workspace for the geodesic MDS paper inside the
`grip` repository.

Active source of truth:

- manuscript source: `dev/papers/geodesic_mds_paper/manuscript/geodesic_mds.tex`
- build script: `dev/papers/geodesic_mds_paper/scripts/build_geodesic_mds_pdf.sh`
- figures: `dev/papers/geodesic_mds_paper/figures/`
- notes: `dev/papers/geodesic_mds_paper/notes/`
- interactive companions: `dev/papers/geodesic_mds_paper/interactive/`

Versioning policy:

- keep exactly one active manuscript source in `manuscript/geodesic_mds.tex`
- use Git history in the `grip` repo as the canonical source history
- keep named manuscript snapshots only in `manuscript/archive/`
- keep build metadata stamped into the compiled PDF through the paper-local
  build script

Legacy locations:

- `dev/papers/gmds/` is a frozen legacy workspace and is not the source of truth
- `/Users/pgajer/current_projects/geodesic_MDS/` is a frozen external legacy
  workspace retained only for migration traceability
