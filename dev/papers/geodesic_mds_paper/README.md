# Geodesic MDS Paper Workspace

This is the canonical workspace for the geodesic MDS paper inside the
`grip` repository.

Active source of truth:

- manuscript source: `dev/papers/geodesic_mds_paper/manuscript/geodesic_mds.tex`
- build script: `dev/papers/geodesic_mds_paper/scripts/build_geodesic_mds_pdf.sh`
- rendered PDF: `dev/papers/geodesic_mds_paper/build/geodesic_mds.pdf`
- figures: `dev/papers/geodesic_mds_paper/figures/`
- notes: `dev/papers/geodesic_mds_paper/notes/`
- interactive companions: `dev/papers/geodesic_mds_paper/interactive/`

Versioning policy:

- keep exactly one active manuscript source in `manuscript/geodesic_mds.tex`
- use Git history in the `grip` repo as the canonical source history
- keep named manuscript snapshots only in `manuscript/archive/`
- keep build metadata stamped into the compiled PDF through the paper-local
  build script
- treat `manuscript/geodesic_mds.pdf` as a local render artifact rather than a
  tracked source file

Active manuscript provenance:

- the active manuscript source is `dev/papers/geodesic_mds_paper/manuscript/geodesic_mds.tex`
- this file was created during the paper-workspace migration in commit
  `b750762`
- its initial content was promoted from the legacy `geodesic_mds_v3.tex`
  manuscript, with path normalization for the canonical in-repo paper layout
- the archived snapshot remains at
  `dev/papers/geodesic_mds_paper/manuscript/archive/geodesic_mds_v3.tex`

Important distinction:

- `archive/geodesic_mds_v3.tex` is a historical snapshot
- `manuscript/geodesic_mds.tex` is the only active manuscript source
- after migration, additional edits were made directly to
  `manuscript/geodesic_mds.tex`, so it is not expected to remain identical to
  archived `v3`

Policy:

- make all current manuscript edits in `manuscript/geodesic_mds.tex`
- preserve numbered `vN` files only as archive snapshots
- if a new archival snapshot is needed, create it explicitly from the current
  active manuscript

Legacy locations:

- `dev/papers/gmds/` is a frozen legacy workspace and is not the source of truth
- `/Users/pgajer/current_projects/geodesic_MDS/` is a frozen external legacy
  workspace retained only for migration traceability
