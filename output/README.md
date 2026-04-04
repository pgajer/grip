# Output Layout

Generated artifacts belong under paper-specific output directories.

Preferred layout:

- `output/rjournal_paper/`
- `output/gkk_lgkk_paper/`
- `output/geodesic_mds_paper/`
- `output/shared/`
- `output/diagnostics/`

Each paper output root can contain:

- `benchmarks/` for benchmark data, figures, tables, and manifests
- `reports/` for rendered PDFs and their sidecars
- `html/` for interactive galleries and widget bundles
- `tmp/` for scratch output that is still worth keeping under the paper root

Migration rule:

- Do not recreate `output/pdf/`.
- Do not recreate the old category-first buckets such as `output/benchmarks/` or `output/reports/`.
- New outputs should be grouped first by paper, then by artifact type and experiment.
- Diagnostics that are not owned by a paper belong under `output/diagnostics/`.
