# Output Layout

Generated artifacts belong under paper-specific output directories.

Preferred layout:

- `output/benchmarks/gkk_lgkk_paper/`
- `output/benchmarks/geodesic_mds_paper/`
- `output/benchmarks/shared/`
- `output/reports/gkk_lgkk_paper/`
- `output/reports/geodesic_mds_paper/`
- `output/reports/diagnostics/`
- `output/galleries/`
- `output/html/`
- `output/tmp/`

Migration rule:

- Do not recreate `output/pdf/`.
- New outputs should be grouped by paper and experiment.
- Use `output/reports/...` for rendered PDFs and report sidecars.
