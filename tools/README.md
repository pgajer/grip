# Tools Layout

Benchmark and reporting scripts should be grouped by paper or shared use.

Preferred layout:

- `tools/benchmarks/gkk_lgkk_paper/`
- `tools/benchmarks/geodesic_mds_paper/`
- `tools/benchmarks/shared/`
- `tools/reports/gkk_lgkk_paper/`
- `tools/reports/geodesic_mds_paper/`
- `tools/reports/rjournal_paper/`
- `tools/figures/shared/`
- `tools/pkg/`
- `tools/utils/`

Migration rule:

- New scripts should not be added at the top level of `tools/` unless they are temporary and about to be rehomed.
