# grip Design Documents

Development-phase design notes, action plans, benchmarks, and prototypes organized by topic.

## Topic directories

- `grip/` — GRIP core, globalrep, cross-family experiments, insertion/refinement
- `weighted-grip/` — weighted GRIP phase plans (MISF, caches, insertion, refinement)
- `lgkk/` — landmark geodesic KK optimizer, bug fixes, test suites
- `gmds/` — geodesic MDS design, implementation plans, cleanup, pathology analysis
- `graph-families/` — synthetic graph family specs, generators, geometry gallery
- `gripui/` — Shiny app design and MVP implementation plans

## Output directories

- `figures/` — exploratory figures from design phase
- `interactive-prototypes/` — HTML/JS interactive visualizations
- `pdf/` — benchmark PDF output directories (by experiment name)
- `tmp/` — transient benchmark run artifacts

## Notes

- Keep reproducible generator scripts in `tools/`.
- Keep generated HTML, PDF, image, and temporary benchmark outputs on a dedicated artifacts branch such as `codex/gmds-artifacts`, not on the merge branch.
- Do not commit generated assets until we explicitly choose final figures.
- When an asset is finalized for user-facing docs, move it to `man/figures/` and reference from README or vignette.
