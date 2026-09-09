# grip Design Documents

Public technical design notes, benchmark specifications, and prototypes,
organized by topic. Internal execution plans and task handoffs are not repository
inputs and are maintained separately.

## Topic directories

- `grip/` — GRIP core, globalrep, cross-family experiments, insertion/refinement
- `weighted-grip/` — weighted GRIP technical notes (MISF, caches, insertion, refinement)
- `graph-families/` — synthetic graph family specs, generators, geometry gallery
- `gripui/` — Shiny app design specifications

## Output directories

- `figures/` — exploratory figures from design phase
- `interactive-prototypes/` — HTML/JS interactive visualizations
- `pdf/` — benchmark PDF output directories (by experiment name)
- `tmp/` — transient benchmark run artifacts

## Notes

- Keep package documentation generators in `tools/pkg/`.
- Manuscript studies, their scripts, and their outputs live in separate private
  manuscript workspaces. Do not store them on public artifact branches.
- Keep disposable package previews in `output/`.
- When an asset is finalized for user-facing docs, move it to `man/figures/` and reference from README or vignette.
