# grip Design Documents

Public technical design notes, benchmark specifications, and prototypes,
organized by topic. Internal execution plans and task handoffs are not repository
inputs and are maintained separately.

## Topic directories

- `grip/` — GRIP core, globalrep, cross-family experiments, insertion/refinement
- `weighted-grip/` — weighted GRIP technical notes (MISF, caches, insertion, refinement)
- `gripui/` — Shiny app design specifications

## Output directories

- `figures/` — exploratory figures from design phase
- `interactive-prototypes/` — HTML/JS interactive visualizations
- `pdf/` — benchmark PDF output directories (by experiment name)
- `tmp/` — transient benchmark run artifacts

## Notes

- Historical graph-family catalogs, development galleries and benchmark-selection
  notebooks are maintained in the separate manuscript workspace. Current public
  examples are in `vignettes/articles/gallery.Rmd`; function contracts are in
  package help. Public builds must not read private historical documents.

- Keep package documentation generators in `tools/pkg/`.
- Manuscript studies, their scripts, and their outputs live in separate private
  manuscript workspaces. Do not store them on public artifact branches.
- Keep disposable package previews in `output/`.
- When an asset is finalized for user-facing docs, move it to `man/figures/` and reference from README or vignette.
