# Repository Instructions

## Scope
- This repository is the source of truth for the `grip` R package, its pkgdown site inputs, vignettes, tests, and package-facing papers/tools that live in-repo.

## Preferred Skills
- Prefer `$r-package-qa` for package QA, release-readiness work, README and vignette drift, and CRAN-style checks.
- Prefer `$manuscript-build-review` for the in-repo paper and manuscript areas under `dev/`, `papers/`, and related report-rendering work.

## Package Source
- Core package code lives in `R/`, `src/`, `man/`, `tests/testthat/`, `vignettes/`, and `inst/`.
- Research notes, design work, and paper assets live under `dev/`, `tools/`, `papers/`, and `output/`.
- Treat `docs/` and `doc/` as generated or publication-facing artifacts unless the task explicitly targets site output.

## R Package Hygiene
- Edit roxygen comments in source files, not generated `.Rd` files.
- When changing exported functions, documentation, or examples, regenerate package docs with `make document`.
- When changing the README, edit `README.Rmd` and regenerate derived files with `make readme`.
- Prefer Makefile targets over ad hoc package commands:
  - `make check-fast` for a quick CRAN-style pass
  - `make check` for a fuller validation pass
  - `make paper-pdf` or related report targets for in-repo paper rendering

## Generated Files
- Do not hand-edit generated files such as `NAMESPACE`, `man/*.Rd`, `README.md`, or pkgdown output unless the task is specifically about generated artifacts.
- Keep local build products such as tarballs, `.Rcheck/`, and transient logs out of functional commits unless the task is release-specific.
- Prefer writing experimental diagnostics and one-off report outputs under `output/` or `tmp/` rather than cluttering package source directories.

## Private Agent Work Products
- Store internal agent-only work products under `~/.codex/private/grip/`, not in the repository. This includes internal audits, agent-to-agent handoffs, intermediate rewrites, working copies of reviewer reports used for agent tasks, historical prompts, and generated review diffs that are not intended as package, manuscript, reproducibility, or submission artifacts.
- Organize private material first by workstream and then, when useful, by artifact type. Use `~/.codex/private/grip/rjournal-paper/` for R Journal paper work, with subdirectories such as `audits/`, `handoffs/`, `drafts/`, and `diffs/`. Use parallel workstream directories such as `cran-release/`, `package-qa/`, or another clearly named scope for other internal work rather than mixing unrelated material into the R Journal directory.
- Maintain a short `README.md` in each workstream directory that identifies the files, their origin or former repository location, purpose, whether they were previously tracked, whether any item is expected to become a formal deliverable, and confirmation that repository builds, tests, manuscript renders, and validation workflows do not require them.
- Keep formal and publication-facing assets in the repository. In particular, do not move manuscript source, bibliography, figures, rendering tools, citation-verification evidence, reproducibility data or scripts, checksums, provenance records, package source, documentation, or tests into the private tree.
- Treat draft responses, internal referee simulations, and agent working copies of received reports as private. If a response-to-reviewers document becomes part of an actual submission, copy its finalized submission version into the appropriate repository submission bundle while retaining working drafts in the private tree.
- Do not make repository builds depend on files under `~/.codex/private/grip/`; the private tree is working context, not a reproducibility input. When retiring a tracked internal file from the repository, preserve its Git history through the repository deletion and record its private destination in the workstream README.
- The private directory is not a credentials store. Do not place passwords, access tokens, private keys, or other authentication secrets there.

## Review and Release Safety
- Preserve public API names and documented behavior unless the user explicitly requests a breaking change.
- For experimental methods, keep the user-facing positioning consistent with the package docs and README.
- When preparing release-oriented changes, verify documentation, tests, and pkgdown-facing material stay aligned.
