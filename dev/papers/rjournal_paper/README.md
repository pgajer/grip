# R Journal Support Workspace

The canonical GRIP package-paper manuscript now lives at:

```text
papers/grip-software-paper/grip-software-paper.Rmd
```

This directory retains auxiliary and historical R Journal material:

- `figures/`: figure-generation source only
- `tables/`: paper tables and curated CSV inputs
- `notes/`: internal writing notes and revision plans
- `submission/`: journal-specific packaging and cover-letter materials

Historical scaffolding also exists under `dev/papers/r-journal/` and
`manuscript/legacy_r_journal_drafts/`; neither is an active manuscript source.

Build and citation-audit commands:

```bash
make paper-all
make paper-citation-check
```

Generated figures belong under `output/rjournal_paper/figures/`; rendered
manuscripts and sidecars belong under `output/rjournal_paper/reports/`.
