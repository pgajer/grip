# H005 Polished Paper Reviews

This directory contains near-public, research-paper-style reviews derived from
the H005 literature-review workflow.

Source hierarchy:

1. The canonical local PDF is the primary source for each review.
2. The Markdown memo and audit memo are archival evidence and structure aids.
3. The final H005 synthesis supplies cross-paper context, not paper-specific
   authority.
4. `../paper_figure_screenshots.yml` is the provenance source for reproduced
   paper-figure screenshots.

Generated HTML/PDF outputs are derived files. Edit each paper's `P##_review.tex`
source and rebuild rather than editing rendered output.

## Writer Expectations

Each writer reads the full paper PDF, including figures and tables, before
drafting. The review should read as an intelligible paper review for a
mathematically capable non-specialist. It should not simply paraphrase the
Markdown memo.

## Auditor Expectations

Each auditor independently reads the full paper PDF, the writer's LaTeX review,
the archival Markdown memo, the original audit memo, and the figure manifest.
The auditor checks scientific accuracy, clarity, self-contained background,
natural flow, figure interpretation, limitations, and SIMODS/gflow relevance.

## Build

Build all available reviews:

```bash
python3 build_paper_reviews.py
```

Build a single review:

```bash
python3 build_paper_reviews.py P03
```

