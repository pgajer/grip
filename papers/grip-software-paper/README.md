# GRIP Software Paper

This directory is the canonical home of the GRIP software-package manuscript.
The source of truth is `grip-software-paper.Rmd`; the `.tex`, PDF, HTML, and
figure sidecars are generated build products unless explicitly tracked.

From the package repository root, render and audit the paper with:

```bash
make paper-all
make paper-citation-check
```

Repository-level rendered deliverables are written under
`output/rjournal_paper/` by the reporting workflow.
