# GRIP Software Paper

This directory is the canonical home of the GRIP software-package manuscript.
The source of truth is `grip-software-paper.Rmd`. Rendered PDF and HTML
versions, knitted and LaTeX intermediates, wrapper files, and figure sidecars
are generated under `build/` unless explicitly tracked. The
`citation_verification.html` page remains at this level as tracked audit
evidence rather than a manuscript render.

From the package repository root, render and audit the paper with:

```bash
make paper-all
make paper-citation-check
make paper-submission-bundle
```

Repository-level rendered deliverables are written under
`output/rjournal_paper/` by the reporting workflow.

The submission target renders stable PDF and HTML files without the internal
draft build stamp, verifies citations, and writes a clean, versioned directory
and ZIP archive under `output/rjournal_paper/submission/`. The archive includes
the manuscript sources, figure sidecars, motivation letter, package list, and
the self-contained reproduction supplement.

The manuscript currently targets grip 0.2.0. Do not submit the archive until
that exact version is publicly visible on CRAN and the final package and paper
checks have been rerun against it.
