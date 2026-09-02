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
make paper-supplement
make paper-citation-check
make paper-submission-bundle
```

Repository-level rendered deliverables are written under
`output/rjournal_paper/` by the reporting workflow.

The submission target renders stable PDF and HTML files without the internal
draft build stamp, verifies citations, and writes a clean, versioned directory
and ZIP archive under `output/rjournal_paper/submission/`. The archive includes
the manuscript sources, figure sidecars, motivation letter, package list, and
the self-contained reproduction materials, Supplement S1 on weighted-GRIP
complexity, and Supplement S2 with the full higher-dimensional Möbius-strip
comparison. Both supplements' canonical LaTeX sources, bibliographies, and
citation evidence live in `supplement/`; PDFs and intermediates are generated
under `build/supplement/`. The regular paper-rendering targets also build and
check both supplements. Their PDFs, sources, and portable build instructions
are included in the submission archive.

The Möbius comparison uses shared functions in
`reproducibility/scripts/mobius-dimension-comparison.R`: the manuscript
regenerates the 300-vertex pair, and S2 shows 150-, 300-, and 1,500-vertex pairs.
The archive includes the full figure and its coordinate/provenance RDS.
This is a modern, single-seed illustration, not a replication of the historical
GRIP experiment or a demonstration of improved graph-distance fidelity.

The complexity bounds in S1 are conditional on its explicit locality and
hierarchy assumptions. They are not empirical scalability results or a
guarantee for arbitrary weighted graphs. The motivation letter requests
mathematical review of the analysis; no independent review is claimed here.

The manuscript currently targets grip 0.2.0. Do not submit the archive until
that exact version is publicly visible on CRAN and the final package and paper
checks have been rerun against it.
