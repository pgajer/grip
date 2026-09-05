# Expanding paraboloid and saddle MDS report

`report.tex` is the canonical comprehensive LaTeX document for the investigation.
It covers both smooth-geodesic experiments, all asymptotic arguments and
antipodal estimates, the ambient-distance comparison, numerical validation,
and the explanation of axis scaling. It includes all 15 distinct figures.

From the repository root:

```sh
make -C tools/reports/geodesic-mds pdf
make -C tools/reports/geodesic-mds verify
```

The build consumes existing verified results; it does not rerun MDS. If those
results are absent, first run the experiment targets listed in Appendix C of
the report. Python with NumPy and Matplotlib, `latexmk`, `pdflatex`, `bibtex`,
and `pdftotext` are required. The local experiment verifiers are invoked by
the report build. No private agent files or personal notes are build inputs.

Outputs are under `output/pdf/geodesic-mds/`:

- `report.pdf`: compiled report, with linked contents and references.
- `geodesic-mds-report.zip`: self-contained LaTeX source, compiled PDF, vector
  figures, generated tables, all CSV results, bibliography, citation evidence,
  and input/output checksums. Compile the extracted bundle with
  `latexmk -pdf report.tex`.
- `manifest.json`: checksums for consumed source/data/figure inputs and the
  deliverables, with an automatically generated Eastern build timestamp.

The figure atlas uses the plots' native sizes on larger pages so dense panel
labels do not become unreadable when fitted into a normal text column. The
axis-scaling illustration is regenerated as a vector plot with clearer ticks;
the other 14 figures are the existing vector exports without alteration.
`build.py` generates tables from the recorded CSVs and audits all citation
keys, caption counts, resolved cross-references, overfull boxes, and checksums.
The final PDF is also inspected visually after rendering.

The report distinguishes proved finite-sample statements, numerical local
fits, and population eigenoperator estimates. It does not claim a universal
nonplanarity theorem for raw-stress saddle MDS or a justified arbitrary joint
sample-size/large-radius limit for the population antipodal calculation.
