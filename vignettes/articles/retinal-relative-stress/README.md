# Relative-distance error in retinal metric MDS

[Article source](main.tex) | [Local PDF](main.pdf)

This article explains why regional sparse SGD, its endpoint-average proxy, and
equally weighted all-pairs relative stress are different. It contains a
reproducible pilot on 1,500 cells from the retinal PC20 matrix: regional methods
with 25 and 200 landmarks versus all-pairs SGD, three matched random starts,
and 500 epochs. Evaluation uses every unordered pair and terminal coordinates,
without post-fit rescaling.

The all-pairs reference uses the explicit-constraint interface with every pair
and unit endpoint multiplicities. It is mathematically an all-pairs
inverse-squared objective even though it uses the sparse-capable SGD kernel.
This matches the existing small-reference approach in
`inst/scripts/sparse-mds-methods-validation.R`.

This is a standalone methodological article, excluded from the CRAN package
archive by the existing `vignettes/articles` build exclusion. The retinal data,
cell IDs, large saved objects, machine-specific provenance, and internal
run records remain outside the public package tree. Compact numeric results
and reproducible analysis scripts accompany the article.

## Reproduce

Use optimized local grip and dgraphs source builds. The input is an RDS list
with `scores` (120,804 by 20) and `cells`, specifically the PC20 input already
used for the retinal demonstrations. No raw expression preprocessing is
performed here. Exact input and source hashes are recorded by the runner.

From the grip root, with `PC20_RDS` pointing to the input and `RESULT_DIR`
pointing to a new private directory:

```sh
env -u R_HOME -u R_LIBS -u R_LIBS_USER -u R_LIBS_SITE \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  Rscript vignettes/articles/retinal-relative-stress/compare.R \
  "$PWD" "$PC20_RDS" "$RESULT_DIR"
Rscript vignettes/articles/retinal-relative-stress/validate.R "$RESULT_DIR"
Rscript vignettes/articles/retinal-relative-stress/summarize.R \
  "$RESULT_DIR" vignettes/articles/retinal-relative-stress
Rscript vignettes/articles/retinal-relative-stress/build.R \
  vignettes/articles/retinal-relative-stress
```

The runner refuses to overwrite an existing protocol. It stores common starts,
subset indices, preparations, checkpoint coordinates, settings, source hashes,
and fit-status accounting. It verifies source stability before returning.
No cache is silently reused. The companion validator independently recomputes
all checkpoint errors. The summarizer and PDF builder do not rerun optimization.

## Vignette Development

The article can become a vignette after separating its methodological
explanation from this private-data reproduction. A small, publicly distributable
example should carry the runnable tutorial. The same contract tests can remain:
all-pairs proxy identity, terminal-coordinate policy, identical common starts,
and exact continuation equivalence.

The current pilot varies landmark count and initialization only. It is not an
ablation of endpoint asymmetry, schedule choice, or preparation randomness,
and is not a CRAN-readiness check.
