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

## Output-Dimension Experiment

The companion `dimension-sweep.py` and `dimension-sweep.R` scripts test output
dimensions 300, 200, 100, 50, 25, 20, 10, 5, 4, and 3 on the same 1,500 cells.
This is a separate experiment: its targets are Euclidean distances in the
saved 1,814-PC matrix, not the earlier 20-PC matrix. The full-population PCA
retains 99.002% of variance in the centered, unscaled transformed gene matrix.

Each dimension receives one seeded regional sparse SGD fit with 200 internally
selected distance-proportional landmarks, ten local neighbors, and a shared
preparation. Neighbor search uses a subset PCA retaining at least 90% of
variance; all target distances still use the complete 1,814-PC input.
The Zheng schedule spans exactly 100 epochs with epsilon 0.1. Initial normal
coordinates share nested columns across dimensions and are centered and
scale-fitted before optimization; the fitted coordinates are not rescaled.
All 1,124,250 pairs are evaluated, separately from the sparse fitting proxy.
An unrescaled subset-PCA projection provides a reference, not an optimizer
of relative-distance error.

The runner verifies the existing projection hash and subset cell order, then
launches ten single-threaded workers. Its two-second watchdog limits aggregate
process-tree resident memory to 50 GiB; this is a sampled stop threshold, not
an operating-system allocation reservation. Individual fitting workspace
allowances are 4 GiB. Completed ten-epoch checkpoints are saved atomically.
It does not load or resume the paused full-population embedding.

```sh
python3 vignettes/articles/retinal-relative-stress/dimension-sweep.py \
  --root "$PWD" --input-dir "$PC1814_INPUT_DIR" \
  --prior-protocol "$PRIOR_RESULT_DIR/protocol.rds" --out "$DIMENSION_RESULT_DIR"
env -u R_HOME -u R_LIBS -u R_LIBS_USER -u R_LIBS_SITE \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  Rscript vignettes/articles/retinal-relative-stress/dimension-validate.R \
  "$DIMENSION_RESULT_DIR"
```

`PC1814_INPUT_DIR` contains `pc1814.npy`, `projection.json`, and
`annotations.csv` from the previously validated projection. The validator
independently recomputes final distance errors through Gram matrices and
writes a summary table and dimension-error figure. One start and one
preparation per dimension do not establish convergence, robustness to random
initialization, or the dataset's intrinsic dimension.

The [completed dimension sweep](dimension-summary.csv) gives all-pairs RMS
relative errors of 30.16% in 3D, 10.00% in 20D, 4.52% in 100D, 5.02% in 200D,
and 2.68% in 300D. These are terminal 100-epoch results from one start per
dimension, not best achievable errors. The increase from 100D to 200D is real
in these fits, even though the sparse proxy falls from 2.37% to 0.83%.
The subset PCA reference and its variance fractions appear in the same table;
those fractions describe this subset, not the full-population variance spectrum.

The article now includes the dimension results and their error plot. A
report-only diagnostic decomposes the 100D and 200D scores into directly
fitted pairs and the remaining pairs. The former improve from 2.545% to
0.896% RMS error; the latter worsen from 5.023% to 5.801%. Remaining pairs
account for 74.18% of all pairs, explaining the higher overall 200D error.
This diagnoses where the error increases, not why a particular optimizer
trajectory reaches that configuration.

Regenerate this material without rerunning fits:

```sh
env -u R_HOME -u R_LIBS -u R_LIBS_USER -u R_LIBS_SITE \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  Rscript vignettes/articles/retinal-relative-stress/dimension-report.R \
  "$DIMENSION_RESULT_DIR" vignettes/articles/retinal-relative-stress \
  "$PC1814_INPUT_DIR/annotations.csv"
Rscript vignettes/articles/retinal-relative-stress/build.R \
  vignettes/articles/retinal-relative-stress
```

The final section proposes a balanced developmental-stage comparison: 1,500
cells per stage, the same global 1,814-PC input coordinates, 200 landmarks,
100 epochs, and dimensions 3, 4, 5, 10, 20, 25, and 50. All ten stages have
enough cells. These 70 stage-specific fits are proposed, not completed.
