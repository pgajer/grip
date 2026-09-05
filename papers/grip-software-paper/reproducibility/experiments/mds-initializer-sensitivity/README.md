# MDS initialization and neighborhood sensitivity on the five saddle clouds

This experiment compares classical scaling and raw-stress metric MDS, with
and without edge-KK, on the existing manuscript samples. It preserves the
historical fits and does not replace manuscript figures. The radius study is
a separate phase. Read `summary/results.md` for the completed numerical results.

## Frozen design

The input is the five area-uniform samples of 1,000 points on the bounded saddle
`z = 0.8*(x^2-y^2)`, `(x,y) in [-1,1]^2`, from the
[original pilot](../two-fidelity-pilot/README.md). Its saved symmetric-union kNN
graphs use ambient chord edge lengths, exact neighbors, no pruning, and
component-MST repair if needed. **MDS receives graph shortest-path distances**;
the independently computed numerical smooth-surface geodesics are references
for assessing graph and end-to-end fidelity, not MDS inputs in this experiment.
The geodesic-input radius experiment remains separate.

The original selected graphs have k = 71, 70, 67, 72, 73. Before any new scientific
fits, `freeze.R` records the rule
`sort(unique(c(32, selected_k-5, selected_k, selected_k+5, 80)))`.
The resulting 25 graphs include both sides of each selected k and a broader
sparse/dense contrast. The five values per cloud are not a complete integer sweep.
The saved 1%-of-minimum calibration plateaus are 71–72, 64–71, 65–68, 70–73,
and 65–73. We record those intervals; we do not redefine them from embedding
outcomes. The nearby ±5 values probe a wider neighborhood than this narrow
plateau. Their calibration errors are within about 7% of each cloud's minimum.
The existing mesh-resolution/source-subset comparisons remain numerical
sensitivity checks, not confidence intervals or certified accuracy bounds.

All layouts are 3D. The candidates are original coordinates (an exact
edge/path control), classical MDS, stress MDS, classical MDS + edge-KK, and
stress MDS + edge-KK. The selected-graph comparison was reviewed before expanding to the other 20
graphs. The final consistency rerun repeats all 25 graphs with the strict input
matrices described below.

* Classical MDS calls `grip::classical.mds()`. At the selected graphs, we rerun
  the classical pipeline and compare its distances and scores with the original
  saved results. This also gives current, matched timing measurements.
* Stress MDS calls `grip::metric.mds()`, using the existing ratio-SMACOF backend
  and uniform all-pair raw distance stress. Each graph receives six independent
  calls: one classical start and five Gaussian starts, each with 1,000 iterations
  and normalized-stress tolerance `1e-8`. Seeds are
  `7300000 + 10000*replicate + start`; the same random configurations are used
  across k within a cloud. The least achieved raw stress selects the primary
  result. A finite local fit or tolerance stop is not a global certificate.
* Both edge-KK branches use the same density-to-uniform continuation
  `(0, .25, .5, .75, 1)`, 200 steps per stage, a profiled target scale, zero
  edge-length stabilizer, and the existing C++ Armijo solver. Input coordinate
  units are retained: SMACOF's normalization is undone by the package adapter.
* At each selected graph only, a separate budget check extends the selected
  stress-MDS fit for 2,000 iterations at tolerance `1e-10`, and separately
  extends each primary edge-KK fit for 1,000 uniform-stiffness iterations.
  These checks never replace primary results. The extended MDS fit is not used
  to initialize the primary edge-KK branch.
* Per-start checkpoints, stopping records, graph identities, versions, settings,
  and source/input hashes prevent old classical caches being mistaken for stress
  fits. Edge-KK traces retain scalar diagnostics but omit bulky coordinate frames.

## Score definitions

All graph diagnostics use the 499,500 unordered vertex pairs and one retained
input-graph shortest path per pair. Routes are fixed across candidates on the
same graph. They can change when k changes. The same 119,744 reference pairs
per cloud are used at every k. No embedded shortest paths are recomputed.

For observed lengths y and target lengths d, define `a = sum(y*d)/sum(d^2)`.
The reported edge and path relative errors are
`sqrt(sum((y-a*d)^2)/sum((a*d)^2))`, with separately fitted scales.
The chord **profiled Stress-1** is
`sqrt(sum((y-a*d)^2)/sum(y^2))`.
Raw chord stress is `sum((y-d)^2)` at the returned coordinate scale; its
target-normalized RMSE and literal unprofiled Stress-1 are also recorded.
Thus a chord score, a fixed-path score, and an optimizer objective have explicit,
different names. Raw-stress-optimal coordinates satisfy a coordinate multiplier
`sum(y*d)/sum(y^2)` close to one. With free scale, squared profiled Stress-1 equals
the minimum raw stress divided by `sum(d^2)`; the verifier checks this identity.

For X→G error, compare graph distances directly with the numerical smooth
reference. For X→Z path error, divide embedded path lengths by the edge-fitted
scale, then compare with that reference. No scale is fitted to surface distances.
We additionally record the G→Z path error at this common edge scale.

Geometric diagnostics use known row correspondence. Rigid alignment permits
translation and rotation/reflection; similarity alignment also fits a uniform
scale. Relative coordinate RMSE divides by the original cloud's RMS radius.
Surface RMS uses the original parameter-space Delaunay triangles for each
embedded mesh and a twice-subdivided, lifted reference mesh over the same sampled
footprint. It averages the two directional squared closest-triangle distances,
using 8,000 area-uniform samples per direction and seed `1901 + replicate`.
Rigid and similarity alignments are scored separately, without another surface
registration. Monte Carlo SE describes sampling only; small surface RMS does
not certify absence of folds or recovery of correspondence. The original mesh
provides a nonzero discretization control. Both smaller singular-value ratios
are recorded without anisotropic display rescaling.

## Reproduction and outputs

Run from the repository root with grip 0.2.0.9000, smacof 2.1-7, Rcpp, geometry,
and the original pilot caches in `papers/grip-software-paper/build/two-fidelity-pilot`.
Plot/report generation additionally uses Python with NumPy, pandas, and
Matplotlib, plus pdflatex and BibTeX. See the original pilot README to regenerate
those inputs. No build depends on
private notes. An optional `GRIP_MDS_SENSITIVITY_OUTPUT` changes only the output
root; use a new directory when changing scientific settings or fit sources.

```sh
Rscript papers/grip-software-paper/reproducibility/experiments/mds-initializer-sensitivity/freeze.R
# Repeat for r = 1,...,5; independent clouds may run concurrently.
Rscript papers/grip-software-paper/reproducibility/experiments/mds-initializer-sensitivity/run.R 1 selected
# Review the five selected-graph results before running the sensitivity stage.
Rscript papers/grip-software-paper/reproducibility/experiments/mds-initializer-sensitivity/run.R 1 sensitivity
Rscript papers/grip-software-paper/reproducibility/experiments/mds-initializer-sensitivity/surface.R
Rscript papers/grip-software-paper/reproducibility/experiments/mds-initializer-sensitivity/collect.R
Rscript papers/grip-software-paper/reproducibility/experiments/mds-initializer-sensitivity/verify.R
python3 papers/grip-software-paper/reproducibility/experiments/mds-initializer-sensitivity/plot.py
```

The build directory `papers/grip-software-paper/build/mds-initializer-sensitivity-strict/` stores full checkpoints and worker logs. The tracked
`summary/` directory contains compact numerical evidence, coordinates for figure
reproduction, protocol/environment records, checksums, validation results, and
the interpretation. `figures/` contains generated static scientific figures.
Plotting uses only the compact tracked summaries and never reruns optimization.

## Numerical preparation finding

The initial selected-graph stage passed, but cloud 3 at k=62 exposed a mismatch
between retained-path preparation and the new strict MDS symmetry validation:
the prepared matrix was asymmetric by `2.674647e-8`. The shortest-path-tree
near-tie rule in `grip.dijkstra.tree()` can retain a slightly longer route and
return tiny directional differences. This is much smaller than the scientific
fitting errors, but an asymmetric matrix is unsuitable as the MDS contract.

Every **final** classical and stress-MDS fit therefore uses the saved strict
`igraph` distance matrix, averaged with its transpose to remove floating-point
directional rounding. The retained flat routes are unchanged and continue to
define path scores. Strict graph distances define raw chord stress and Stress-1.
The original preparation's route targets and strict distances are checked to
agree within `1e-7`. All primary and additional fits were rerun in a new output
directory under protocol `saddle-initializers-k-v2`. The partial preliminary run
is excluded from the summaries; no symmetry check was weakened. The historical
classical-distance comparison uses `1e-7`, accommodating the observed `1.9e-9`
cloud-3 difference. The final score discrepancy is reported in validation.json.

This identifies a package integration limitation: `metric.mds(prepared=...)`
can reject a near-tied path-prepared matrix. The experiment supplies a verified
symmetric matrix explicitly. A package-level repair of general preparation is
separate from these scientific fits and must preserve retained-route semantics.
The package itself and the original pilot inputs were not modified in this phase.

Component timings are wall-clock measurements from concurrent cloud workers,
not a controlled performance benchmark. Warm-up, preparation-cache reuse at
selected graphs, and the six-start budget must be considered when comparing them.

`reproduce-near-tie.R` reduces the integration limitation to a three-vertex graph:
one direct edge has length `2 + 1e-8`, while the two-edge route has length `2`.
The path-prepared call is rejected; the distance-only metric-MDS call succeeds.
This diagnostic makes the pending package issue independently reproducible.

After collecting and verifying the complete experiment, generate the PDF readout
and rerun its citation gate with:

```sh
make -C papers/grip-software-paper/reproducibility/experiments/mds-initializer-sensitivity report
make -C papers/grip-software-paper/reproducibility/experiments/mds-initializer-sensitivity citation-check
```

`make artifact-check` in this directory verifies the source, summary, figure,
and PDF checksums written by the report build. The full R verifier additionally
requires the original input and final fit/path caches.
