# MDS initialization and neighborhood sensitivity: results

Across all 25 evaluated graphs, edge-KK reduces fixed-path error in 25/25 classical-initialized and 25/25 stress-initialized comparisons. It increases chord profiled Stress-1 in 25/25 and 25/25 comparisons, respectively. Stress initialization gives the lower final path error on 23/25 graphs (14/15 at the selected and nearby choices). This is an empirical preference, not uniform superiority over graph choices or global minimizers.

At the selected graphs, median path error is 3.451% for classical MDS and 2.264% for stress MDS; after edge-KK it is 0.245% and 0.164%, respectively. The corresponding median chord scores increase from 3.459% to 4.715% and from 2.803% to 5.213%.

Geometric recovery is less consistent. The stress-initialized edge-KK coordinate error at selected graphs ranges from 4.60–27.54% after similarity alignment. Good fixed-path agreement therefore does not establish correspondence or surface recovery. Across the selected and nearby choices, its largest coordinate error is 28.11%.

[Full PDF report](../report.pdf) · [Protocol and reproduction](../README.md)

## Scope

Five saved area-uniform clouds, n=1,000, on z=0.8(x²−y²) over [-1,1]²; 25 saved
symmetric-union kNN graphs; four fitted 3D candidates and the original-coordinate
control. MDS input is the strict graph shortest-path matrix built from ambient
chord-weighted edges. Numerical smooth-surface geodesics are evaluation references,
not the MDS inputs in this phase. Increasing-radius/geodesic-input experiments
remain a separate phase.

## Selected graphs

All entries below are medians across five paired clouds. Edge and path scores
use separately fitted target scales; chord is profiled Stress-1. Coordinate and
surface scores use similarity alignment. Surface RMS is in coordinate units.

| Method | Path (%) | Edge (%) | Chord Stress-1 (%) | X→Z path (%) | Coordinate (%) | Surface RMS |
| --- | --- | --- | --- | --- | --- | --- |
| Classical MDS | 3.451 | 7.133 | 3.459 | 4.553 | 29.190 | 0.240 |
| Stress MDS | 2.264 | 5.081 | 2.803 | 3.020 | 23.204 | 0.195 |
| Classical + edge-KK | 0.245 | 0.382 | 4.715 | 0.253 | 30.671 | 0.215 |
| Stress + edge-KK | 0.164 | 0.227 | 5.213 | 0.224 | 22.324 | 0.154 |

![Selected comparison](../figures/selected-comparison.png)

Exact per-cloud scores, in percent where appropriate: [selected-scores-percent.csv](selected-scores-percent.csv).
All returned-scale raw stresses, raw target-normalized RMSE, literal Stress-1,
rigid coordinate error, and singular-value ratios: [scores.csv](scores.csv).

## Neighborhood sensitivity

The selected-graph X→G errors are 0.203–0.217%. They are 0.697–0.814% at k=32 and 0.220–0.243% at k=80. The original-coordinate control has essentially zero graph-to-layout edge/path error at every k, while retaining the graph/reference error. This demonstrates why graph-to-layout fidelity cannot substitute for X→G validation.

The final-path advantage of stress over classical initialization reverses at:
**cloud 5, k=78; cloud 5, k=80**. These are retained results. Each branch still has its own
initializer-to-refinement comparison in the full tables. The experiment does
not establish monotonic improvement with k or a universally optimal initializer.

All 25 evaluated graphs have connectivity and degree/length statistics
in [graph-statistics.csv](graph-statistics.csv). Components before repair range
from 1 to 1; the
total number of added bridges over these graphs is 0.
The full calibration plateau and reference-sensitivity fields are retained in
[graph-selection.csv](graph-selection.csv). The five k values per cloud do not
cover every integer in the plausible region.

![Distance sensitivity](../figures/k-distance-sensitivity.png)

## Geometric agreement

![Geometric sensitivity](../figures/k-geometry-sensitivity.png)

Coordinate correspondence and closest-surface agreement answer different
questions from edge/path fidelity. Surface scoring uses the original parameter
triangulation, a twice-subdivided lifted reference, and 8,000 samples per direction.
Original-mesh surface RMS is 0.0073–0.0471;
this discretization control is not subtracted from candidate scores. Both rigid
and similarity alignment and surface Monte Carlo SE are recorded in
[surface-scores.csv](surface-scores.csv). Surface RMS is not a fold, topology,
or injectivity certificate.

![First-cloud illustration](../figures/selected-shapes.png)

The visual example uses cloud 1 by a fixed first-cloud rule. All panels have equal
units and the same limits after similarity alignment. Color identifies the
original height of corresponding vertices. No anisotropic display scaling is used.

## Optimization and additional budgets

At selected graphs, the worst of six achieved raw stresses exceeds the best by 2.37–9.72%. Across the full sweep, 28/150 primary starts reach the iteration limit; 0/25 selected starts do so. The other primary starts stop at the backend stress-change tolerance. The selected-graph 2,000-step tighter-tolerance continuation lowers raw stress by 0.040–0.397%. All 260 recorded edge-KK stages reach their iteration caps. The achieved finite-budget fits must not be called global optima.

![Optimizer sensitivity](../figures/optimizer-sensitivity.png)

The additional MDS result is a continuation of the selected fit, not a new
multi-start search, and does not replace the primary initializer. Each extra
edge-KK fit extends its own primary fit. Small raw-stress changes can accompany
noticeable coordinate changes; an objective tolerance is not a recovery guarantee.
See [starts.csv](starts.csv), [additional-budget-scores.csv](additional-budget-scores.csv),
[additional-mds-starts.csv](additional-mds-starts.csv), and
[edge-optimizer-status.csv](edge-optimizer-status.csv).
Elapsed times in [timings.csv](timings.csv) are concurrent-worker component
measurements, not a controlled performance benchmark.

## Numerical validation and input limitation

The complete verifier passed: 25 graphs, 125
primary score rows, 150 starts, 280
surface rows, and 15 additional fits. Independent
recomputation, control realizations, raw/profile scale identities, least-stress
selection, descent within edge-KK stages, and input/result checksums were checked.
The maximum historical score discrepancy is 8.26e-11;
all original pilot inputs remain unchanged. See [validation.json](validation.json).

An initial partial run exposed near-tie asymmetry in grip's retained-path
preparation (2.67e-8 on one graph), rejected by the metric-MDS symmetry validator.
Every final fit was rerun with the saved symmetric strict graph distances under
protocol v2. Retained routes still define path diagnostics; strict distances
define the chord objective. The partial run is excluded. The package validator
was not weakened and no package algorithm was modified. A general package repair
for path-prepared near-tie inputs remains a separate integration issue; the
current experiment supplies its verified symmetric inputs explicitly.

Reference-geodesic numerical validation is inherited from the original pilot.
Its sampled discrepancies are not certified all-pairs bounds or statistical
confidence intervals. These five samples and finite optimization budgets do not
establish population convergence or uniqueness of an embedding.

## Recommendation for the R Journal manuscript

Include a concise initializer comparison and plausible-k sensitivity result in
the main paper. The supported message is that edge-KK improves retained-path
fidelity under both MDS initializers, with a chord-fidelity tradeoff; evaluating
the graph against a reference remains a separate requirement. Describe the
stress-initializer advantage as frequent in these samples, not uniform. Avoid
claiming that low path error recovers the generating saddle or its pointwise
geometry. Put the full k grid, starts, graph statistics, surface/coordinate
checks, and budget sensitivity in supplementary reproducibility material.

The study uses the existing ratio-SMACOF implementation described by
[Mair, Groenen, and de Leeuw (2022)](https://doi.org/10.18637/jss.v102.i10);
its substantive citation evidence is in [citation_verification.html](../citation_verification.html).
The algebraic scale identity and numerical results here are independently checked.

Phase 2 does not edit or replace the main manuscript, old fits, or public package
API. Radius studies and manuscript integration remain later phases.
