# Classical versus stress MDS: migration and experiment status

From grip 0.2.0.9000, the former `metric.mds()` implementation is named
`classical.mds()`. The new `metric.mds()` minimizes raw distance stress through
ratio SMACOF. Existing manuscript reproduction calls have been migrated to
`classical.mds()` to preserve their algorithm. The edge-KK default likewise
remains classical, with explicit name `init = "classical_mds"`.

Phase 1 preserved all saved fits. Phase 2 now supplies a separate
[initializer and neighborhood comparison](experiments/mds-initializer-sensitivity/summary/results.md);
it does not replace those historical results. Older labels such as `Metric MDS`, `metric-MDS`,
and serialized fields such as `metric_mds` refer to classical scaling in these
existing artifacts. The main manuscript and supplement will receive a complete
figure/label revision when the paired baseline comparison is incorporated.

| Existing workflow | Phase 1 action | Later comparison |
|---|---|---|
| Five-cloud two-fidelity saddle pilot | Preserve graphs, references, selected classical fits, and edge-KK results; migrate `fit-layouts.R` | Completed separately in Phase 2: five clouds, 25 graphs, both MDS initializers and their own edge-KK refinements |
| Main panel-E workflow | Call `classical.mds()` explicitly | Update panel inputs if the new baseline is selected |
| HMP precomputation | Preserve classical calls and saved field identities | Refit both initialization and refinement if switching the actual baseline |
| Weighted-saddle and random-saddle studies | Preserve classical calls | Rerun only if their scientific baseline is changed |
| Saddle-reference and single-saddle checks | Preserve classical calls | Keep reference computation separate from optimizer comparison |
| Supplement S3 paired examples | Preserve classical calls | Relabel displays alongside the manuscript revision |

New comparison datasets must use distinct algorithm identifiers such as
`classical_mds` and `metric_stress_mds`, and record the graph identity, grip and
backend versions, objective, scale convention, initialization, iteration budget,
and selected start. Do not let a cache key based only on the old function name
reuse a classical result for stress MDS. Changing an MDS initializer requires
rerunning downstream edge-KK; unchanged graph and geodesic references may be
reused after verification.

The Phase 1 tests check agreement with direct SMACOF calls, restoration of input
units, independently calculated raw and normalized scores, Euclidean recovery,
changes of units, multiple starts, RNG preservation, failure/termination
reporting, and initializer dispatch. They do not establish global optimality
or constitute the new manuscript experiment.


## Phase 2 experiment checkpoint

The separate [report](experiments/mds-initializer-sensitivity/report.pdf) and
[reproduction protocol](experiments/mds-initializer-sensitivity/README.md)
compare five saved clouds at the selected k, k ± 5, k=32, and k=80. Each graph
has classical MDS, six-start stress MDS, both downstream edge-KK fits, and the
original-coordinate control. The results include fixed-path and chord scores,
reference-geodesic errors, coordinate/surface diagnostics, start spread, graph
statistics, timings, and selected-graph additional-budget checks.

Edge-KK improves path error and worsens chord Stress-1 in all 25 comparisons
under each initializer. Stress initialization has lower final path error on
23 of 25 graphs; this is not uniform superiority or a recovery theorem. The
original selected classical scores reproduce within 8.3e-11. Main-manuscript
integration and radius experiments remain later phases.

The experiment also exposed tiny asymmetric distances from near-tie path
preparation. All final MDS fits use independently saved strict symmetric graph
distances; retained routes still define path scores. A three-vertex reproducer
and the general package integration limitation are documented in the experiment
README. No package code was changed in Phase 2.
