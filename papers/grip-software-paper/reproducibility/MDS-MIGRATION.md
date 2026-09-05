# Classical versus stress MDS: Phase 1 migration

From grip 0.2.0.9000, the former `metric.mds()` implementation is named
`classical.mds()`. The new `metric.mds()` minimizes raw distance stress through
ratio SMACOF. Existing manuscript reproduction calls have been migrated to
`classical.mds()` to preserve their algorithm. The edge-KK default likewise
remains classical, with explicit name `init = "classical_mds"`.

This migration does not replace saved fits or claim that a stress-MDS comparison
has already been performed. Older labels such as `Metric MDS`, `metric-MDS`,
and serialized fields such as `metric_mds` refer to classical scaling in these
existing artifacts. The main manuscript and supplement will receive a complete
figure/label revision when the paired baseline comparison is incorporated.

| Existing workflow | Phase 1 action | Later comparison |
|---|---|---|
| Five-cloud two-fidelity saddle pilot | Preserve graphs, references, selected classical fits, and edge-KK results; migrate `fit-layouts.R` | Add stress MDS and its own edge-KK refinement, then nearby-k sensitivity |
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
