# Non-Oracle Graph Selection Method Inventory

Built 2026-05-15. Revised 2026-05-15 after audit.

## Purpose

This inventory separates the non-oracle graph-selection methods that have
appeared across SIMODS planning notes, `gflow` examples, CT-clearance workflows,
cell-cycle handoffs, and symptoms/ASV scripts.  Its purpose is to give the
Report Agent and Experiment Engineer a source-driven map of what each method
selects, what criterion it uses, which inputs it needs, which scripts already
demonstrate it, and where it belongs in the SIMODS manuscript.

The key distinction is that all of these methods select graph parameters
without using the oracle surface geodesic target.  They do not all answer the
same question.  Some select graphs that make selected feature signals smooth,
some select graphs that are structurally stable across adjacent parameters,
some select graphs for a named downstream response, and the planned synthetic
feature benchmark tests whether these criteria recover oracle-good graph
parameters on controlled surfaces.

## Summary Table

| Method family | Selects | Mathematical criterion | Main inputs | Current readiness | Recommended placement |
|---|---|---|---|---|---|
| Top-variable-feature graph-smoothing GCV | Graph parameter for smooth high-variance or high-abundance feature panels | Robust aggregate of per-feature minimum GCV over top variable/abundant features | Candidate graphs, feature matrix, top-feature rule, `fit.rdgraph.regression()` / `refit.rdgraph.regression()` or equivalent smoother | Implemented in symptoms; related CT scripts include a faster matrix-smoother variant | Primary candidate for synthetic non-oracle benchmark; main text only after selector-vs-oracle evidence |
| Gene-set/task-aligned GCV | Graph parameter for a named gene/feature set | \(J(\theta)=|F|^{-1}\sum_{g\in F}\min_\eta \mathrm{GCV}_g(\theta,\eta)\), preferably normalized and robustly aggregated | Candidate graphs, named feature set \(F\), expression/features, graph smoother | Collaborator-ready cell-cycle handoff exists; needs SIMODS-specific framing | Supplement or future-work/application method; mention in main text as objective-aligned variant |
| PC-score / latent-coordinate GCV | Graph parameter that smooths dominant latent coordinates | Aggregate GCV over PC1--PCq or other latent coordinates | Candidate graphs, PCA/latent embedding, PC response matrix | Demonstrated in CT-clearance prerun/production logic | Supplement/background; useful sensitivity check, not primary SIMODS criterion |
| Structural graph-stability selectors | Stable region of graph-parameter sequence | Jensen-Shannon divergence between graph-summary PMFs; edge symmetric difference/edit distance; composite stability scores | Ordered graph sequence, graph summaries, optional labels | Implemented in `gflow` helpers and examples | Main-text diagnostic or supplement, not sole accuracy claim |
| Connectivity and MST-repair burden | Feasible graph region / guardrail | Component counts, largest-component fraction, repair edge count or burden | Candidate graphs and repair metadata | Present in SIMODS benchmark plan and several workflows | Guardrail in main text; detailed criteria in supplement |
| Task/outcome conditional-expectation GCV | Graph parameter for a specific response/outcome | Minimum GCV for a response \(y\), often semi-supervised and component/labeled-count weighted | Candidate graphs, outcome labels, labeled vertices, graph smoother | Implemented in CT-clearance and symptoms VAG_ODOR scripts | Application/supplement; do not present as universal geodesic selector |
| Synthetic-feature GCV benchmark design | Test design for selector-vs-oracle recovery | Generate smooth/noisy/distractor features, run non-oracle selectors, compare \(\hat\theta\) to \(\theta_{\mathrm{oracle}}\) | Controlled surface data, oracle metrics, synthetic feature matrix, candidate graph grid | Planned in H003 benchmark plan | Main bridge experiment if run; currently a design, not result evidence |

## Conductance Rule and Smoother Semantics

The current `fit.rdgraph.regression()` low-pass operator is not a direct
metric-length graph smoother.  For precomputed graphs, supplied `weight.list`
values are positive edge lengths used to order and truncate local
neighborhoods.  The current mass-symmetrized spectral path then uses
Riemannian-complex edge masses \(\rho_1(e)\) derived from local-neighborhood
overlaps, with conductance
\[
  c_e^\rho = \frac{1}{\max(\rho_1(e),10^{-10})}.
\]
It should not be described as using the standard length conductance
\[
  c_e^{\mathrm{len}} = \frac{1}{\ell_e+\epsilon}.
\]

This distinction matters for non-oracle graph selection.  The existing
`rdgraph` smoother may prefer graph parameters that produce favorable
overlap-density/Riemannian-complex geometry, whereas a length-conductance
smoother would more directly test metric edge-length geometry.  The synthetic
selector-vs-oracle benchmark should therefore include a planned comparator:
current `rdgraph` overlap-density smoothing versus a direct length-conductance
graph-signal smoother.  The comparison should ask whether each non-oracle GCV
criterion selects graph parameters that agree with oracle geodesic-isometry
performance on synthetic quad benchmarks.

For source details, use
`notes/gflow_low_pass_gcv_operator/gflow_low_pass_gcv_operator.tex`, especially
the sections on the recent GCV graph-selection regime, neighborhood-overlap
edge masses, and `L.c1` conductance assignment.

## Per-Script Parameter Regimes

The source scripts do not all use the same operator regime.  Some are exact
`fit.rdgraph.regression()` / `refit.rdgraph.regression()` examples; others use
a faster row-normalized adjacency smoother as a GCV proxy.  Manuscript claims
about the exact `gflow` low-pass operator should primarily cite the exact
`fit`/`refit` scripts and the low-pass appendix, not the matrix-smoother
diagnostic scripts.

For the intended SIMODS synthetic benchmark, do not rely on package defaults in
the methods prose.  State the intended exact setting explicitly:
`use.counting.measure = TRUE`, `response.penalty.exp = 0`,
`filter.type = "heat_kernel"`, and `refit.rdgraph.regression(per.column.gcv =
TRUE)` for multi-feature scoring, unless a branch is deliberately labeled as a
different comparator.

| Script / workflow | Selector family | `use.counting.measure` | `response.penalty.exp` | `filter.type` | `per.column.gcv` | Exact `fit`/`refit` or approximation | Notes |
|---|---|---:|---:|---|---|---|---|
| `symptoms/R/11_asv_hv_k_gcv_sweep.R` | top-variable ASV GCV | default `TRUE` | default `0` | default `heat_kernel` | `TRUE` for refit columns | exact `fit` seed + exact `refit` | Strongest direct top-feature implementation for the current operator. |
| `symptoms/R/17_fit_vag_odor_gcv_over_asv_graphs.R` | outcome GCV | default `TRUE` | default `0` | default `heat_kernel` | not used | exact `fit` only | Outcome-specific selector; useful as application evidence, not universal geometry evidence. |
| `cell_cycle/share/yang_k_selection_handoff/run_gene_set_condexp_gcv_k_sweep.R` | gene-set/task-aligned GCV | default `TRUE` | default `0` | default `heat_kernel` | `TRUE` for refit genes | exact `fit` seed + exact `refit` | Strong task-aligned `fit`/`refit` example. |
| `CT_clearance/scripts/59_run_gvf_rdgraph_gene_smoothing_prerun.R` graph scoring loop | PC-score and top-gene GCV | not applicable | not applicable | not applicable | not applicable | faster matrix smoother | Uses `make.smoother.matrix(adj.list)` and `score.responses()` for graph scoring; supports GCV-selection idea, not exact `rdgraph` operator claims. |
| `CT_clearance/scripts/59_run_gvf_rdgraph_gene_smoothing_prerun.R` timing fire test | exact refit timing check | `FALSE` in historical script | default `0` | `heat_kernel` | `TRUE` for 100-gene refit | exact `fit` + exact `refit`, historical/stale regime | Current `gflow` precomputed-graph wrapper requires `use.counting.measure = TRUE`, so this path should be treated as historical unless rerun against a compatible source revision. |
| `CT_clearance/scripts/61_run_gvf_rdgraph_gene_smoothing_production.R` graph scoring loop | PC-score and top-gene GCV | not applicable | not applicable | not applicable | not applicable | faster matrix smoother | Production `k` scoring uses row-normalized adjacency smoother, not exact `rdgraph` spectral smoothing. |
| `CT_clearance/scripts/61_run_gvf_rdgraph_gene_smoothing_production.R` anchor/refit | production all-gene smoothing after selected `k` | enforced `TRUE` | default `0` | `heat_kernel` | `TRUE` for all-gene refit | exact `fit` + exact `refit` after matrix-smoother selection | Useful for exact `fit`/`refit` application after selection; not the selection criterion itself. |
| `CT_clearance/scripts/66_run_ct_phylotype_iknn_1x1_absorb_dcst.R` | top-phylotype GCV and structural stability | not applicable | not applicable | not applicable | not applicable | faster matrix smoother | Important real-data graph-selection diagnostic, but it should not be cited as exact `fit.rdgraph.regression()` evidence. |
| `CT_clearance/scripts/68_fit_ct_clearance_condexp_over_phylotype_iknn_graphs.R` | outcome conditional-expectation GCV | default argument `FALSE` | default `0` | default `heat_kernel` | not used | exact component-wise `fit`, historical/stale regime | Default cached-graph regime conflicts with current wrapper guard; treat as historical unless updated/rerun. |
| `CT_clearance/scripts/73_fit_ct_clearance_condexp_ct_only_phylotype_iknn_graphs.R` | outcome conditional-expectation GCV | `TRUE` | `0` | `heat_kernel` | not used | exact semi-supervised `fit` | Current-regime outcome-specific example with `y.vertices`. |

## Method Families

### 1. Top-Variable-Feature Graph-Smoothing GCV

**What it is trying to select.**  This method selects a graph parameter
\(\theta\), such as \(k\), radius, adaptive-radius scale, or cKNN parameter, by
asking which graph best supports smooth reconstruction of a panel of highly
variable or highly abundant measured features.  It is the closest implemented
ancestor of the SIMODS synthetic-feature GCV criterion.

**Mathematical criterion.**  For a candidate graph \(G_\theta\), feature
responses \(y_g\), and graph-spectral smoother \(S_{\theta,\eta}\),
the core per-feature score is
\[
  \mathrm{GCV}_g(\theta)
  =
  \min_\eta
  \frac{\|y_g-S_{\theta,\eta}y_g\|_2^2}
       {(n-\mathrm{tr}(S_{\theta,\eta}))^2}.
\]
The symptoms workflow normalizes per-feature GCV across \(k\) using robust
per-feature scaling and aggregates by mean, median, and trimmed mean.  A
manuscript-facing version should be written as
\[
  J_{\mathrm{top}m}(\theta)
  =
  \operatorname{median}_{g\in F_{\mathrm{top}m}}
  z_g(\theta)
  \quad\text{or}\quad
  \operatorname{trimmedmean}_{g\in F_{\mathrm{top}m}} z_g(\theta),
\]
where \(z_g(\theta)\) is a per-feature normalized GCV curve.

The exact `gflow` low-pass operator is documented in
`notes/gflow_low_pass_gcv_operator/gflow_low_pass_gcv_operator.tex`: the
precomputed `weight.list` is an edge-length list; it is not a conductance list,
and it is not directly converted as \(1/(\ell_{ij}+\epsilon)\).  In the current
operator, edge lengths influence neighborhood ordering; the spectral
conductance is derived from overlap-density edge mass as
\(c_e^\rho=1/\max(\rho_1(e),10^{-10})\).

**Required inputs.**

- Ordered candidate graph family \(\{G_\theta\}\) with adjacency and positive
  edge-length lists.
- A feature matrix with aligned rows/vertices.
- A feature-panel rule, typically top 20/top 30/top 50 by variance or abundance.
- A connected-graph policy: skip disconnected graphs, restrict to the largest
  connected component, or record nearest-connected fallback.
- A per-feature normalization and aggregation policy.

**Existing source scripts/reports.**

- `/Users/pgajer/current_projects/symptoms/R/11_asv_hv_k_gcv_sweep.R`
  implements top-variance ASV panels, fits one seed feature with
  `fit.rdgraph.regression()`, refits remaining features with
  `refit.rdgraph.regression(per.column.gcv = TRUE)`, normalizes GCV per
  feature, aggregates by mean/median/trimmed mean, and bootstraps selected
  `k`.
- `/Users/pgajer/current_projects/symptoms/R/12_plot_asv_hv_k_curves.R`
  plots `mean.norm` and `median.norm` against `k`.
- `/Users/pgajer/current_projects/symptoms/R/14_asv_full_graph_hv_criteria_k_selection.R`
  repeats the high-variance criterion on full-graph families.
- `/Users/pgajer/current_projects/symptoms/R/README.md` summarizes this
  workflow and its outputs.
- `/Users/pgajer/current_projects/CT_clearance/scripts/66_run_ct_phylotype_iknn_1x1_absorb_dcst.R`
  computes top 20/30/50 phylotype GCV summaries, robust-normalized summaries,
  bootstrap `k` stability, and graph stability diagnostics.  This script uses a
  simple row-normalized adjacency smoothing matrix helper, so it supports the
  general idea of GCV-based graph selection but is not evidence for the exact
  `fit.rdgraph.regression()` spectral operator.

**Known limitations and failure modes.**

- Top-variable or top-abundance features are not a neutral view of geometry.
  They can overrepresent high-abundance, high-variance, or noisy signals.
- Selection can become feature-panel-dependent.  Disagreement between top20,
  top30, and top50 should be treated as a stability diagnostic, not hidden.
- Raw feature averaging is unsafe; normalization and robust aggregation are
  required.
- GCV tests smoothness of chosen signals on the graph, not graph-geodesic
  recovery directly.
- Existing CT top-feature scripts include fast matrix-smoothing approximations;
  manuscript claims about the exact low-pass operator should cite the
  `gflow_low_pass_gcv_operator` note and the symptoms/cell-cycle exact
  `fit`/`refit` implementations.  CT 66 can support a claim that
  GCV-style graph-selection diagnostics have been used in real workflows, not a
  claim about the exact `rdgraph` overlap-density spectral smoother.

**Manuscript placement recommendation.**  Treat this as the primary candidate
GCV selector for the synthetic non-oracle benchmark, because it matches the
SIMODS scope note's top20/top30 language.  It becomes a main-text selector only
after the selector-vs-oracle benchmark is run and shows useful recovery.  In
the main text, describe the criterion conditionally and report
selector-vs-oracle recovery when available.  Put sensitivity to top20 versus
top30/top50, normalization, and bootstrap stability in the supplement.

### 2. Gene-Set / Task-Aligned GCV

**What it is trying to select.**  This method selects graph parameters for the
specific feature set or biological task that the analysis cares about.  It is
not the same as top-variable-feature GCV: \(F\) is supplied by the scientific
question, not by variance ranking.

**Mathematical criterion.**  The prior Codex session
`/Users/pgajer/.codex/sessions/2026/02/26/rollout-2026-02-26T17-42-24-019c9c1e-18af-7a83-83fd-e449fe553537.jsonl`
states the intended criterion:
\[
  J_F(k)=\frac{1}{|F|}\sum_{g\in F}\min_{\eta_g}
  \mathrm{GCV}_g(k,\eta_g),
  \qquad
  k^\star=\arg\min_k J_F(k).
\]
For manuscript use this should be robustified:
\[
  J_F(k)=\operatorname{median}_{g\in F} z_g(k)
  \quad\text{or}\quad
  \operatorname{trimmedmean}_{g\in F} z_g(k),
\]
where \(z_g(k)\) is a normalized per-gene GCV curve.  If the curve is flat, use
a one-standard-error-style rule and choose the simpler or more stable graph.

**Required inputs.**

- Candidate graph sequence and edge-length lists.
- Expression/feature matrix aligned to graph vertices.
- A named gene/feature set \(F\) with documented provenance.
- A seed gene or response for the initial `fit.rdgraph.regression()` call.
- Per-gene GCV extraction from `refit.rdgraph.regression(per.column.gcv=TRUE)`.

**Existing source scripts/reports.**

- `/Users/pgajer/current_projects/cell_cycle/share/yang_k_selection_handoff/README.md`
  is a collaborator-ready handoff for gene-set conditional-expectation GCV
  `k` selection.
- `/Users/pgajer/current_projects/cell_cycle/share/yang_k_selection_handoff/run_gene_set_condexp_gcv_k_sweep.R`
  loads a target gene set, fits a seed gene for each `k`, refits remaining
  genes with per-column GCV, writes raw/normalized per-gene GCV tables, and
  selects `k` by mean or median normalized GCV.
- `/Users/pgajer/current_projects/cell_cycle/share/yang_k_selection_handoff/plot_gene_set_condexp_gcv_k_sweep.R`
  creates per-`k` normalized-GCV histograms and summary curves.
- The CT-clearance conditional expectation scripts are response-aligned
  relatives, but they belong more directly to method family 6 below.

**Known limitations and failure modes.**

- This is objective-aligned, not universally geometry-aligned.  Different
  scientifically valid feature sets can select different graph parameters.
- It can overfit the chosen feature panel unless checked by bootstrap,
  feature-subset resampling, held-out/cross-fit evaluation, or a simple-graph
  rule.
- Comparisons across graph families require the same smoother, feature panel,
  normalization, eta policy, and graph feasibility rules.
- A very small or highly coexpressed feature set may produce an unstable or
  overly narrow criterion.

**Manuscript placement recommendation.**  Mention as an objective-aligned
variant of graph-smoothing GCV, but keep it in supplement/application/future
work unless a real-data SIMODS case study is explicitly organized around a
named feature set.  Do not collapse it into top20/top30 variable-feature GCV.

### 3. PC-Score / Latent-Coordinate GCV

**What it is trying to select.**  This method selects a graph that smoothly
represents dominant latent coordinates, usually PC scores.  It asks whether the
candidate graph supports the main variance axes of the data matrix.

**Mathematical criterion.**  For PC score columns
\(\mathrm{PC}_1,\ldots,\mathrm{PC}_q\),
\[
  J_{\mathrm{PC}}(\theta)
  =
  \operatorname{median}_{r=1}^{q}
  z_{\mathrm{PC}_r}(\theta),
\]
or the mean-normalized analogue used in source scripts.

**Required inputs.**

- Candidate graph sequence.
- PCA or latent embedding with documented transform, variance target, and cap.
- A choice of response PCs, usually PC1--PC5.
- Connectivity/major-component policy.

**Existing source scripts/reports.**

- `/Users/pgajer/current_projects/CT_clearance/scripts/59_run_gvf_rdgraph_gene_smoothing_prerun.R`
  computes PCA responses, scores candidate iKNN graphs with PC1--PC5 and top
  log-coverage gene panels, and records selected `k`.
- `/Users/pgajer/current_projects/CT_clearance/manifests/runs/gvf_rdgraph_gene_smoothing_prerun_2026-04-20.md`
  records that PC1--PC5, top20, and top30 selectors were evaluated; in that run
  they agreed for key branches, while PCA variance targets did not resolve
  under the 100-PC cap.
- `/Users/pgajer/current_projects/CT_clearance/scripts/61_run_gvf_rdgraph_gene_smoothing_production.R`
  continues this production logic with fixed-PC and variance-target branches.
- `/Users/pgajer/current_projects/CT_clearance/manifests/runs/gvf_rdgraph_gene_smoothing_production_2026-04-20.md`
  records selected `k`, nearest-connected fallbacks, and major-component
  policies.

**Known limitations and failure modes.**

- PC scores are derived from the same data used to build the graph, so the
  criterion can be circular if presented as independent validation.
- PC-score GCV tends to select graphs preserving dominant variance directions,
  which may not equal intrinsic geodesic-distance recovery.
- PCA caps and unresolved variance targets make interpretation fragile.  The
  CT prerun explicitly warns that unresolved PCA targets should be treated as
  fixed-cap sensitivity branches.
- PC score smoothness may be too easy on synthetic examples if the graph was
  built from the same coordinates.

**Manuscript placement recommendation.**  Treat as a sensitivity/diagnostic
method in supplement.  It can be useful in real-data workflows, but it should
not be the main SIMODS non-oracle graph-selection argument.

### 4. Structural Graph-Stability Selectors

**What they are trying to select.**  These selectors choose a stable region of a
candidate graph sequence without using feature responses.  They are graph-only
diagnostics of whether nearby parameter values produce similar graph summaries
or edge sets.

**Mathematical criteria.**  The current `gflow` implementation supports graph
summary PMFs \(P_\theta\), including degree distribution, edge-weight
distribution, component-size distribution, and neighborhood-label distribution.
For adjacent parameters,
\[
  \mathrm{JS}_{\mathrm{summary}}(\theta_j)
  =
  \mathrm{JS}(P_{\theta_j},P_{\theta_{j+1}}),
\]
where `jensen.shannon.divergence()` uses base-2 Jensen-Shannon divergence.

For edge-set stability, the iKNN helper computes the unweighted edge symmetric
difference
\[
  \mathrm{Edit}(\theta_j)
  =
  |E_{\theta_j}\setminus E_{\theta_{j+1}}|
  +
  |E_{\theta_{j+1}}\setminus E_{\theta_j}|.
\]
The older general `graph.edit.distance()` also supports weighted edit distance
with edge insertion/deletion cost and a weight-difference cost.

`find.optimal.k.from.stability()` combines scaled badness in edit distance and
JS divergence with scaled edge-count stability into a score where larger is
better.  The SIMODS benchmark plan additionally proposes selecting the smallest
eligible parameter within a tolerance of the minimum JS/edit score.

**Required inputs.**

- Ordered graph sequence with comparable vertex sets.
- Parameter ordering and, for graph-family-general use, a common way to compare
  adjacent parameters by native scale or edge count.
- Optional labels when using neighborhood-label distribution.

**Existing source scripts/reports.**

- `/Users/pgajer/current_projects/gflow/R/graph_summary_divergence.R` defines
  `compute.graph.summary.pmf()`, `graph.summary.divergence()`, and
  `compute.graph.summary.stability()`.
- `/Users/pgajer/current_projects/gflow/R/divergences.R` defines
  `jensen.shannon.divergence()`.
- `/Users/pgajer/current_projects/gflow/R/iknn_graphs.R` defines
  `compute.stability.metrics()`, `compute.edit.distances()`, and
  `find.optimal.k()`.
- `/Users/pgajer/current_projects/gflow/R/graph_edit_distance.R` defines the
  older weighted `graph.edit.distance()` and `calculate.edit.distances()`.
- `/Users/pgajer/current_projects/gflow_examples/nyc_taxi_regimes/scripts/taxi_workflow.R`
  uses degree-distribution JS stability as the primary package-level
  graph-selection criterion, with borough-label divergence, graph edit
  distance, and component checks as diagnostics.
- `/Users/pgajer/current_projects/gflow_examples/retina_cell_cycle/scripts/retina_workflow.R`
  uses degree-distribution stability with age ordering, pseudotime continuity,
  subtype coherence, graph edit distance, and connectivity guardrails.
- `/Users/pgajer/current_projects/CT_clearance/scripts/66_run_ct_phylotype_iknn_1x1_absorb_dcst.R`
  reports adjacent-`k` graph edit distance and degree-distribution JS
  diagnostics.

**Known limitations and failure modes.**

- Structural stability is not geodesic accuracy.  A graph sequence can become
  stable only after it is too dense or too sparse.
- Degree-distribution JS ignores where edges occur.
- Edge symmetric difference can be overly sensitive to local edge swaps and
  parameter-grid spacing.
- Neighborhood-label stability requires trusted labels and can become a
  supervised criterion.
- Structural criteria can disagree with GCV for valid reasons.  The benchmark
  should report disagreement rather than force a single universal score.

**Manuscript placement recommendation.**  Use structural stability as a
main-text diagnostic alongside GCV if the selector-vs-oracle results support
it; otherwise put detailed curves in supplement.  It should not be presented as
a proven proxy for geodesic accuracy without oracle comparison.

### 5. Connectivity and MST-Repair Burden

**What it is trying to select.**  Connectivity and repair criteria identify
feasible graph regions and graph-construction failure modes.  They are
guardrails, not accuracy criteria by themselves.

**Mathematical criteria.**

- Component count: \(c(G_\theta)\).
- Largest-component fraction:
  \(\mathrm{LCC}(\theta)=|V_{\max}(G_\theta)|/|V|\).
- MST-repair burden: number or fraction of repair edges added to make a graph
  connected.
- Ineligibility rule such as \(c(G_\theta)>1\) or
  \(\mathrm{LCC}(\theta)<0.99\), depending on the workflow.

**Required inputs.**

- Candidate graphs before and after repair, or graph-construction metadata
  recording repair edges.
- Vertex count and component assignments.
- A declared policy: skip, repair, largest-component fit, nearest-connected
  fallback, or major-component-qualified selection.

**Existing source scripts/reports.**

- `/Users/pgajer/current_projects/grip/papers/data_to_geodesic_embedding_paper/notes/non_oracle_graph_selection_benchmark_plan.md`
  treats connectivity and MST-repair burden as a guardrail for selector
  eligibility and failure-mode reporting.
- `/Users/pgajer/current_projects/gflow_examples/nyc_taxi_regimes/scripts/taxi_workflow.R`
  applies a largest-component guardrail of at least 0.90 when selecting by
  degree stability.
- `/Users/pgajer/current_projects/gflow_examples/retina_cell_cycle/scripts/retina_workflow.R`
  records component stability, `n_components`, and `lcc_frac` as selected-`k`
  guardrails.
- `/Users/pgajer/current_projects/CT_clearance/manifests/runs/gvf_rdgraph_gene_smoothing_production_2026-04-20.md`
  records major-component-qualified selected `k` and nearest-connected `k`
  used for production fitting when the selected graph was disconnected.

**Known limitations and failure modes.**

- Connectivity is necessary for many smoothers but not sufficient for geometry.
- MST repair can dominate sparse graphs; if many repair edges are added, the
  selected graph may mostly reflect the repair policy.
- Largest-component fitting changes the estimand by dropping vertices.
- Component policies can affect comparability across graph families.

**Manuscript placement recommendation.**  Include as a main-text guardrail in
the non-oracle benchmark design and report repair burden in selector tables.
Detailed component summaries belong in supplement.

### 6. Task / Outcome Conditional-Expectation GCV

**What it is trying to select.**  This method selects the graph parameter that
best supports a specific outcome or response, such as CT clearance or VAG_ODOR.
It is useful for applications, but it is not a universal graph-geodesic
criterion.

**Mathematical criterion.**  For a response \(y\), possibly observed only on
labeled vertices \(L\),
\[
  J_y(\theta)=\min_\eta \mathrm{GCV}_L(y;G_\theta,\eta).
\]
In CT-clearance scripts, component-level GCV values may be weighted by labeled
sample count:
\[
  J_y(\theta)
  =
  \frac{\sum_c n_{L,c}\,\mathrm{GCV}_{c}(\theta)}
       {\sum_c n_{L,c}}.
\]

**Required inputs.**

- Candidate graph sequence.
- A response vector with label mask or complete observed response.
- Component/largest-component policy for disconnected graphs.
- Fixed `fit.rdgraph.regression()` parameters and GCV extraction policy.

**Existing source scripts/reports.**

- `/Users/pgajer/current_projects/CT_clearance/scripts/68_fit_ct_clearance_condexp_over_phylotype_iknn_graphs.R`
  estimates \(E[\mathrm{CT\ clearance}\mid G]\) over CT+VIRGO2 phylotype iKNN
  graphs, fits component-wise graph regression, and writes weighted GCV by `k`.
- `/Users/pgajer/current_projects/CT_clearance/scripts/73_fit_ct_clearance_condexp_ct_only_phylotype_iknn_graphs.R`
  performs direct CT-only CT-clearance fits over `k=3:15`, restricting
  disconnected graphs to a major component before calling
  `fit.rdgraph.regression()`.
- `/Users/pgajer/current_projects/CT_clearance/scripts/67_generate_ct_phylotype_iknn_1x1_weighted_grip_dcst_html.R`
  includes the CT-clearance GCV plot in an inspection gallery when available.
- `/Users/pgajer/current_projects/CT_clearance/scripts/70_compare_ct_phylotype_graph_dcst_congruence.R`
  compares dCST congruence, edge homophily, and CT-clearance weighted GCV
  across relative-abundance and Hellinger graph branches.
- `/Users/pgajer/current_projects/symptoms/R/17_fit_vag_odor_gcv_over_asv_graphs.R`
  fits `fit.rdgraph.regression()` for `VAG_ODOR` over ASV graph families and
  `k` values.
- `/Users/pgajer/current_projects/symptoms/R/18_plot_vag_odor_gcv_vs_k.R`
  plots VAG_ODOR GCV by `k` and writes minima by family and overall.

**Known limitations and failure modes.**

- The selected graph is response-optimal, not necessarily geometry-optimal.
- Semi-supervised or component-wise fitting can make scores hard to compare
  when labeled counts differ by component or graph parameter.
- A response-specific optimum may conflict with top-feature GCV or structural
  stability for valid reasons.
- It can encourage overfitting to one outcome if not checked on held-out labels
  or alternate outcomes.

**Manuscript placement recommendation.**  Present as an application-oriented
extension or real-data case-study criterion, not as the core SIMODS bridge from
oracle surfaces to real data.  It can be used in supplement to show that
non-oracle graph selection can be adapted to specific downstream tasks.

### 7. Synthetic-Feature GCV Benchmark Design

**What it is trying to select.**  This is not yet a completed result.  It is the
planned benchmark design for testing whether GCV selectors recover
oracle-good graph parameters on controlled quadratic surfaces.

**Mathematical criterion.**  Generate feature signals \(Y\) without using the
oracle distance matrix, run a GCV selector over candidate graphs, and compare
the selected parameter to the oracle-best parameter:
\[
  \hat\theta_{\mathrm{GCV}}
  =
  \arg\min_\theta J_Y(\theta),
  \qquad
  \theta_{\mathrm{oracle}}
  =
  \arg\min_\theta \mathrm{Err}_{\mathrm{oracle}}(G_\theta).
\]
Evaluation should report oracle rank, oracle error ratio
\[
  \mathrm{Err}_{\mathrm{oracle}}(G_{\hat\theta})/
  \min_\theta \mathrm{Err}_{\mathrm{oracle}}(G_\theta),
\]
near-oracle success within 5% and 10%, selected edge-count ratio, repair
burden, and seed stability.

A required methodological fork is the smoothing conductance rule:
\[
  c_e^\rho=\frac{1}{\max(\rho_1(e),10^{-10})}
  \quad\text{versus}\quad
  c_e^{\mathrm{len}}=\frac{1}{\ell_e+\epsilon}.
\]
The first is the current `fit.rdgraph.regression()` overlap-density /
Riemannian-complex convention.  The second is a planned direct
length-conductance graph-signal comparator.  They answer related but distinct
questions: overlap-density smoothing may select graphs with stable local
neighborhood-overlap geometry, while length-conductance smoothing more directly
tests whether the metric edge lengths support smooth graph signals.  The
benchmark should report whether these two non-oracle selectors agree with each
other and with oracle geodesic-isometry performance.

**Required inputs.**

- Controlled quadratic-surface datasets with oracle graph-geodesic metrics.
- Candidate graph families and parameter grids.
- Synthetic feature matrix containing ambient coordinates, smooth polynomial
  signals, noisy smooth signals, and distractor/noise features.
- A top20/top30 variable-feature selection rule used as a stress test.
- Both an exact/current `rdgraph` overlap-density smoother and a planned direct
  length-conductance smoother, or an explicit note if the length-conductance
  comparator is deferred.
- Output schema joining oracle metrics, non-oracle selector scores, and
  selector choices.

**Existing source scripts/reports.**

- `/Users/pgajer/current_projects/grip/papers/data_to_geodesic_embedding_paper/notes/non_oracle_graph_selection_benchmark_plan.md`
  defines the current planned design, including feature blocks, selector-vs-
  oracle metrics, smoke/full runs, and missing implementation pieces.
- `/Users/pgajer/current_projects/grip/papers/data_to_geodesic_embedding_paper/codex_project/project_status/handoffs/H003_experiment-engineer_non-oracle-selection-benchmark-plan.md`
  records the Experiment Engineer's answered plan and missing implementation
  pieces.
- `/Users/pgajer/current_projects/grip/papers/data_to_geodesic_embedding_paper/notes/simods_data_geodesic_mds_paper_scope.md`
  identifies non-oracle graph selection as the major missing bridge from oracle
  benchmarks to real data.

**Known limitations and failure modes.**

- Ambient coordinates may make the task too easy.
- GCV may select graphs that smooth synthetic signals well but have poor
  graph-geodesic recovery.
- The current `rdgraph` smoother and a length-conductance smoother can
  disagree for principled reasons, because they test overlap-density geometry
  versus edge-length metric geometry.
- Noise and distractor levels can change selected parameters.
- Oracle error curves may be flat; near-oracle success rates are more
  informative than exact parameter matches.
- Selector disagreement is expected and should be reported as evidence, not
  treated as a bookkeeping error.

**Manuscript placement recommendation.**  If implemented, this is the main
manuscript bridge experiment for non-oracle graph parameter selection.  Until
then it is planning evidence, not result evidence.  Main-text wording should
call top-variable/synthetic-feature GCV a primary candidate selector, not a
validated selector, until the selector-vs-oracle benchmark has been run.

## Cross-Method Limitations

- All GCV variants measure signal smoothability on the candidate graph, not
  geodesic accuracy directly.
- Feature-panel choice is part of the method.  Top-variable, gene-set,
  PC-score, and outcome-specific panels are distinct selectors.
- Per-feature GCV must be normalized before aggregation unless all responses
  are on a controlled common scale.
- Disconnected graphs require explicit policy, because `fit.rdgraph.regression()`
  expects connected precomputed graphs.
- Graph-family comparisons are only meaningful if the same smoother, eta-grid
  policy, feature panel, graph stage, and feasibility rules are used.
- Exact `fit.rdgraph.regression()` results and faster matrix-smoother GCV
  diagnostics are not interchangeable.  Use symptoms/cell-cycle exact
  `fit`/`refit` scripts and the low-pass appendix for claims about the current
  operator; use CT 59/61/66 matrix-smoother branches only as evidence for the
  broader GCV-selection design pattern.
- The current `rdgraph` smoother uses overlap-density conductance, not direct
  length conductance.  Any claim about metric edge-length smoothing requires a
  separate length-conductance comparator or a clearly marked planned method.
- Structural stability and GCV can disagree for scientifically meaningful
  reasons; the benchmark should compare both to oracle error, not force a
  single proxy to win everywhere.

## Recommended Manuscript Placement

For the SIMODS main text, the cleanest story is:

1. Define the exact `gflow` low-pass GCV smoother and precomputed edge-length
   convention, including the overlap-density conductance caveat.
2. Run top-variable/synthetic-feature GCV as a primary candidate non-oracle
   selector on controlled surfaces.
3. Compare GCV-selected parameters to oracle-best parameters.
4. Report JS degree stability, edge symmetric-difference/edit stability, and
   connectivity/MST repair burden as diagnostics and guardrails.
5. If feasible, compare the current overlap-density `rdgraph` smoother with a
   direct length-conductance smoother as a planned metric-geometry comparator.

The supplement can then hold:

- top20/top30/top50 sensitivity;
- gene-set/task-aligned GCV;
- PC-score GCV;
- outcome-specific CT-clearance and VAG_ODOR examples;
- detailed structural-stability and component-policy tables.

The manuscript should not claim that any non-oracle criterion is a universal
surrogate for surface geodesic accuracy until the selector-vs-oracle benchmark
has been run.

## Source Asset Index

| Source asset | Status | Method families supported | Notes |
|---|---|---|---|
| `/Users/pgajer/current_projects/grip/papers/data_to_geodesic_embedding_paper/notes/non_oracle_graph_selection_benchmark_plan.md` | exists | synthetic-feature design, top-variable GCV, structural stability, connectivity | Current H003 benchmark plan and missing implementation list |
| `/Users/pgajer/current_projects/grip/papers/data_to_geodesic_embedding_paper/notes/simods_data_geodesic_mds_paper_scope.md` | exists | manuscript placement, non-oracle bridge | Defines non-oracle graph selection as essential missing bridge |
| `/Users/pgajer/current_projects/grip/papers/data_to_geodesic_embedding_paper/notes/gflow_low_pass_gcv_operator/gflow_low_pass_gcv_operator.tex` | exists | all GCV methods | Exact `fit.rdgraph.regression()` / `refit.rdgraph.regression()` operator note |
| `/Users/pgajer/current_projects/gflow/R/riem_dcx_regression.R` | exists | graph smoothing GCV | Current source for `fit.rdgraph.regression()` |
| `/Users/pgajer/current_projects/gflow/R/refit_rdgraph_regression.R` | exists | multi-feature refit and per-column GCV | Current source for `refit.rdgraph.regression()` |
| `/Users/pgajer/current_projects/gflow/R/data_smoother.R` | exists | package-level smoothing/GCV workflow | Uses GCV over `k` and refit on full matrix |
| `/Users/pgajer/current_projects/gflow/R/graph_summary_divergence.R` | exists | structural stability | Graph-summary PMFs and JS stability |
| `/Users/pgajer/current_projects/gflow/R/iknn_graphs.R` | exists | iKNN structural stability | `compute.stability.metrics()`, edit distance, `find.optimal.k()` |
| `/Users/pgajer/current_projects/gflow/R/graph_edit_distance.R` | exists | structural graph edit | Older weighted edit-distance helper |
| `/Users/pgajer/current_projects/gflow/R/rdgraph_regression.R` | missing | stale handoff path | No file by this name found; use `riem_dcx_regression.R` |
| `/Users/pgajer/current_projects/gflow/R/rdgraph_smoothing.R` | missing | stale handoff path | No file by this name found; use `data_smoother.R` and refit source |
| `/Users/pgajer/current_projects/symptoms/R/11_asv_hv_k_gcv_sweep.R` | exists | top-variable GCV | Strongest direct top-feature `fit`/`refit` example |
| `/Users/pgajer/current_projects/symptoms/R/17_fit_vag_odor_gcv_over_asv_graphs.R` | exists | outcome GCV | VAG_ODOR response-specific GCV |
| `/Users/pgajer/current_projects/CT_clearance/scripts/59_run_gvf_rdgraph_gene_smoothing_prerun.R` | exists | PC-score GCV, top-gene GCV | Graph scoring uses faster matrix smoother; timing fire test uses exact `fit`/`refit` with historical `use.counting.measure=FALSE` |
| `/Users/pgajer/current_projects/CT_clearance/scripts/61_run_gvf_rdgraph_gene_smoothing_production.R` | exists | PC-score GCV production branch | Matrix-smoother selection followed by exact `fit`/`refit` production smoothing with counting measure TRUE |
| `/Users/pgajer/current_projects/CT_clearance/scripts/66_run_ct_phylotype_iknn_1x1_absorb_dcst.R` | exists | top-feature GCV, structural stability | Faster matrix-smoother diagnostic; do not cite as exact `rdgraph` operator evidence |
| `/Users/pgajer/current_projects/CT_clearance/scripts/68_fit_ct_clearance_condexp_over_phylotype_iknn_graphs.R` | exists | outcome GCV | Component-wise CT-clearance exact fits; default `use.counting.measure=FALSE` is historical/stale relative to current wrapper guard |
| `/Users/pgajer/current_projects/CT_clearance/scripts/73_fit_ct_clearance_condexp_ct_only_phylotype_iknn_graphs.R` | exists | outcome GCV | Direct major-component CT-only fits |
| `/Users/pgajer/current_projects/cell_cycle/share/yang_k_selection_handoff/README.md` | exists | gene-set GCV | Collaborator-ready task-aligned `k` selection handoff |
| `/Users/pgajer/current_projects/cell_cycle/share/yang_k_selection_handoff/run_gene_set_condexp_gcv_k_sweep.R` | exists | gene-set GCV | Run script for target gene-set sweep |
| `/Users/pgajer/current_projects/gflow_examples/nyc_taxi_regimes/scripts/taxi_workflow.R` | exists | structural stability | Example primary criterion: degree-distribution stability with guardrails |
| `/Users/pgajer/current_projects/gflow_examples/retina_cell_cycle/scripts/retina_workflow.R` | exists | structural stability, biological guardrails | Degree-stability primary, age/pseudotime/subtype checks |

## Open Questions / Needs Audit

- The missing `R/rdgraph_regression.R` and `R/rdgraph_smoothing.R` rows in the
  source asset index are retained only as provenance for stale handoff paths.
  Report Agent should not cite them as live `gflow` sources.
- The Experiment Engineer should choose one canonical GCV interface for the
  synthetic benchmark: exact `fit.rdgraph.regression()` plus
  `refit.rdgraph.regression(per.column.gcv=TRUE)` is preferred for manuscript
  consistency, but a faster Laplacian smoother may be needed for smoke tests.
- The Experiment Engineer should decide whether the length-conductance
  comparator is implemented in the first synthetic benchmark or marked as a
  planned follow-up.  Without that comparator, the benchmark tests the current
  overlap-density smoother, not direct metric-length graph smoothing.
- The benchmark needs a graph-family-general way to order adjacent parameters
  for JS/edit stability.  Native parameter order may not be comparable across
  adaptive-radius, cKNN, sKNN, and iKNN.
- The resource handoff listed two stale `gflow` R source paths:
  `R/rdgraph_regression.R` and `R/rdgraph_smoothing.R`; these should not be used
  as canonical citations.
- The CT scripts include exact `gflow` graph regression, historical cached-graph
  regimes with `use.counting.measure = FALSE`, and simplified smoother-matrix
  GCV variants.  The inventory separates them, but a final methods section must
  not conflate them.
- The Report Agent should not write result claims for the synthetic-feature
  selector benchmark until experiments exist and are joined to oracle metrics.
