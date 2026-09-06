# SIMODS Paper Scope And Outline

Build timestamp: 2026-05-14

## Proposed Paper Home

The proposed SIMODS paper should live in:

```text
/Users/pgajer/current_projects/geodesicMDS
```

The reason is conceptual and practical.  The paper's endpoint is a geodesic
MDS method, and the existing manuscript/report assets for that method already
live in `geodesicMDS`.  The companion project

```text
/Users/pgajer/current_projects/geodesic_data_geometry
```

should be treated as the source of the data-to-graph reconstruction layer:
graph construction methods, oracle quadratic-surface benchmarks, non-oracle
graph-selection diagnostics, and data geometry reports.  The resulting paper is
not a merger of repositories.  It is a single SIMODS manuscript whose canonical
manuscript home is `geodesicMDS`, with `geodesic_data_geometry` supplying a
major part of the experimental and methodological foundation.

A natural manuscript path is:

```text
/Users/pgajer/current_projects/grip/papers/data_to_geodesic_embedding_paper/
```

This note is a scope document and should be superseded by a manuscript-specific
planning file once the SIMODS manuscript directory is created.

## Working Title

Possible titles:

- Data Geodesic Geometry And Geodesic Multidimensional Scaling
- From Data Graphs To Geodesic MDS
- Geodesic Quasi-Isometric Embeddings From Data-Derived Graphs
- Edge-Isometric Graph Layouts For Data Geodesic Geometry

The first title is currently the strongest because it signals that the paper
contains both pieces of the pipeline:

```text
data -> data-derived weighted graph -> graph geodesic metric -> geodesic MDS embedding
```

## Central Thesis

Classical MDS takes a dissimilarity matrix and constructs an embedding whose
ambient Euclidean distances approximate those dissimilarities.  In manifold and
structured-data settings, the relevant geometry is often not an ambient metric
but a geodesic metric.  A data-analysis workflow therefore needs two linked
problems:

1. **Data geodesic geometry reconstruction**: construct a weighted graph
   `G(X)` from a finite dataset `X` such that graph geodesic distances
   approximate intrinsic/geodesic distances on the sampled object.

2. **Geodesic MDS**: construct an embedding of the vertices of a weighted graph
   such that the graph geodesic metric is represented by embedded graph path
   lengths, not merely by all-pairs ambient chords.

The combined pipeline gives a route from data to a geodesic quasi-isometric
embedding:

```text
X subset R^p
  -> weighted graph G(X)
  -> graph geodesic distances d_G
  -> embedding Z with d_G^Z approximately d_G
```

This framing is important.  A standalone geodesic-MDS paper can look detached
from data analysis because its input is a graph.  A standalone data-graph paper
can look incomplete because it stops at a graph metric.  The SIMODS paper
should present the two as one scientific workflow.

## Mathematical Problems

### Data Geodesic Geometry Reconstruction

Let `X = {x_1, ..., x_n}` be a finite sample from a structured geometric object
`Gamma`, with intrinsic/geodesic distance `d_Gamma`.  A weighted graph
construction produces

```text
G(X) = (V, E, ell)
```

where `V = {1, ..., n}` and edge lengths `ell_ij > 0` are assigned for
`(i,j) in E`.  The induced graph geodesic distance is

```text
d_ij^G = shortest path distance in G(X).
```

The goal is to make `d_ij^G` close to the sampled intrinsic distances
`d_Gamma(x_i, x_j)`, up to a controlled notion of relative distortion.

For oracle examples such as quadratic surfaces, `d_Gamma` is approximated by a
high-resolution reference geodesic oracle.  For real data, `d_Gamma` is
unknown, so graph parameters must be selected by non-oracle stability,
smoothness, and graph-statistical criteria.

### Geodesic MDS

Given a weighted graph `G = (V, E, ell)`, seek an embedding

```text
Z = (z_1, ..., z_n)^T in R^{n x q}
```

usually with `q = 2` or `3`.  Embedded edge lengths are

```text
d_ij^Z = ||z_i - z_j||_2,   (i,j) in E.
```

Embedded graph path lengths induce a graph-on-layout distance

```text
d_ij^{G,Z} = shortest path distance in G using edge lengths d_ij^Z.
```

The geodesic MDS target is

```text
d_ij^{G,Z} approximately d_ij^G.
```

If all edge lengths are preserved, `d_ij^Z = ell_ij` for all edges, then all
graph path lengths and all graph geodesic distances are preserved exactly.
This is the key bridge from graph drawing to metric representation: edge
length fidelity is not merely aesthetic.  It is sufficient for graph-geodesic
isometry.

## Main Algorithmic Story

The current practical algorithms are:

```text
weighted GRIP -> edge-KK
metric MDS    -> edge-KK
```

where edge-KK is the edge-restricted Kamada--Kawai repair operator.  Metric MDS
uses the graph geodesic distance matrix as an initialization.  Weighted GRIP is
a graph-native initialization that avoids all-pairs dense MDS.

The edge-KK objective can be written schematically as

```text
sum_{(i,j) in E} k_ij (||z_i-z_j||_2 - ell_ij)^2,
```

with density-aware stiffness schedules that first emphasize common edge-length
ranges and then mix toward uniform edge weights.  This is equivalent to the
quadratic edge-barrier objective under the convention

```text
a_ij = k_ij ell_ij^2
```

when the barrier is written in terms of ratios

```text
s_ij(Z) = d_ij^Z / ell_ij.
```

The current timing experiments indicate that edge-KK itself is fast with the
C++ backend.  For large `n`, the expensive part of the metric-MDS pipeline is
the all-pairs graph distance computation, not the edge repair.  Weighted GRIP
therefore becomes the more scalable default initialization unless a specific
experiment requires metric-MDS initialization.

## Proposed Paper Outline

### 1. Introduction

Motivate structured data whose geometry is not well represented by ambient
Euclidean distances.  Explain the gap between classical MDS and the desired
geodesic quasi-isometric embedding.

Key point: the paper is about a full pipeline from data to geodesic embedding,
not only about a graph drawing method.

### 2. Related Work

Include:

- classical and metric MDS;
- Isomap and graph-geodesic manifold learning;
- PHATE and diffusion geometry;
- kNN, mutual kNN, symmetric/union kNN, adaptive-radius, continuous kNN;
- graph drawing, Kamada--Kawai, Fruchterman--Reingold, GRIP;
- graph signal smoothing and graph parameter selection;
- manifold reconstruction and geodesic distance approximation.

### 3. Data Geodesic Geometry Reconstruction

Formalize the problem:

```text
X -> G(X),    d_G(X) approximately d_Gamma on X.
```

Define graph families:

- symmetric kNN / union kNN;
- mutual kNN;
- intersection kNN if retained as a baseline;
- adaptive-radius graphs;
- continuous kNN / geometric-mean adaptive-radius graphs;
- Delaunay 1-skeleton reference graphs where applicable;
- PHATE graph as a related adaptive affinity graph, not necessarily a primary
  geodesic graph construction.

Explain MST connectivity repair and geometric pruning.  If pruning is not a
main paper result, present it as a sparsification tool that preserves geometry
when the pruning criterion is satisfied.

### 4. Oracle Benchmarks On Quadratic Surfaces

Use paraboloid and saddle surfaces as controlled examples.  Include curvature
families and sample-size sweeps.  Report graph geodesic distortion relative to
surface geodesic oracles.

Likely message:

- adaptive-radius and cKNN are the strongest current graph families;
- fixed-radius and mKNN are weaker in current tests and can be sparse
  baselines rather than central methods;
- pruning mostly reduces edge count without substantially changing isometry;
- sample-oracle target was diagnostic but biased enough that surface target is
  the cleaner headline metric.

### 5. Non-Oracle Graph Parameter Selection

This section is essential for the SIMODS paper.  Oracle graph selection is
useful for method development, but real data do not provide intrinsic
geodesic-distance truth.

Candidate non-oracle selection criteria include gflow methods already developed
for iKNN and adaptable to all graph families:

- generalized cross-validation (GCV) of graph smoothing on the top 20 or 30
  most variable features;
- Jensen--Shannon distance between degree distributions of consecutive graphs;
- edit-distance or edge-symmetric-difference criteria between consecutive
  graphs;
- stability of graph statistics across nearby parameter values;
- connectivity diagnostics and MST-repair burden;
- graph signal smoothness under the graph Laplacian;
- stability of layout diagnostics under subsampling;
- agreement between independently motivated graph families such as
  adaptive-radius and cKNN.

The GCV criterion is especially attractive because it introduces a graph
smoothing operator.  That makes the graph-selection section stronger: the graph
is selected not only by combinatorial stability, but also by how well it
supports smooth representation of measured features on the data geometry.

### 6. Geodesic MDS

Define the graph-to-layout problem.  Distinguish:

- classical/metric MDS on `d_G`, which approximates graph distances by ambient
  chords;
- GMDS, which evaluates graph distances induced by embedded edge lengths;
- edge-KK, which repairs edge lengths and thereby preserves graph geodesics
  when successful.

Present:

```text
metric MDS -> edge-KK
weighted GRIP -> edge-KK
```

as the two main algorithms.

### 7. GMDS Diagnostics And Stopping Rules

Define:

- edge relative RMS error;
- q90/q95 edge stretch residuals;
- graph-geodesic stress using embedded edge lengths;
- signed bias and shortcut fraction;
- spread measures such as covariance trace, covariance volume, pairwise spread;
- convergence/stopping criteria for edge-KK.

The stopping rule should be practical:

```text
stop if edge_rRMSE <= tolerance
or if relative improvement over a window is below threshold
or if max iterations are reached.
```

The paper should explicitly separate:

- edge fidelity as the primary isometry criterion;
- spread/readability as a secondary visualization criterion;
- runtime and memory as scalability criteria.

### 8. Experiments

Minimum experiment groups:

1. Oracle graph reconstruction on quadratic surfaces.
2. Non-oracle graph parameter selection on the same oracle examples, showing
   how well selection criteria recover oracle-good parameters.
3. GMDS layout comparison:
   `metric MDS`, `weighted GRIP`, `metric MDS -> edge-KK`,
   `weighted GRIP -> edge-KK`.
4. Scaling/timing comparison through at least `n = 3200` on adaptive-radius and
   cKNN.
5. One or more real data case studies.

Potential real datasets already available or plausible:

- American Gut Project (AGP);
- cell-cycle single-cell data;
- vaginal microbiome 16S Valencia 13k training set;
- VIRGO2 2600 metagenomic dataset;
- any existing preliminary cell-cycle or microbiome examples with interpretable
  trajectories/gradients.

### 9. Discussion

Emphasize the conceptual bridge:

- MDS wants quasi-isometric embeddings of data;
- graph drawing wants readable drawings of graphs;
- this paper connects them through graph-geodesic isometry.

Graph drawing is placed on firmer metric ground: a layout is not only visually
pleasing if it preserves bars/edges; it becomes a representation of the graph's
geodesic metric.

## What Still Needs To Be Added

### A. Decide Manuscript Asset Structure

Use the SIMODS manuscript directory:

```text
/Users/pgajer/current_projects/grip/papers/data_to_geodesic_embedding_paper/
```

Recommended substructure:

```text
data_to_geodesic_embedding_paper/
  simods_data_geodesic_mds.tex
  references.bib
  figures/
  tables/
  scripts/
  notes/
```

The manuscript should link to generated assets from both `geodesicMDS` and
`geodesic_data_geometry`, but publication figures should eventually be copied
or regenerated into the manuscript-local `figures/` directory for reproducible
builds.

### B. Consolidate Graph Reconstruction Results

Need a concise result set from `geodesic_data_geometry` showing:

- which graph families were tested;
- which parameter grids were used;
- which examples were used;
- oracle surface-target performance;
- final recommendation that adaptive-radius and cKNN are the leading graph
  families for the next paper experiments;
- why fixed-radius and mKNN are not central methods in the paper.

The paper should not show every exploratory report.  It should distill them
into a small number of decisive tables/figures.

### C. Implement And Test Non-Oracle Graph Parameter Selection Across Families

The existing gflow k-selection machinery developed for iKNN should be audited
and generalized so it applies uniformly to:

- sKNN / union kNN if retained;
- adaptive-radius;
- cKNN;
- Delaunay reference graphs where relevant;
- possibly PHATE-style graph support as a comparison.

Selection criteria to implement or standardize:

1. GCV of graph smoothing on top variable features.
2. Jensen--Shannon distance between degree distributions of consecutive graphs.
3. Graph edit distance / edge symmetric difference between consecutive graphs.
4. Stability of selected graph statistics across parameter sweeps.
5. Connectivity/MST repair burden.
6. Optional layout-stability criteria.

The GCV path needs a clearly defined graph smoothing operator.  Candidate
operators:

- graph Laplacian Tikhonov smoothing:
  `min_f ||y - f||^2 + lambda f^T L f`;
- heat-kernel smoothing;
- Markov diffusion smoothing;
- local averaging over graph neighborhoods.

The manuscript should describe one primary operator and perhaps list others as
variants.

### D. Show Oracle Recovery Of Non-Oracle Choices

On quadratic surfaces, compare selected parameters to oracle-best parameters.
This is the bridge between controlled geometry and real-data usability.

Questions:

- Does GCV select parameters near oracle-optimal graph geodesic distortion?
- Does JS-degree stabilization select similar parameters?
- Are adaptive-radius and cKNN choices stable across sample size and curvature?
- Are selected parameters robust to noise?

### E. Define Final edge-KK Stopping Criteria

Current experiments use fixed budgets or target edge-rRMSE in some scripts.
The paper needs one coherent stopping rule.

Candidate rule:

```text
Stop when edge_rRMSE <= tau_edge,
or median relative improvement over the last W accepted iterations <= tau_improve,
or max_iter is reached.
```

Recommended diagnostics to store:

- edge rRMSE;
- q95 absolute relative edge residual;
- maximum edge stretch;
- energy;
- gradient norm;
- accepted step size;
- number of iterations;
- elapsed time.

Need experiments showing that the stopping rule behaves sensibly as `n`
increases.

### F. Decide Whether Repulsive Unfolding Belongs In This Paper

The current visual examination suggests that

```text
weighted GRIP -> edge-KK
metric MDS -> edge-KK
```

already give satisfactorily spread embeddings on two-dimensional quadratic
surfaces.  Repulsive unfolding may therefore be too much for the main paper
unless it solves a clear failure mode.

Recommendation:

- Keep repulsive unfolding as a discussion/future-method section unless a
  simple experiment shows a decisive benefit.
- Focus the paper on data graph reconstruction plus edge-KK GMDS.

### G. Real Data Case Studies

Choose one or two real datasets with interpretable geometry.

Candidates:

1. Cell-cycle single-cell data:
   likely strongest if the expected geometry is circular/periodic.

2. Vaginal microbiome 16S Valencia 13k:
   attractive for biological interpretability and state/gradient structure.

3. VIRGO2 2600 metagenomic dataset:
   smaller and possibly easier to use in a manuscript-scale experiment.

4. American Gut Project:
   large and public, but likely messier; may be better as a scalability or
   exploratory example than the first biological headline.

For each real dataset:

- define preprocessing;
- define feature space/distance;
- run non-oracle graph selection;
- compute GMDS layout;
- compare to metric MDS, PHATE/UMAP if appropriate, and graph-only GRIP;
- quantify graph/layout diagnostics;
- interpret known metadata gradients or clusters.

### H. Clarify Software Ownership

Current code spans:

- `gflow`: graph construction, geodesic data geometry utilities, graph
  selection machinery;
- `grip`: graph drawing and GMDS/edge-KK layout functions;
- `geodesicMDS`: manuscript, experiments, reports;
- `geodesic_data_geometry`: data-to-graph experiment reports and plans;
- `gmdsui`: visual frame/suite inspection.

For the paper this is acceptable, but the reproducibility story must be clean.
The manuscript should include a software map saying which package owns which
functionality.

### I. Decide Submission Target And Paper Length

Target: SIAM Journal on Mathematics of Data Science (SIMODS).

Implication:

- emphasize mathematical formulation and reproducible experiments;
- avoid presenting the paper as only a visualization/graph drawing method;
- keep biological examples as demonstrations, not the only source of evidence;
- put broad exploratory plots in supplements.

## Immediate Next Actions

1. Populate the SIMODS manuscript directory under
   `geodesicMDS/data_to_geodesic_embedding_paper/`.
2. Create a first LaTeX skeleton using the outline above.
3. Ask the `geodesic_data_geometry` orchestrator to produce a distilled graph
   reconstruction result package for the paper.
4. Audit gflow's existing k/parameter-selection methods and write a short
   inventory of what can be reused for adaptive-radius and cKNN.
5. Select the first real dataset case study, probably cell cycle or VIRGO2.
6. Define the canonical edge-KK stopping rule and retrofit the main scripts to
   record it.

## Current Position

The work is publishable in direction, but the paper is not yet submission-ready.
The strongest version is not "a new graph layout algorithm" and not "a graph
construction benchmark."  It is a unified mathematical data-science paper:

```text
recover data geodesic geometry by graph construction,
then represent that graph geometry by geodesic MDS.
```

That is a coherent SIMODS story.
