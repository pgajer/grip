# Graph Hessian Literature Review Note

Built May 20, 2026.

This note records a short literature review answering the question:

> Were there attempts to define a graph version of the Hessian, understood as
> the matrix of partial second derivatives?

## Short Answer

Yes. There have been several attempts, but there is **no single canonical
"graph Hessian"** for an arbitrary abstract graph in the same way there is a
canonical graph Laplacian.

The main reason is that the Hessian is not just "second differences"; it is a
**matrix/tensor of second partial derivatives along coordinate directions**. A
plain graph gives neighborhoods and edge connectivity, but not coordinate axes,
tangent bases, or mixed partial directions. So most graph Hessian constructions
add extra structure: a point-cloud embedding, local tangent coordinates,
voxel/grid axes, a manifold assumption, or a graph-derived local coordinate
system.

## Most Relevant Attempts

### 1. Hessian Eigenmaps / Hessian LLE

Donoho and Grimes, 2003, *Hessian Eigenmaps*.

This is probably the cleanest early answer. They replace the Laplacian
Eigenmaps quadratic form with a Hessian-based quadratic form. For data points
sampled from a manifold, they estimate local tangent coordinates from
neighborhoods and define a local Hessian estimator. The global operator is not
usually stored as a per-node Hessian matrix, but as a sparse quadratic form
penalizing Hessian energy.

Key idea:

```math
H(f)=\int_M \| \operatorname{Hess} f(m)\|_F^2\,dm.
```

The paper explicitly says it substitutes a Hessian-based quadratic form for a
Laplacian-based one, and defines the Hessian through local tangent coordinates.

Source: [Donoho and Grimes PDF](https://www.math.ucdavis.edu/~saito/data/high-dimensions/HessianEigenmaps.pdf).

### 2. Hessian Energy for Semi-Supervised Regression

Kim, Steinke, and Hein, NeurIPS 2009, *Semi-supervised Regression using
Hessian energy*.

This is highly relevant if the target is graph/manifold regularization. They
argue graph-Laplacian regularization biases solutions toward constants and
propose a second-order Hessian energy that prefers functions varying linearly
along the manifold. They estimate the Hessian energy from samples, typically
via neighborhood graphs and local coordinates.

This is a graph/manifold Hessian regularizer more than a literal per-node
Hessian matrix.

Sources:

- [NeurIPS abstract](https://papers.neurips.cc/paper/3741-semi-supervised-regression-using-hessian-energy-with-an-application-to-semi-supervised-dimensionality-reduction)
- [PDF](https://papers.neurips.cc/paper_files/paper/2009/file/f4552671f8909587cf485ea990207f3b-Paper.pdf)

### 3. Discrete Hessian Eigenmaps

Ye and Zhi, 2015, *Discrete Hessian Eigenmaps method for dimensionality
reduction*.

This explicitly develops a discrete version of Hessian Eigenmaps and introduces
a "discrete Hessian operator." It is still in the manifold-learning /
point-cloud-neighborhood setting, not a universal abstract-graph Hessian.

Source: [ScienceDirect page](https://www.sciencedirect.com/science/article/pii/S0377042714004075).

### 4. Hessian Regularization Follow-On Work

There is a small family of machine-learning papers using Hessian
regularization, often for semi-supervised learning, image annotation, NMF, or
GCN-style models. Examples include:

- *Hessian regularization by patch alignment framework*
- *HesGCN: Hessian graph convolutional networks for semi-supervised classification*
- multiview Hessian regularization papers

These generally inherit the Hessian-energy idea rather than defining a clean
graph Hessian tensor from first principles.

Source example: [HesGCN ScienceDirect page](https://www.sciencedirect.com/science/article/pii/S0020025519310643).

### 5. Voxel / Mesh / Discrete Differential Geometry Operators

A more literal "matrix of partial second derivatives" appears when the graph
comes from a voxel complex, mesh, or embedded discretized domain. In that case
there are coordinate axes or local cells, so one can define gradient,
divergence, Laplacian, and Hessian-like operators algebraically.

A recent example is Nourian and Azadi, 2024, *Voxel graph operators*, which
derives algebraic differential operators for voxel connectivity graphs and
explicitly discusses a Hessian obtained from graph/voxel gradient operators.

Source: [ScienceDirect open-access page](https://www.sciencedirect.com/science/article/pii/S0965997824001297).

### 6. High-Order Graph Regularization / Hodge Laplacians

Zhou and Burges, 2008, *High-Order Regularization on Graphs*, develops
higher-order graph regularization using a discrete Laplace-de Rham operator.
This is related, but it is not exactly a Hessian matrix of second partial
derivatives. It is better viewed as Hodge/Laplace-type higher-order graph
calculus.

Source: [Microsoft Research page](https://www.microsoft.com/en-us/research/publication/high-order-regularization-graphs/).

### 7. Bethe Hessian

The "Bethe Hessian" is important in spectral clustering and community
detection, but it is **not** a Hessian of a scalar graph signal in the calculus
sense. It is a symmetric deformed Laplacian related to the Bethe free energy /
non-backtracking operator.

Source: [NeurIPS 2014 page](https://papers.neurips.cc/paper/5520-spectral-clustering-of-graphs-with-the-bethe-hessian).

## Bottom Line

There are real attempts, especially under names like:

- Hessian Eigenmaps
- Hessian LLE
- Hessian energy
- Hessian regularization
- discrete Hessian operator
- voxel graph Hessian / discrete differential operators

But for a **plain weighted graph**, there is no universally accepted analogue
of the full Hessian matrix. The graph Laplacian is the standard second-order
scalar operator, roughly like the trace of the Hessian. To get a full
Hessian-like object, one usually needs:

- local tangent coordinates;
- an embedding;
- edge directions;
- a voxel/grid/mesh structure;
- or a learned/local coordinate system.

For SIMODS-style graph smoothing or graph-selection work, the most relevant
branch is probably **Hessian energy regularization** rather than Bethe Hessian
or Hodge Laplacians. It would let one compare Laplacian smoothing, which favors
constants, against a second-order regularizer whose null space includes locally
linear or geodesic-linear functions.

