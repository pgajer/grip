# Native SGD backend for metric MDS

## Scope and compatibility

`metric.mds(backend = "sgd")` minimizes the existing uniform all-pairs raw
distance stress. `backend = "sgd"` is the default; explicitly select
`backend = "smacof"` to reproduce a previous SMACOF analysis. Both new arguments,
`backend` and `sgd_control`, are appended to the public signature. Inputs,
distance normalization, initialization, output units, diagnostic definitions,
and multiple-start selection retain the existing contract. There is no Python
runtime dependency, automatic fallback, weighted objective, sparse surrogate,
or lazy graph-distance calculation. Refinement methods retain their default
classical initializer; explicitly
selecting `init = "metric_mds"` now follows SGD. Supply coordinates explicitly
to preserve a SMACOF initializer.

## Frozen references and provenance

* Zheng et al., *Graph Drawing by Stochastic Gradient Descent* (2018):
  <https://doi.org/10.1109/TVCG.2018.2859997>.
  `jxz12/s_gd2` revision `52ab0a5bee183b45eed60061cabe8e63f006e8a0`,
  `cpp/s_gd2/layout.cpp` SHA256
  `37fd3f7d3452eb58017adc8ed555c5b48b23945d7d8b0fa95e22dec87386cfa0`.
* Hangan et al., *Bridging Graph Drawing and Dimensionality Reduction with
  Stochastic Stress Optimization* (2026): <https://arxiv.org/abs/2605.00641>.
  `DanielHangan01/scikit-learn` revision
  `7c4301d490c9ec42400ac90ffd69308bfbf4863e`.
  `_sgd_mds.py` SHA256
  `077c1bf33f1a924605144f19efffac26419e21eb5589c755df913767a234df65`;
  `_sgd_mds_cython.pyx` SHA256
  `e14d6a70f0895200f5706bdc46ff857fcad5b004918fd536ac27eedcb55b4b71`.

The new C++ implementation adapts the pair update and the exponential/hybrid
schedules. It does not vendor the Python estimator, upstream PRNG, graph-distance
code, or sparse methods. MIT and BSD notices are retained in `inst/COPYRIGHTS`.
Algorithm identity: `grip-sgd-mds-v1`. Matching upstream random seeds does not
imply identical trajectories.

## Exact optimization rules

Distances are divided by their root mean square as in the SMACOF wrapper.
Pairs are enumerated `(0,1), (0,2), ..., (n-2,n-1)`, which matches R's condensed
`dist` order. Every epoch resets that order and applies descending Fisher-Yates
shuffling with an unbiased rejection-sampled integer from `std::mt19937_64`.
The seed for each start is drawn from the caller's R RNG (recorded in results).
With a supplied public seed, the R RNG state is restored on return, including
errors. Native shuffling consumes no further R random numbers.

For a pair with separation vector `v = z_i-z_j`, distance `r`, and target `d`,
set `mu = min(eta, 1)` and update both endpoints symmetrically by
`mu * (r-d) / 2 * (v/r)`. Compute the unit vector before multiplying to avoid
overflow in `(r-d)/r`. A coincident pair with positive target uses a seeded
unit vector obtained by normalizing independent uniform coordinates in [-1,1).
This exceptional direction is not rotation-equivariant. A coincident pair with
zero target does nothing; zero targets are not replaced by a positive floor.
This explicitly differs from the upstream Cython kernel's skipped collisions.
Without collisions, updates preserve the initial affine span; deficient starts
are not silently perturbed. Wholly collapsed public starts remain invalid.

Schedules use zero-based epoch t and planned epoch count T. Exponential:
`eta(t) = eta0 * (eta_final/eta0)^(t/T)`. Hybrid uses
`tau = floor(switch_ratio*T)` and `eta_mid = eta0/10`: exponential decay from
eta0 toward eta_mid for t < tau, followed by
`eta_mid / (1 + b*(t-tau))`, where
`b = (eta_mid/eta_final - 1)/(T-tau)` when eta_mid > eta_final, and zero
otherwise. If tau is zero, eta_mid is eta0. These are the pinned reference's
indexing conventions: for exponential decay, and for an active hybrid harmonic
segment with eta_mid > eta_final, eta_final is the planned boundary value at
t=T, not necessarily the last applied value. If eta_mid <= eta_final, the
harmonic segment stays constant at eta_mid instead. If switch_ratio=1 there
is no harmonic segment and eta_final does not affect the schedule. Only
positive, nonincreasing schedules are
accepted (eta_final <= eta0).

Provisional calibration choices: hybrid schedule, eta0=0.5, eta_final=0.01,
switch_ratio=0.4, checkpoint_every=1, maximum native workspace 256 MiB.
The public `max_iter` retains its existing default of 1000; explicitly request
30 for an initial short SGD trial. An epoch is not a SMACOF iteration.
`eps` belongs to SMACOF; an explicitly supplied `eps` with SGD is rejected.
SGD performs its prescribed schedule and reports `iteration_limit`, never a
convergence certificate. No automatic plateau criterion is implemented.

At epoch zero, each checkpoint, and the terminal epoch, independently recompute
raw and optimally scale-profiled stress by summing all pair residuals. This
does not rescale the ongoing search. Retain the configuration with smallest
profiled stress; exact ties retain the earlier configuration. The R wrapper
independently profiles/scales and scores the selected configuration again.
Trace stress is in normalized distance units, with the normalization recorded.
The terminal trace remains available even when an earlier checkpoint wins.

Pair counts use checked size_t arithmetic. Native workspace estimates cover
the pair array, shuffled indices, two coordinate buffers, and checkpoint
storage; the R distance matrix, R copies, and returned trace incur additional
memory. Reject excessive workspace before allocating those buffers. Check user
interrupts during pair construction, updates, and scoring. Numerical failures
fail the affected start; the wrapper records them and continues other starts.
User interrupts propagate instead of being converted into failed starts.

## Validation and evidence gates

1. Compare fixed-order pair updates, full epochs, schedule values, and checkpoint
   selection against a small independent R implementation.
2. Verify exact Euclidean targets, imperfect starts, duplicate observations,
   rank-deficient starts, dimensionality, unit changes, RNG restoration, invalid
   inputs, absent SMACOF, multiple starts, and failure isolation.
3. Obtain an independent implementation audit before performance comparisons.
4. Freeze a bounded benchmark manifest before running it. Match targets and
   supplied starts, report achieved stress versus elapsed time and pair work,
   include scoring costs, and retain every failed or interrupted fit.
5. Run package checks and report platform coverage honestly. A local check is
   not evidence of checks on other operating systems. Cross-platform publication
   readiness must be established separately.

Frozen geometry reference fixtures and unrelated working-tree changes must not
be modified. Generated evidence belongs under ignored `output/`, not `dev/`.

## Default-backend decision

The initial implementation and pilot kept SMACOF as the default. Following the
completed pilot and user direction, the public default changes to SGD. The
native kernel, optimization objective, iteration budget, tuning controls and
stopping semantics are unchanged. At the predeclared 0.3-second threshold, SGD
had lower stress in 42 of 50 paired cases; SMACOF had lower stress in eight.
These were small Euclidean-cloud, weighted-graph and karate-club cases, with
two and three dimensions and five matched starts, evaluated retrospectively
over three separately restarted schedules. The comparison does not validate
geometry recovery, larger datasets, or the tuning defaults.
