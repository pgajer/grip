# Radius study: results

[Full PDF](../report.pdf) · [Protocol and reproduction](../README.md)

The graph regime materially changes the geometric answer. At $r=64$, $k=32$, the fixed-scale stress-MDS-to-edge-KK branch gives the following base-disk medians across independent samples. Coordinate error uses similarity alignment; path/reference error compares edge-calibrated retained paths with smooth geodesics.

| Surface | Graph regime | s₂/s₁ | s₃/s₂ | Coordinate error (%) | Path/reference error (%) |
|---|---|---:|---:|---:|---:|
| Paraboloid | geodesic | 0.038 | 0.846 | 1.46 | 0.34 |
| Paraboloid | ambient | 0.028 | 0.915 | 0.88 | 0.38 |
| Saddle | geodesic | 0.901 | 0.896 | 80.38 | 0.36 |
| Saddle | ambient | 0.019 | 0.917 | 0.90 | 37.00 |

Across 1176 primary $n=240$ graphs, profiled-scale edge-KK lowers fixed-path error in 1016 stress-initialized comparisons and increases it in 76; fixed-scale edge-KK lowers it in 1075 and increases it in 17. Changes smaller than $10^{-8}$ in absolute relative-error units are counted as ties. Good path fidelity is therefore neither a universal refinement improvement nor evidence of surface recovery.

The pilot also exposed a scale issue: the profiled, unnormalized edge objective can decrease by uniform contraction. Fixed-scale controls prevent this degeneracy and are essential when interpreting absolute loss, stopping values, and spatial extent.

![Graph fidelity](../figures/graph-reference.png)

The experiment contains 1216 primary graphs, 16416 scored primary candidates, and 3648 MDS starts. 74 graphs require MST augmentation; the largest retained-route/strict-distance discrepancy is 6.67e-08 in normalized units. Independent edge/path calculations agree with grip to 1.25e-14; explicit route sums agree to 8.88e-16.

All 88 numerical input matrices pass symmetry, all-pair chord lower bounds, explicit surface-path upper bounds, and exhaustive sampled-vertex triangle checks at tolerance $10^{-7}$ times RMS distance. The largest normalized triangle violation is 8.49e-09. Independent high-radius validation covers 256 actual pairs in 16 cases, including short/far and random pairs. The largest sampled relative discrepancy is 5.16e-08. Saddle checks use tighter shooting and independent SciPy collocation/integration; paraboloid checks use independent quadrature. This sampled maximum is not an all-pair accuracy bound. A separate reconstruction of graph distances and diagnostics from 16928 saved coordinate arrays agrees to 2.59e-12 after scaling dimensional discrepancies. Primary profiled fits reach edge-calibration factors as small as 4.07e-10, whereas the fixed-scale controls range from 0.99339 to 1.02624. The contraction occurs at intermediate radii too, so it is not an asymptotic property of the surface.

25/3648 MDS starts reach the iteration cap; the largest between-start difference in target-normalized RMSE is 1.530 percentage points. Iteration limits and small stress changes do not certify global optimality. The additional controls contain 512 scored candidates on 32 graphs. Under profiled scale, the largest $s_3/s_2$ change is 0.0003 for tiny perturbations, 0.0124 for extra steps, and 0.7807 for independent random starts. Under fixed scale, the largest $s_3/s_2$ change is 0.0003 for tiny perturbations, 0.0263 for extra steps, and 0.8208 for independent random starts. These distinguish local stability around MDS starts from dependence on the choice of a different starting configuration.

The results support presenting edge-KK as a graph-edge refinement with separately assessed path and geometric fidelity. They do not support a general claim that it reverses MDS flattening or recovers the generating manifold. The graph/reference discrepancy, neighborhood locality, and scale policy must be visible alongside the final graph-to-layout score. The complete ambient control is exact Euclidean realization, not evidence of intrinsic-geodesic recovery.

All factor combinations are in [scores.csv](scores.csv); paired initializer and sample-size comparisons are in [initializer-comparison.csv](initializer-comparison.csv) and [size-comparison.csv](size-comparison.csv).
