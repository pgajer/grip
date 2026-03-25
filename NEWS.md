# grip 0.1.0

* Initial CRAN submission.
* Implements the GRIP multiscale force-directed graph layout algorithm in 2D
  and 3D via C++ (Rcpp).
* Provides edge-list generators for common graph families: paths, cycles,
  meshes, cylinders, tori, k-ary trees, Sierpinski triangles, tetrahedra,
  and carpets.
* Layout quality scoring with `grip.score.layout()` and multi-candidate
  comparison with `grip.compare.layouts()`.
* Interactive Shiny-based layout explorer (`gripui`) for browsing, filtering,
  and visualizing layout catalogs.
