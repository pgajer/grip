# Graph examples from the SuiteSparse Matrix Collection

The six matrices in this directory are unchanged Matrix Market text files,
compressed with gzip, downloaded on 2026-09-23 from the University of Florida
mirror of the SuiteSparse Matrix Collection. The Les Miserables vertex-name
file is also unchanged. Each matrix retains its original metadata header.

| Source | Original authors | Editors |
|---|---|---|
| [HB/dwt_66](https://sparse.tamu.edu/HB/dwt_66) | G. Everstine, D. Taylor | I. Duff, R. Grimes, J. Lewis |
| [Newman/lesmis](https://sparse.tamu.edu/Newman/lesmis) | D. Knuth | M. Newman |
| [HB/dwt_307](https://sparse.tamu.edu/HB/dwt_307) | G. Everstine, D. Taylor | I. Duff, R. Grimes, J. Lewis |
| [HB/494_bus](https://sparse.tamu.edu/HB/494_bus) | D. Tylavsky | I. Duff, R. Grimes, J. Lewis |
| [HB/dwt_1005](https://sparse.tamu.edu/HB/dwt_1005) | G. Everstine, D. Taylor | I. Duff, R. Grimes, J. Lewis |
| [HB/1138_bus](https://sparse.tamu.edu/HB/1138_bus) | D. Tylavsky | I. Duff, R. Grimes, J. Lewis |

Source archives: `https://www.cise.ufl.edu/research/sparse/MM/GROUP/NAME.tar.gz`,
using the group/name in the table. The current Texas A&M HTTPS file server was
unreachable during this import; the collection's original Florida mirror was used.

## License and attribution

SuiteSparse distributes the matrices under
[Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/),
as stated in the [collection's license policy](https://sparse.tamu.edu/about).
This license applies to these matrices, their metadata and the derived
`zheng.graphs` dataset, separately from grip's software license.

Please credit the original authors above and cite Davis, T. A. and Hu, Y.
(2011), *The University of Florida Sparse Matrix Collection*, ACM Transactions
on Mathematical Software 38(1), Article 1,
https://doi.org/10.1145/2049662.2049663. For Les Miserables, also cite Knuth,
D. E. (1993), *The Stanford GraphBase: A Platform for Combinatorial Computing*,
Addison-Wesley. The complete source-specific notes remain in each matrix header.

The examples were selected from Zheng, J. X., Pawar, S. and Goodman, D. F. M.
(2018), *Graph Drawing by Stochastic Gradient Descent*,
https://doi.org/10.1109/TVCG.2018.2859997.

## Derived graphs

`data-raw/zheng_graphs.R` rebuilds `data/zheng.graphs.rda` offline from these
files. Derived names end in `_unweighted_graph`. The conversion removes diagonal
entries, discards explicit zeros, combines symmetric/duplicate edges, and assigns
unit traversal lengths. All vertices and original row indices are retained.
All six graphs are connected; no component was dropped or repaired. Original
matrix values are preserved in the archived files, not treated as distances.
The source names and per-file checksums are retained in each graph's metadata.
