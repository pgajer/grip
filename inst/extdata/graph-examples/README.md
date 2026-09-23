# Graph examples gallery

`layouts.rds` records 13 graphs and their 3D metric-MDS and edge-KK layouts:
six generated graphs, karate club, and six SuiteSparse examples. Each input
retains its edge lengths and, when available, its generating coordinates.
The PNG files show ivue views of the same saved coordinates.

Every graph uses three random MDS starts, seed 2026, and 80 SGD passes per start;
the lowest-stress start is selected by metric.mds(). Edge-KK starts from that
layout with uniform stiffness, fixed target lengths (scale_mode = "identity"),
and up to 100 iterations. These are illustration choices, not tuned settings.
Both graph-distance and edge-length errors are independently recomputed from
returned coordinates. Timings exclude shared graph preparation. Alignment
permits translation and rotation/reflection only, with common display limits.

The saved protocol identifies software and fitting sources. Start summaries,
warnings and unsuccessful outcomes are retained; a stopped graph does not
prevent later graphs from running. The recorded source commit precedes these
new scripts; their checksums identify the actual new source files used.

The final saddle and paraboloid use diagonal as well as horizontal/vertical
connections. An earlier exploratory run used orthogonal connections, which
produce identical edge lengths for these two lifts. That earlier run remains
separate and is not mixed into this gallery.

Reproduction scripts are documented in tools/README.md. Raw logs, intermediate
attempts and browser checks are development records, not installed data.
SuiteSparse sources and CC BY 4.0 attribution are in ../zheng-graphs/PROVENANCE.md
and in each corresponding graph's graph_info field. The original matrices and
Les Miserables character labels are archived in that directory.
