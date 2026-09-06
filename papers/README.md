# Papers associated with grip

Each paper has a separate scope and one editable source home.

| Directory | Purpose and status | Source and build |
| --- | --- | --- |
| [grip-software-paper](grip-software-paper/) | R Journal software article | See its README; `make paper-pdf` from the grip root |
| [gmds_manuscript](gmds_manuscript/) | Methods paper on path fidelity, scale and geometric recovery; author-review draft | `geodesic_mds.tex` at that directory root; `make gmds-paper-pdf` |
| [data_to_geodesic_embedding_paper](data_to_geodesic_embedding_paper/) | Deferred broader data-to-graph and embedding research agenda, with scientific planning and literature assets | No principal manuscript exists yet; see its README |

Edit manuscript masters, never dated review copies. Disposable compilation files
belong in each paper's `build/`; immutable circulation packages belong in its
`review-bundles/`. Evidence and reproducibility inputs remain with their papers.
Internal editorial working notes live outside the repository.

The package's `.Rbuildignore` excludes this entire directory from R source packages.
