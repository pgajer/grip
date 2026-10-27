Vignette and website-article snapshot — 23 September 2026

Source repository commit: 07dc145c89ff31b399a7e4c0b8e1b86f391fc08b
Source directory: grip/vignettes/

This snapshot precedes the proposed documentation restructuring. It contains
all eight package vignette sources, all five website-only article sources,
the bibliography, and existing local HTML/image assets under that directory.
The live sources have not been modified or removed.

Archived Markdown/Rmd filenames have an underscore prefix so pkgdown does not
discover them as articles. Their contents are byte-for-byte copies. Original
paths, archived paths, and SHA-256 checksums are recorded in manifest.tsv.
The entire archives directory is excluded from R package builds.

Source copies, bibliography, and this record are version controlled. Existing
local generated HTML and GIF copies are retained locally and excluded from
Git. Those HTML files are copied as found, not rebuilt or claimed to be the
latest rendered previews. The main current previews remain under
~/.codex/private/grip/vignettes/.

This is a document-source snapshot, not a self-contained historical package:
Rmd sources still depend on package functions, installed scripts, and saved
results elsewhere in the repository. Use the recorded source commit to
recover those dependencies. Existing relative links in the copied HTML have
not been rewritten.

Checksum repair — 27 October 2026

The 30 September commit d213fe48de81953a96f8fb578bf8590aa7356ed1
("Add offline SuiteSparse graph catalogue and replace six-example dataset")
also edited _synthetic-graph-families.Rmd inside this archive. It replaced
zheng.graphs examples and related wording with suitesparse.graph calls,
without changing the original snapshot's manifest entry.

The original _synthetic-graph-families.Rmd has now been restored byte-for-byte
from the source commit recorded above; it matches its original SHA-256.
The later content is preserved byte-for-byte at:
  revisions/d213fe48/_synthetic-graph-families.Rmd
Its separate manifest row records the same original guide path and a distinct
archive path and checksum. Both versions are checked by check-documentation.py.
No archived content or original checksum was discarded. The revision is a
historical copy, not a maintained or independently validated runnable guide.

The preview location above records the September snapshot's environment.
Current maintained sources are under ~/current_projects/grip/; current local
previews are under output/vignette-previews/ in that checkout. Keep frozen
snapshot files unchanged; preserve any later variants separately with their
provenance and a new manifest entry.
