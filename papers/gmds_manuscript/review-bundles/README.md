# Immutable author-review packages

Dated packages preserve versions circulated for review. Their TeX lives inside source archives; never edit an extracted review copy as a new master. `make review-bundle` creates a new date/time/commit-named directory and refuses to overwrite an existing one. Build intermediates belong only in `../build/`.
