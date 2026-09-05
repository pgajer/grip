ROOT := $(abspath ../../../../..)
REL := papers/grip-software-paper/reproducibility/experiments/mds-initializer-sensitivity
.PHONY: collect verify figures report citation-check artifact-check
collect:
	cd $(ROOT) && Rscript $(REL)/collect.R
verify:
	cd $(ROOT) && Rscript $(REL)/verify.R
figures:
	cd $(ROOT) && python3 $(REL)/plot.py
report: figures
	cd $(ROOT) && python3 $(REL)/build-report.py
citation-check:
	cd $(ROOT) && python3 $(REL)/check-citations.py

artifact-check:
	cd $(ROOT) && python3 $(REL)/check-artifacts.py
