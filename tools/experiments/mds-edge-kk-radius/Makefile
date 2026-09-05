ROOT := ../../..
HERE := tools/experiments/mds-edge-kk-radius
PYTHON ?= python3
.PHONY: verify report collect citation-check refit
verify:
	cd $(ROOT) && $(PYTHON) $(HERE)/restore-inputs.py
	cd $(ROOT) && $(PYTHON) $(HERE)/verify.py
	cd $(ROOT) && Rscript $(HERE)/verify-coordinates.R
	$(PYTHON) check-citations.py
	$(PYTHON) check-artifacts.py
report:
	cd $(ROOT) && $(PYTHON) $(HERE)/plot.py
	$(PYTHON) build-report.py
	$(PYTHON) check-citations.py
	$(PYTHON) freeze-artifacts.py
collect:
	cd $(ROOT) && Rscript $(HERE)/collect.R final
	cd $(ROOT) && Rscript $(HERE)/diagnostics.R
	cd $(ROOT) && $(PYTHON) $(HERE)/comparisons.py
	cd $(ROOT) && Rscript $(HERE)/verify-objective.R
	cp $(ROOT)/output/mds-edge-kk-radius/objective-validation.csv summary/objective-validation.csv
	cd $(ROOT) && $(PYTHON) $(HERE)/verify.py
citation-check:
	$(PYTHON) check-citations.py
refit:
	cd $(ROOT) && $(PYTHON) $(HERE)/restore-inputs.py
	cd $(ROOT) && $(PYTHON) $(HERE)/run-grid.py main --workers 6
	cd $(ROOT) && $(PYTHON) $(HERE)/run-grid.py size --workers 2
	cd $(ROOT) && $(PYTHON) $(HERE)/run-grid.py extra --workers 4
