PYTHON ?= python3
.PHONY: figures pdf verify citation-check review-bundle
figures:
	$(PYTHON) scripts/check_coordinate_components.py --refresh
	$(PYTHON) scripts/check_focused_numbers.py
	$(PYTHON) scripts/build_focused_figures.py
pdf:
	zsh scripts/build_geodesic_mds_pdf.sh
verify:
	$(PYTHON) scripts/verify_focused_evidence.py
	$(PYTHON) scripts/check_focused_mathematics.py
	$(PYTHON) scripts/check_coordinate_components.py
	$(PYTHON) scripts/check_focused_numbers.py
	$(PYTHON) scripts/check_focused_citations.py
citation-check:
	$(PYTHON) scripts/check_focused_citations.py

review-bundle: pdf verify
	$(PYTHON) scripts/build_review_bundle.py
