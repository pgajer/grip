REPORT := S1-weighted-grip-complexity
BUILD_DIR ?= ../build/supplement
CHECKER ?= $(if $(wildcard check-citation-verification.py),check-citation-verification.py,../../../tools/reports/check-citation-verification.py)

.PHONY: all pdf citation-check
all: pdf citation-check

pdf:
	@mkdir -p "$(BUILD_DIR)"
	latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir="$(BUILD_DIR)" "$(REPORT).tex"

citation-check:
	python3 "$(CHECKER)" --tex "$(REPORT).tex" --bib "$(REPORT).bib" --html citation_verification.html --log "$(BUILD_DIR)/$(REPORT).log"
