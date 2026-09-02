REPORT := S1-weighted-grip-complexity
MOBIUS_REPORT := S2-mobius-comparison
MOBIUS_SCRIPT := ../reproducibility/scripts/mobius-dimension-comparison.R
BUILD_DIR ?= ../build/supplement
CHECKER ?= $(if $(wildcard check-citation-verification.py),check-citation-verification.py,../../../tools/reports/check-citation-verification.py)

.PHONY: all pdf citation-check argument-check
all: pdf citation-check argument-check

argument-check:
	python3 check-complexity-arguments.py

$(BUILD_DIR)/mobius-dimension-comparison.pdf: $(MOBIUS_SCRIPT)
	@mkdir -p "$(BUILD_DIR)"
	Rscript "$(MOBIUS_SCRIPT)" "$(BUILD_DIR)"

pdf: $(BUILD_DIR)/mobius-dimension-comparison.pdf
	@mkdir -p "$(BUILD_DIR)"
	latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir="$(BUILD_DIR)" "$(REPORT).tex"
	latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir="$(BUILD_DIR)" "$(MOBIUS_REPORT).tex"

citation-check:
	python3 "$(CHECKER)" --tex "$(REPORT).tex" --bib "$(REPORT).bib" --html citation_verification.html --log "$(BUILD_DIR)/$(REPORT).log"
	python3 "$(CHECKER)" --tex "$(MOBIUS_REPORT).tex" --bib "$(MOBIUS_REPORT).bib" --html S2-citation_verification.html --log "$(BUILD_DIR)/$(MOBIUS_REPORT).log"
