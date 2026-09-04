REPORT := S1-weighted-grip-complexity
MOBIUS_REPORT := S2-mobius-comparison
MOBIUS_SCRIPT := ../reproducibility/scripts/mobius-dimension-comparison.R
CONTROLLED_REPORT := S3-controlled-examples
EDGE_WORKFLOW_SCRIPT := ../reproducibility/scripts/edge-kk-workflow.R
BUILD_DIR ?= ../build/supplement
CHECKER ?= $(if $(wildcard check-citation-verification.py),check-citation-verification.py,../../../tools/reports/check-citation-verification.py)

.PHONY: all pdf citation-check argument-check
all: pdf citation-check argument-check

argument-check:
	python3 check-complexity-arguments.py

$(BUILD_DIR)/$(CONTROLLED_REPORT).pdf: $(EDGE_WORKFLOW_SCRIPT)

$(BUILD_DIR)/mobius-dimension-comparison.pdf: $(MOBIUS_SCRIPT)
	@mkdir -p "$(BUILD_DIR)"
	Rscript "$(MOBIUS_SCRIPT)" "$(BUILD_DIR)"

$(BUILD_DIR)/$(CONTROLLED_REPORT).pdf: $(CONTROLLED_REPORT).Rmd $(CONTROLLED_REPORT).bib render-S3.R ../reproducibility/precomputed/two-fidelity-saddle.rds ../reproducibility/scripts/plot-two-fidelity.R ../reproducibility/precomputed/vs_alternatives/benchmark_results.rds ../reproducibility/scripts/plot-weighted-saddle.R ../reproducibility/precomputed/saddle-reference-diagnostics.rds ../reproducibility/precomputed/single-saddle-reference.csv ../reproducibility/experiments/two-fidelity-pilot/plot-saddle-reference.R $(wildcard ../reproducibility/figures/saddle/*.png)
	Rscript render-S3.R "$(BUILD_DIR)"

pdf: $(BUILD_DIR)/mobius-dimension-comparison.pdf $(BUILD_DIR)/$(CONTROLLED_REPORT).pdf
	@mkdir -p "$(BUILD_DIR)"
	latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir="$(BUILD_DIR)" "$(REPORT).tex"
	latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir="$(BUILD_DIR)" "$(MOBIUS_REPORT).tex"

citation-check:
	python3 "$(CHECKER)" --tex "$(REPORT).tex" --bib "$(REPORT).bib" --html citation_verification.html --log "$(BUILD_DIR)/$(REPORT).log"
	python3 "$(CHECKER)" --tex "$(MOBIUS_REPORT).tex" --bib "$(MOBIUS_REPORT).bib" --html S2-citation_verification.html --log "$(BUILD_DIR)/$(MOBIUS_REPORT).log"
	python3 "$(CHECKER)" --tex "$(BUILD_DIR)/$(CONTROLLED_REPORT).tex" --bib "$(CONTROLLED_REPORT).bib" --html S3-citation_verification.html
