REPORT := weighted-grip-complexity
BUILD_DIR ?= ../../../../output/pdf/weighted-grip-complexity
CHECKER ?= ../../../../tools/reports/check-citation-verification.py

.PHONY: all pdf citation-check
all: pdf citation-check

pdf:
	@mkdir -p "$(BUILD_DIR)"
	@printf '\\renewcommand{\\reportbuilddatetime}{%s}\n' "$$(TZ=America/New_York date '+%Y-%m-%d %H:%M:%S %Z')" > "$(BUILD_DIR)/$(REPORT)-build-info.tex"
	TEXINPUTS="$(abspath $(BUILD_DIR)):" latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir="$(BUILD_DIR)" "$(REPORT).tex"

citation-check:
	python3 "$(CHECKER)" --tex "$(REPORT).tex" --bib "$(REPORT).bib" --html citation_verification.html --log "$(BUILD_DIR)/$(REPORT).log"
