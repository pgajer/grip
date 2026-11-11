.PHONY: report citation-check
report:
	@test -n "$(BUILD_DIR)" || (echo 'Set BUILD_DIR to an output directory outside the source tree.'; exit 1)
	python3 build_report.py --pdf --output-dir "$(BUILD_DIR)"
citation-check:
	@test -n "$(CITATION_EVIDENCE)" || (echo 'Set CITATION_EVIDENCE to the citation-verification HTML.'; exit 1)
	@test -n "$(BUILD_DIR)" || (echo 'Set BUILD_DIR to the rendered note directory.'; exit 1)
	python3 check_citations.py --tex methods.tex --bib references.bib --html "$(CITATION_EVIDENCE)" --log "$(BUILD_DIR)/methods.log"
