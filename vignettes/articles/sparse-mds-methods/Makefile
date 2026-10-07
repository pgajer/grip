.PHONY: report citation-check
report:
	python3 build_report.py --pdf
citation-check:
	python3 /Users/pgajer/.codex/notes/agent_instructions/reports/scripts/check_citation_verification.py --tex methods.tex --bib references.bib --html citation_verification.html --log methods.log
