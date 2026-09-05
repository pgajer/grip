PYTHON ?= python3

.PHONY: run verify
run:
	OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 $(PYTHON) experiment.py

verify:
	$(PYTHON) verify.py
