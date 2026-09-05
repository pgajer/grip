PYTHON ?= python3

.PHONY: run geodesic verify
run:
	OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 $(PYTHON) experiment.py

geodesic:
	OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 $(PYTHON) geodesic_experiment.py

verify:
	$(PYTHON) verify.py
