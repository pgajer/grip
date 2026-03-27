.PHONY: clean distclean attrs document build build-verbose build-log check check-clean check-fast check-examples check-dir install readme readme-assets readme-html readme-render paper-pdf paper-html paper-all rchk winbuilder-release winbuilder-devel

PKGNAME := grip
VERSION := $(shell awk '/^Version:/ { print $$2 }' DESCRIPTION)
TARBALL := $(PKGNAME)_$(VERSION).tar.gz
LOGDIR := .claude

clean:
	find src -type f \( -name "*.o" -o -name "*.so" -o -name "*.dll" \) -delete
	rm -rf $(PKGNAME).Rcheck ..Rcheck *.Rcheck
	rm -f $(TARBALL) ../$(TARBALL)
	rm -f $(LOGDIR)/*.log

distclean: clean
	rm -rf $(LOGDIR)

attrs:
	@mkdir -p $(LOGDIR)
	@echo "Running Rcpp::compileAttributes()..."
	@R -q -e "Rcpp::compileAttributes()" > $(LOGDIR)/$(PKGNAME)_rcppattrs.log 2>&1
	@echo "RcppExports regenerated (log: $(LOGDIR)/$(PKGNAME)_rcppattrs.log)"

document: attrs
	@mkdir -p $(LOGDIR)
	@echo "Running devtools::document()..."
	@R -q -e "devtools::document()" > $(LOGDIR)/$(PKGNAME)_document.log 2>&1
	@echo "Documentation generated (log: $(LOGDIR)/$(PKGNAME)_document.log)"

build: clean document
	@mkdir -p $(LOGDIR)
	@echo "Building package..."
	@R CMD build . > $(LOGDIR)/$(PKGNAME)_build.log 2>&1
	@echo "Package built successfully (log: $(LOGDIR)/$(PKGNAME)_build.log)"

build-verbose: clean document
	R CMD build .

build-log: clean document
	@mkdir -p $(LOGDIR)
	R CMD build . > $(LOGDIR)/$(PKGNAME)_build.log 2>&1
	@echo "Build output saved to $(LOGDIR)/$(PKGNAME)_build.log"

check: build
	R CMD check --as-cran $(TARBALL)

check-clean: build
	env R_MAKEVARS_USER=/dev/null R CMD check --as-cran $(TARBALL)

check-fast: build
	R CMD check --as-cran --no-examples --no-tests --no-manual $(TARBALL)

check-examples: build
	R CMD check --as-cran --examples $(TARBALL)

check-dir:
	R CMD check --as-cran .

install: build
	R CMD INSTALL $(TARBALL)

readme-assets:
	Rscript tools/generate-readme-assets.R

readme-render:
	Rscript tools/render-readme.R

readme: readme-assets readme-render

readme-html: readme-assets
	Rscript tools/render-readme.R --html

paper-pdf:
	Rscript tools/render-paper.R

paper-html:
	Rscript tools/render-paper.R --html-only --html

paper-all:
	Rscript tools/render-paper.R --all

rchk:
	@tools/check_rchk.sh

winbuilder-release: build
	Rscript tools/check-win-builder.R release

winbuilder-devel: build
	Rscript tools/check-win-builder.R devel
