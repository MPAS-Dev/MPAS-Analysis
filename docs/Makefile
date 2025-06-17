# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = MPAS-Analysis
SOURCEDIR     = .
BUILDDIR      = _build

# Build into a versioned subdirectory
versioned-html:
	@echo "Building version: $(DOCS_VERSION)"
	$(SPHINXBUILD) -b html "$(SOURCEDIR)" "$(BUILDDIR)/html/$(DOCS_VERSION)"
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html/$(DOCS_VERSION)."
	@echo "Setting up shared version switcher for local preview..."
	mkdir -p _build/html/shared
	cp shared/version-switcher.js _build/html/shared/version-switcher.js
	python generate_versions_json.py --local

# Override html target to include local setup
html:
	$(SPHINXBUILD) -b html "$(SOURCEDIR)" "$(BUILDDIR)/html"
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."


# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

clean:
	rm -rf users_guide/*obs_table.rst developers_guide/generated users_guide/obs
	@$(SPHINXBUILD) -M clean "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

clean-versioned-html:
	rm -rf $(BUILDDIR)/html/*
	@echo "Cleaned versioned HTML builds."

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
