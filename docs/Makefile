SPHINXOPTS    ?= -W
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

html:
	$(SPHINXBUILD) -b html $(SPHINXOPTS) $(SOURCEDIR) $(BUILDDIR)/html

clean:
	rm -rf $(BUILDDIR)

.PHONY: html clean
