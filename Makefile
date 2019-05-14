# Make a local package for import in api/app.py and Sphinx documentation.
# You can set these variables from the command line.
VERSION := 0.0.2

SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = PID_control
SOURCEDIR     = .
BUILDDIR      = _build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

all:
	flake8 PID_control.py
	sudo python3 setup.py sdist

docs:
	make latexpdf
	cp _build/latex/PID_control.pdf docs/PID_control_docs.pdf

clean:
	sudo rm -rf _build build dist PID_pendulum.egg-info

.PHONY: help all docs clean Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
