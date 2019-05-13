# Make a local package for import in api/app.py and Sphinx documentation.
# You can set these variables from the command line.
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
	sudo -H pip3 install -e .
	sudo python3 setup.py bdist_wheel
	sudo -H pip3 install dist/PID_pendulum-0.0.1-py3-none-any.whl

# don't use this target
docs:
	make latexpdf
	cp _build/latex/PID_control.pdf docs/

.PHONY: help all docs Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
