## makefile for corrfitter library
#
# Created by G. Peter Lepage (Cornell University) on 2011-09-22.
# Copyright (c) 2011-2017 G. Peter Lepage.

PYTHON = python
PIP = $(PYTHON) -m pip
VERSION = `cd ..; python -c 'import corrfitter; print corrfitter.__version__'`

DOCFILES :=  $(shell ls doc/source/*.{rst,py,png} examples/*.png)
SRCFILES := $(shell ls setup.py src/*.py)

install-user:
	$(PIP) install . --user

#	$(PYTHON) setup.py install --user	--record files-corrfitter.$(PYTHON)

install install-sys:
	$(PIP) install .

#	$(PYTHON) setup.py install	--record files-corrfitter.$(PYTHON)

uninstall:
	$(PIP) uninstall corrfitter

#	- cat files-corrfitter.$(PYTHON) | xargs rm -rf
#	- cat files-dataset.$(PYTHON) | xargs rm -rf

install-dataset:
	$(PIP) install dataset

#	$(PYTHON) dataset-setup.py install --user --record files-dataset.$(PYTHON)

uninstall-dataset:
	$(PIP) uninstall dataset

try:
	$(PYTHON) setup.py install --user --record files-corrfitter

untry:
	- cat files-corrfitter | xargs rm -rf

.PHONY : tests

tests test-all:
	@echo 'N.B. Some tests involve random numbers and so fail occasionally'
	@echo '     (less than 1 in 100 times) due to multi-sigma fluctuations.'
	@echo '     Run again if any test fails.'
	@echo ''
	cd tests; $(PYTHON) -m unittest discover; cd ..

run run-examples:
	$(MAKE) -C examples run-all

time:
	$(MAKE) -C examples time

sdist:			# source distribution
	$(PYTHON) setup.py sdist

corrfitter.tz:	# everything distribution
	make clean
	tar -C .. --exclude '\.svn' -z -c -v -f corrfitter.tz corrfitter

doc-html:
	make doc/html/index.html

doc/html/index.html : $(DOCFILES) $(SRCFILES)
	rm -rf doc/html; sphinx-build -b html doc/source doc/html

doc-pdf:
	make doc/corrfitter.pdf

doc/corrfitter.pdf : $(DOCFILES) $(SRCFILES)
	rm -rf doc/corrfitter.pdf
	sphinx-build -b latex doc/source doc/latex
	cd doc/latex; make corrfitter.pdf; mv corrfitter.pdf ..

doc-zip doc.zip:
	cd doc/html; zip -r doc *; mv doc.zip ../..

doc-all: doc-html doc-pdf

%.so : %.pyx
	$(PYTHON) $*-setup.py build_ext --inplace

upload-pypi:
	python setup.py sdist upload

upload-twine:
	twine upload dist/corrfitter-$(VERSION).tar.gz

upload-git:
	make doc-html doc-pdf
	git diff --exit-code
	git diff --cached --exit-code
	git push origin master


test-download:
	-$(PIP) uninstall corrfitter
	$(PIP) install corrfitter --no-cache-dir

tag-git:
	echo  "version $(VERSION)"
	git tag -a v$(VERSION) -m "version $(VERSION)"
	git push origin v$(VERSION)

clean:
	cd tests; make clean; cd ..
	cd examples; make clean; cd ..
	rm -f *.pyc *.tmp corrfitter*.tz doc.zip
	rm -rf __pycache__
	rm -rf dist
	rm -rf build
	cd doc/source; make clean



