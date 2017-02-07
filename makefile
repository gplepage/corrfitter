## makefile for corrfitter library
#
# Created by G. Peter Lepage (Cornell University) on 2011-09-22.
# Copyright (c) 2011-2013 G. Peter Lepage.

PYTHON = python
PIP = $(PYTHON) -m pip
VERSION = `python -c 'import corrfitter; print corrfitter.__version__'`

install:
	$(PIP) install . --user

#	$(PYTHON) setup.py install --user	--record files-corrfitter.$(PYTHON)

install-sys:
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
	$(PYTHON) setup.py install --user --record files-corrfitter.$(PYTHONVERSION)

untry:
	- cat files-corrfitter.$(PYTHONVERSION) | xargs rm -rf

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

doc-html:		# html version of documentation (in doc/html)
	sphinx-build -b html doc/source doc/html

doc-zip doc.zip:
	cd doc/html; zip -r doc *; mv doc.zip ../..

doc-pdf:		# pdf version of documentation (in doc/)
	sphinx-build -b latex doc/source doc/latex
	cd doc/latex; make corrfitter.pdf; mv corrfitter.pdf ..

doc-all: doc-html doc-pdf

%.so : %.pyx
	$(PYTHON) $*-setup.py build_ext --inplace

upload-pypi:
	python setup.py sdist upload

upload-git:
	make doc-all
	git commit -a -m "prep for upload"
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



