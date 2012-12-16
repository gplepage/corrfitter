## makefile for corrfitter library
#
# Created by G. Peter Lepage on 2011-09-22.
# Copyright (c) 2011-2012 Cornell University. All rights reserved.

PYTHON = python

install:
	$(PYTHON) setup.py install --user	--record files-corrfitter.$(PYTHON)

install-sys:
	$(PYTHON) setup.py install	--record files-corrfitter.$(PYTHON)

uninstall :			# not sure this works --- be careful
	- cat files-corrfitter.$(PYTHON) | xargs rm -rf
	- cat files-dataset.$(PYTHON) | xargs rm -rf

install-dataset:
	$(PYTHON) dataset-setup.py install --user --record files-dataset.$(PYTHON)

.PHONY : tests

tests test-all:
	@echo 'N.B. Some tests involve random numbers and so fail occasionally'
	@echo '     (less than 1 in 100 times) due to multi-sigma fluctuations.'
	@echo '     Run again if any test fails.'
	@echo ''
	$(PYTHON) -m unittest discover

run-examples:
	$(MAKE) -C examples run-all

test-corrfitter:
	$(MAKE) -C tests test-corrfitter

sdist:			# source distribution
	$(PYTHON) setup.py sdist

corrfitter.tz:	# everything distribution
	make clean
	tar -C .. --exclude '\.svn' -z -c -v -f corrfitter.tz corrfitter

doc-html:		# html version of documentation (in doc/html)
	-mkdir doc/source/_build		# probably exists
	-mkdir doc/source/_static		# probably exists
	-mkdir doc/source/_templates	# probably exists
	rm -rf doc/html; cd doc/source; make html; mv _build/html ..

doc-pdf:		# pdf version of documentation (in doc/)
	-mkdir doc/source/_build    	# probably exists
	-mkdir doc/source/_static   	# probably exists
	-mkdir doc/source/_templates	# probably exists
	cd doc/source; make latex; cd _build/latex;  make all-pdf
	mv doc/source/_build/latex/corrfitter.pdf doc/corrfitter.pdf

doc-all: doc-html doc-pdf

%.so : %.pyx
	$(PYTHON) $*-setup.py build_ext --inplace

clean:
	cd tests; make clean; cd ..
	rm -f *.pyc *.tmp corrfitter*.tz
	rm -rf dist
	rm -rf build
	cd doc/source; make clean



