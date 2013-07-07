## makefile for corrfitter library
#
# Created by G. Peter Lepage (Cornell University) on 2011-09-22.
# Copyright (c) 2011-2013 G. Peter Lepage.

PYTHON = python

install:
	$(PYTHON) setup.py install --user	--record files-corrfitter.$(PYTHON)

install-sys:
	$(PYTHON) setup.py install	--record files-corrfitter.$(PYTHON)

uninstall:			# not sure this works --- be careful
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

clean:
	cd tests; make clean; cd ..
	cd examples; make clean; cd ..
	rm -f *.pyc *.tmp corrfitter*.tz doc.zip
	rm -rf __pycache__
	rm -rf dist
	rm -rf build
	cd doc/source; make clean



