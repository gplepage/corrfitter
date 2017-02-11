PYTHON = python
PLOT = True

run run-all run-examples:
	-$(PYTHON) etas.py > etas.tmp
	-diff -w -I 't=' -I svdcut etas.tmp etas.out
	-$(PYTHON) etas-Ds.py $(PLOT) > etas-Ds.tmp
	-diff -w -I 't=' -I svdcut etas-Ds.tmp etas-Ds.out
	-$(PYTHON) etas-Ds-marginalize.py > etas-Ds-marginalize.tmp
	-diff -w -I 't=' -I svdcut etas-Ds-marginalize.tmp etas-Ds-marginalize.out
	-$(PYTHON) etas-Ds-chained.py > etas-Ds-chained.tmp
	-diff -w -I 't=' -I svdcut etas-Ds-chained.tmp etas-Ds-chained.out
	-$(PYTHON) Ds-Ds.py > Ds-Ds.tmp
	-diff -w -I 't=' -I svdcut Ds-Ds.tmp Ds-Ds.out
	-$(PYTHON) etab.py > etab.tmp
	-diff -w -I 't=' -I svdcut etab.tmp etab.out
	-$(PYTHON) etab-alt.py > etab-alt.tmp
	-diff -w -I 't=' -I svdcut etab-alt.tmp etab-alt.out
	-$(PYTHON) etab-stab.py > etab-stab.tmp
	-diff -w -I 't=' -I svdcut etab-stab.tmp etab-stab.out

time:
	time $(MAKE) _time PLOT=False

_time:
	-$(PYTHON) etas.py > /dev/null
	-$(PYTHON) etas-Ds.py $(PLOT) > /dev/null
	-$(PYTHON) etas-Ds-marginalize.py > /dev/null
	-$(PYTHON) etas-Ds-chained.py > /dev/null
	-$(PYTHON) Ds-Ds.py > /dev/null
	-$(PYTHON) etab.py > /dev/null
	-$(PYTHON) etab-alt.py > /dev/null
	-$(PYTHON) etab-stab.py > /dev/null

update-all:
	-$(PYTHON) etas.py > etas.out
	-$(PYTHON) etas-Ds.py > etas-Ds.out
	-$(PYTHON) etas-Ds-marginalize.py > etas-Ds-marginalize.out
	-$(PYTHON) etas-Ds-chained.py > etas-Ds-chained.out
	-$(PYTHON) Ds-Ds.py > Ds-Ds.out
	-$(PYTHON) etab.py > etab.out
	-$(PYTHON) etab-alt.py > etab-alt.out
	-$(PYTHON) etab-stab.py > etab-stab.out

run-extras:
	-$(PYTHON) extra-fits.py > extra-fits.tmp
	-diff -w -I 't=' -I svdcut extra-fits.tmp extra-fits.out

clean:
	rm -f *.pyc *.p *.tmp *.prof