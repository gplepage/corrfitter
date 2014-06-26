PYTHON = python

run run-all run-examples:
	-$(PYTHON) etas.py > etas.tmp	
	-diff -w -I 't=' -I svdcut etas.tmp etas.out
	-$(PYTHON) etas-Ds.py > etas-Ds.tmp	
	-diff -w -I 't=' -I svdcut etas-Ds.tmp etas-Ds.out
	-$(PYTHON) etas-Ds-marginalize.py > etas-Ds-marginalize.tmp	
	-diff -w -I 't=' -I svdcut etas-Ds-marginalize.tmp etas-Ds-marginalize.out
	-$(PYTHON) etas-Ds-chained.py > etas-Ds-chained.tmp	
	-diff -w -I 't=' -I svdcut etas-Ds-chained.tmp etas-Ds-chained.out
	-$(PYTHON) etab.py > etab.tmp	
	-diff -w -I 't=' -I svdcut etab.tmp etab.out

update-all:
	-$(PYTHON) etas.py > etas.out	
	-$(PYTHON) etas-Ds.py > etas-Ds.out	
	-$(PYTHON) etas-Ds-marginalize.py > etas-Ds-marginalize.out	
	-$(PYTHON) etas-Ds-chained.py > etas-Ds-chained.out	
	-$(PYTHON) etab.py > etab.out	

clean:
	rm -f *.pyc *.p *.tmp *.prof