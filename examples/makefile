PYTHON = python

run run-all run-examples:
	### illustrate corrfitter: Corr2 and Corr3
	-$(PYTHON) example.py > example.tmp	# few sec
	-diff -w -I 't=' example.tmp example.out
	-$(PYTHON) example-marginalize.py > example-marginalize.tmp	# few sec
	-diff -w -I 't=' example-marginalize.tmp example-marginalize.out
	-$(PYTHON) example-chained.py > example-chained.tmp	# few sec
	-diff -w -I 't=' example-chained.tmp example-chained.out
	-$(PYTHON) corr2corr3.py > corr2corr3.tmp	# few sec
	-diff -w -I 't=' corr2corr3.tmp corr2corr3.out
	-$(PYTHON) corr2corr3-chd.py > corr2corr3-chd.tmp	# few sec
	-diff -w -I 't=' corr2corr3-chd.tmp corr2corr3-chd.out
	-$(PYTHON) marginalization.py > marginalization.tmp	# few sec
	-diff -w -I 't=' marginalization.tmp marginalization.out
	-$(PYTHON) marginalization-chd.py > marginalization-chd.tmp	# few sec
	-diff -w -I 't=' marginalization-chd.tmp marginalization-chd.out
	-$(PYTHON) bootstrap.py > bootstrap.tmp	# 6 sec
	-diff -w -I 't=' bootstrap.tmp bootstrap.out
	-$(PYTHON) transposedV.py > transposedV.tmp	# few sec
	-diff -w -I 't=' transposedV.tmp transposedV.out
	-$(PYTHON) symmetricV.py > symmetricV.tmp	# few sec
	-diff -w -I 't=' symmetricV.tmp symmetricV.out
	-$(PYTHON) symmetricV-chd.py > symmetricV-chd.tmp	# few sec
	-diff -w -I 't=' symmetricV-chd.tmp symmetricV-chd.out
	-$(PYTHON) matrix-correlator.py > matrix-correlator.tmp # 10 sec
	-diff -w -I 't=' matrix-correlator.tmp matrix-correlator.out
	-$(PYTHON) matrix-correlator-chd.py > matrix-correlator-chd.tmp # 10 sec
	-diff -w -I 't=' matrix-correlator-chd.tmp matrix-correlator-chd.out
	-$(PYTHON) antiperiodic.py > antiperiodic.tmp # 1 sec
	-diff -w -I 't=' antiperiodic.tmp antiperiodic.out
	### DONE

update-all:
	### illustrate corrfitter: Corr2 and Corr3
	-$(PYTHON) example.py > example.out	# few sec
	-$(PYTHON) corr2corr3.py > corr2corr3.out	# few sec
	-$(PYTHON) corr2corr3-chd.py > corr2corr3-chd.out	# few sec
	-$(PYTHON) marginalization.py > marginalization.out	# few sec
	-$(PYTHON) marginalization-chd.py > marginalization-chd.out	# few sec
	-$(PYTHON) bootstrap.py > bootstrap.out	# 6 sec
	-$(PYTHON) transposedV.py > transposedV.out	# few sec
	-$(PYTHON) symmetricV.py > symmetricV.out	# few sec
	-$(PYTHON) symmetricV-chd.py > symmetricV-chd.out	# few sec
	-$(PYTHON) matrix-correlator.py > matrix-correlator.out # 10 sec
	-$(PYTHON) matrix-correlator-chd.py > matrix-correlator-chd.out # 10 sec
	-$(PYTHON) antiperiodic.py > antiperiodic.out # 1 sec
	### DONE

clean:
	rm -f *.pyc *.p *.tmp *.prof