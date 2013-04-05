PYTHON = python

run run-all run-examples:
	### illustrate corrfitter: Corr2 and Corr3
	-$(PYTHON) example.py > example.tmp	# few sec
	-diff -w -I 't=' example.tmp example.out
	-$(PYTHON) example-marginalize.py > example-marginalize.tmp	# few sec
	-diff -w -I 't=' example-marginalize.tmp example-marginalize.out
	-$(PYTHON) example-chained.py > example-chained.tmp	# few sec
	-diff -w -I 't=' example-chained.tmp example-chained.out
	-$(PYTHON) test-1.py > test-1.tmp	# few sec
	-diff -w -I 't=' test-1.tmp test-1.out
	-$(PYTHON) test-1-chd.py > test-1-chd.tmp	# few sec
	-diff -w -I 't=' test-1-chd.tmp test-1-chd.out
	-$(PYTHON) test-1m.py > test-1m.tmp	# few sec
	-diff -w -I 't=' test-1m.tmp test-1m.out
	-$(PYTHON) test-1m-chd.py > test-1m-chd.tmp	# few sec
	-diff -w -I 't=' test-1m-chd.tmp test-1m-chd.out
	-$(PYTHON) test-1b.py > test-1b.tmp	# 6 sec
	-diff -w -I 't=' test-1b.tmp test-1b.out
	-$(PYTHON) test-1t.py > test-1t.tmp	# few sec
	-diff -w -I 't=' test-1t.tmp test-1t.out
	-$(PYTHON) test-1s.py > test-1s.tmp	# few sec
	-diff -w -I 't=' test-1s.tmp test-1s.out
	-$(PYTHON) test-1s-chd.py > test-1s-chd.tmp	# few sec
	-diff -w -I 't=' test-1s-chd.tmp test-1s-chd.out
	-$(PYTHON) test-ups.py > test-ups.tmp # 10 sec
	-diff -w -I 't=' test-ups.tmp test-ups.out
	-$(PYTHON) test-ups-chd.py > test-ups-chd.tmp # 10 sec
	-diff -w -I 't=' test-ups-chd.tmp test-ups-chd.out
	-$(PYTHON) test-antiperiodic.py > test-antiperiodic.tmp # 1 sec
	-diff -w -I 't=' test-antiperiodic.tmp test-antiperiodic.out
	### DONE

update-all:
	### illustrate corrfitter: Corr2 and Corr3
	-$(PYTHON) example.py > example.out	# few sec
	-$(PYTHON) test-1.py > test-1.out	# few sec
	-$(PYTHON) test-1-chd.py > test-1-chd.out	# few sec
	-$(PYTHON) test-1m.py > test-1m.out	# few sec
	-$(PYTHON) test-1m-chd.py > test-1m-chd.out	# few sec
	-$(PYTHON) test-1b.py > test-1b.out	# 6 sec
	-$(PYTHON) test-1t.py > test-1t.out	# few sec
	-$(PYTHON) test-1s.py > test-1s.out	# few sec
	-$(PYTHON) test-1s-chd.py > test-1s-chd.out	# few sec
	-$(PYTHON) test-ups.py > test-ups.out # 10 sec
	-$(PYTHON) test-ups-chd.py > test-ups-chd.out # 10 sec
	-$(PYTHON) test-antiperiodic.py > test-antiperiodic.out # 1 sec
	### DONE

clean:
	rm -f *.pyc *.p *.tmp *.prof