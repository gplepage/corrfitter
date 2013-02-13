PYTHON = python

run run-all run-examples:
	### illustrate corrfitter: Corr2 and Corr3
	-$(PYTHON) test-1.py > test-1.tmp	# few sec
	-diff -w -I 't=' test-1.tmp test-1.out
	-$(PYTHON) test-1m.py > test-1m.tmp	# few sec
	-diff -w -I 't=' test-1m.tmp test-1m.out
	-$(PYTHON) test-1b.py > test-1b.tmp	# 6 sec
	-diff -w -I 't=' test-1b.tmp test-1b.out
	-$(PYTHON) test-1t.py > test-1t.tmp	# few sec
	-diff -w -I 't=' test-1t.tmp test-1t.out
	-$(PYTHON) test-1s.py > test-1s.tmp	# few sec
	-diff -w -I 't=' test-1s.tmp test-1s.out
	-$(PYTHON) test-ups.py > test-ups.tmp # 10 sec
	-diff -w -I 't=' test-ups.tmp test-ups.out
	-$(PYTHON) test-antiperiodic.py > test-antiperiodic.tmp # 1 sec
	-diff -w -I 't=' test-antiperiodic.tmp test-antiperiodic.out
	### DONE


clean:
	rm -f *.pyc *.p *.tmp *.prof