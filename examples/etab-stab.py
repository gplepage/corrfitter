from __future__ import print_function   # makes this work for python2 and 3

# use etab.py but with different make_prior
import etab
main = etab.main
etab.DISPLAYPLOTS = False       # display plots at end of fits?

def make_prior(N, basis):
    return basis.make_prior(nterm=N, keyfmt='etab.{s1}', states=[0, 1, 2])

etab.make_prior = make_prior

if __name__ == '__main__':
    import gvar as gv
    gv.ranseed(1)
    main()