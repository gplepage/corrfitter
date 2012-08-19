#!/usr/bin/env python
# encoding: utf-8
"""
test-1.py

Created by Peter Lepage on 2010-11-26.
Copyright (c) 2010/2011 Cornell University. All rights reserved.
"""

import os
from corrfitter import Corr2,Corr3,CorrFitter
from gvar import gvar,log,exp,evalcov,BufferDict,fmt_errorbudget
from gvar.dataset import Dataset,avg_data
from numpy import array,arange,dot

import lsqfit

DISPLAYPLOTS = False         # display plots at end of fitting
try: 
    import matplotlib
except ImportError:
    DISPLAYPLOTS = False

TEST = True         # testing mode? (True, False, or "dump")

if TEST:
    NEXP_LIST = [3]
    TEST_FILENAME = 'test-1s.testp'
    P0_TEST = lsqfit.nonlinear_fit.load_parameters(TEST_FILENAME)
else:
    NEXP_LIST = [2,3,4,5,6]

def main():
    # dfile = "coarse_3pt_etas_etas_kinCth1.10.576.bint8"  # data file
    dfile = "coarse_3pt_1link_Ds_and_etas_kinCth1.10.576.bint8"
    data = avg_data(Dataset(dfile,keys=['DsDsT18','DsDsT15','Ds']))
    fitter = CorrFitter(models=build_models())
    pfile = "test-1s.p" # last fit
    p0 = P0_TEST if TEST else pfile
    for nexp in NEXP_LIST:
        print '========================== nexp =',nexp
        prior = build_prior(nexp)
        fit = fitter.lsqfit(data=data,prior=prior,p0=p0,debug=True)
        print_results(fit,prior,data)
        print '\n\n'
        if TEST == 'dump':
            fit.dump_pmean(TEST_FILENAME)
    if DISPLAYPLOTS:
        fitter.display_plots()
##
      
def print_results(fit,prior,data):
    """ print out additional results from the fit """
    dE = exp(fit.p['logdE'])
    print 'dE =',fmtlist(dE[:3])
    print ' E =',fmtlist([sum(dE[:i+1]) for i in range(3)])
    V = fit.p['Vnn'][0]
    print 'm->V->m =',V.fmt()
    print
    outputs = {'V':V}
    inputs = {'stat':data}          # all data statistical errors
    inputs.update(prior)            # all parts of the prior, separately
    print fmt_errorbudget(outputs,inputs)
##

def fmtlist(x):
    return '  '.join([xi.fmt() for xi in x])
##
    
def build_prior(nexp):
    """ build prior """
    prior = BufferDict()
    ## etas ##
    prior['a'] = [gvar(0.01,1.0) for i in range(nexp)]
    prior['logdE'] = [log(gvar(0.5,0.5)) for i in range(nexp)]
    prior['logdE'][0] = log(gvar(1.40,0.2))
    prior['ao'] = [gvar(0.01,1.0) for i in range(nexp)]
    prior['logdEo'] = [log(gvar(0.5,0.5)) for i in range(nexp)]
    prior['logdEo'][0] = log(gvar(1.59,0.2))
    ##
    ## V ##
    nV = (nexp*(nexp+1))/2
    prior['Vnn'] = [gvar(0.01,1.0) for i in range(nV)]
    prior['Voo'] = [gvar(0.01,1.0) for i in range(nV)]
    prior['Vno'] = [[gvar(0.01,1.0) for i in range(nexp)] for j in range(nexp)]
    return prior
##

def build_models():
    """ build models """
    tmin = 2
    tdata = range(64)
    tfit = range(tmin,64-tmin) # all ts
    
    tp = 64 # periodic
    models = [
        Corr2(datatag='Ds',tp=tp,tdata=tdata,tfit=tfit,
            a=('a','ao'),b=('a','ao'),dE=('logdE','logdEo'),s=(1.,-1.)),
            
        Corr3(datatag='DsDsT18',tdata=range(19),T=18,tfit=range(tmin,18-tmin),
            a=('a','ao'),dEa=('logdE','logdEo'),tpa=tp,sa=(1.,-1),
            b=('a','ao'),dEb=('logdE','logdEo'),tpb=tp,sb=(1.,-1.),
            Vnn='Vnn',Voo='Voo',Vno='Vno',Von='Vno',symmetric_V=True),
        #
        Corr3(datatag='DsDsT15',tdata=range(16),T=15,tfit=range(tmin,15-tmin),
            a=('a','ao'),dEa=('logdE','logdEo'),tpa=tp,sa=(1.,-1),
            b=('a','ao'),dEb=('logdE','logdEo'),tpb=tp,sb=(1.,-1.),
            Vnn='Vnn',Voo='Voo',Vno='Vno',Von='Vno',symmetric_V=True)
                    
    ][:]
    return models
##    


if __name__ == '__main__':
    main()

