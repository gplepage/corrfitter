#!/usr/bin/env python
# encoding: utf-8
"""
test-1.py  --- standard usage (including error budgets and plots)

Created by Peter Lepage on 2010-11-26.
Copyright (c) 2010-2012 Cornell University. All rights reserved.
"""

from __future__ import print_function   # makes this work for python2 and 3

import os
import json
from corrfitter import Corr2,Corr3,CorrFitter
from gvar import gvar,log,exp,BufferDict,fmt_errorbudget
from gvar.dataset import Dataset,avg_data
from numpy import array,arange
import lsqfit


DISPLAYPLOTS = True         # display plots at end of fitting
try: 
    import matplotlib
except ImportError:
    DISPLAYPLOTS = False

TEST = True        # testing mode? (True, False, or "dump")

if TEST:
    NEXP_LIST = [6]
    TEST_FILENAME = 'test-1.testp'
    # P0_TEST = lsqfit.nonlinear_fit.load_parameters(TEST_FILENAME)
    with open(TEST_FILENAME,"r") as f:
        P0_TEST = BufferDict.load(f, use_json=True)
else:
    TEST_FILENAME = 'test-1.testp'
    P0_TEST = None
    NEXP_LIST = [2,3,4,5,6]
    
def main():
    dfile = "coarse_3pt_etas_Ds.bint8"   # data file
    data = avg_data(Dataset(dfile)) # compute avg and cov for data
    pfile = "test-1.p"                   # last fit stored in file test-1.p
    fitter = CorrFitter(models=build_models())
    p0 = P0_TEST if TEST else pfile
    for nexp in NEXP_LIST:
        print('========================== nexp =',nexp)
        prior = build_prior(nexp)
        fit = fitter.lsqfit(data=data,prior=prior,p0=p0) # pfile)
        print_results(fit,prior,data)
        print('\n\n')
    if TEST == 'dump':
        with open(TEST_FILENAME,"w") as f:
            fit.pmean.dump(f, use_json=True)
    if DISPLAYPLOTS:
        fitter.display_plots()
##
      
def print_results(fit,prior,data):
    """ print out additional results from the fit """
    ## etas parameters ##
    dEetas = exp(fit.p['log(etas:dE)'])
    print('etas:dE =',fmtlist(dEetas[:3]))
    print(' etas:E =',fmtlist([sum(dEetas[:i+1]) for i in range(3)]))
    aetas = exp(fit.p['log(etas:a)'])
    print(' etas:a =',fmtlist(aetas[:3]))
    print()
    ##
    ## Ds parameters ##
    dEDs = exp(fit.p['log(Ds:dE)'])
    print('Ds:dE =',fmtlist(dEDs[:3]))
    print(' Ds:E =',fmtlist([sum(dEDs[:i+1]) for i in range(3)]))
    aDs = exp(fit.p['log(Ds:a)'])
    print(' Ds:a =',fmtlist(aDs[:3]))
    print()
    ##
    ## vertices ##
    print('etas->V->Ds  =',fit.p['Vnn'][0,0].fmt())
    print('etas->V->Dso =',fit.p['Vno'][0,0].fmt())
    print()
    ##
    ## error budget ##
    outputs = dict(Vnn=fit.p['Vnn'][0,0],Vno=fit.p['Vno'][0,0],
                 Eetas=dEetas[0],EDs=dEDs[0])
    inputs = {'stat.':data, 'svd':fit.svdcorrection} # statistical errors in data
    inputs.update(prior)            # all entries in prior, separately
    print(fmt_errorbudget(outputs,inputs,ndigit=3))
    ##
##

def fmtlist(x):
    """ Make neat lists for printing. """
    return '  '.join([xi.fmt() for xi in x])
##
    
def build_prior(nexp):
    """ build prior """
    prior = BufferDict()
    ## etas ##
    prior['log(etas:a)'] = [log(gvar(0.3,0.3)) for i in range(nexp)]
    prior['log(etas:dE)'] = [log(gvar(0.5,0.5)) for i in range(nexp)]
    prior['log(etas:dE)'][0] = log(gvar(0.4,0.2))
    ##
    ## Ds ##
    prior['log(Ds:a)'] = [log(gvar(0.3,0.3)) for i in range(nexp)]
    prior['log(Ds:dE)'] = [log(gvar(0.5,0.5)) for i in range(nexp)]
    prior['log(Ds:dE)'][0] = log(gvar(1.2,0.2))
    prior['log(Ds:ao)'] = [log(gvar(0.1,0.1)) for i in range(nexp)]
    prior['log(Ds:dEo)'] = [log(gvar(0.5,0.5)) for i in range(nexp)]
    prior['log(Ds:dEo)'][0] = log(exp(prior['log(Ds:dE)'][0])+gvar(0.3,0.3))
    ##
    ## V ##
    prior['Vnn'] = [[gvar(0.01,1.0) for i in range(nexp)] for j in range(nexp)]
    prior['Vno'] = [[gvar(0.01,1.0) for i in range(nexp)] for j in range(nexp)]
    return prior
##

def build_models():
    """ build models """
    tmin = 5
    tdata = range(64)
    tfit = range(tmin,64-tmin) # all ts
    
    tp = 64 # periodic
    models = [
        Corr2(datatag='etas',tp=tp,tdata=tdata,tfit=tfit,
            a='log(etas:a)',b='log(etas:a)',dE='log(etas:dE)'),
            
        Corr2(datatag='Ds',tp=tp,tdata=tdata,tfit=tfit,
            a=('log(Ds:a)','log(Ds:ao)'),b=('log(Ds:a)','log(Ds:ao)'),
            dE=('log(Ds:dE)','log(Ds:dEo)'),s=(1.,-1.)),

        Corr3(datatag='3ptT15',tdata=range(16),T=15,tfit=range(tmin,16-tmin),
            a='log(etas:a)',dEa='log(etas:dE)',tpa=tp,
            b=('log(Ds:a)','log(Ds:ao)'),dEb=('log(Ds:dE)','log(Ds:dEo)'),tpb=tp,sb=(1,-1.),
            Vnn='Vnn',Vno='Vno',transpose_V=False),
            
        Corr3(datatag='3ptT16',tdata=range(17),T=16,tfit=range(tmin,17-tmin),
            a='log(etas:a)',dEa='log(etas:dE)',tpa=tp,
            b=('log(Ds:a)','log(Ds:ao)'),dEb=('log(Ds:dE)','log(Ds:dEo)'),tpb=tp,sb=(1,-1.),
            Vnn='Vnn',Vno='Vno',transpose_V=False)
    ]
    return models
##    
    
if __name__ == '__main__':
    if True:
        main()
    else:
        import hotshot, hotshot.stats
        prof = hotshot.Profile("lsqfit-test1.prof")
        prof.runcall(main)
        prof.close()
        stats = hotshot.stats.load("lsqfit-test1.prof")
        stats.strip_dirs()
        stats.sort_stats('time', 'calls')
        stats.print_stats(140)
