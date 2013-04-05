#!/usr/bin/env python
# encoding: utf-8
"""
corr2corr3.py  --- standard usage (including error budgets and plots)

Created by Peter Lepage on 2010-11-26.
Copyright (c) 2010-2013 G. Peter Lepage.
"""

from __future__ import print_function   # makes this work for python2 and 3

import os
import json
import collections
from corrfitter import Corr2, Corr3, CorrFitter,  fastfit
from gvar import gvar, log, exp, BufferDict, fmt_errorbudget
from gvar.dataset import Dataset, avg_data
from numpy import array, arange
import lsqfit
import copy


DISPLAYPLOTS = False         # display plots at end of fitting
try: 
    import matplotlib
except ImportError:
    DISPLAYPLOTS = False

TEST = True             # testing mode? (True,  False,  or "dump")

FASTFIT = True          # compute effective Es

if TEST:
    NTERMLIST = [6]
    TEST_FILENAME = 'corr2corr3.testp'
    try:
        with open(TEST_FILENAME, "r") as f:
            P0_TEST = BufferDict.load(f,  use_json=True)
    except (IOError,  EOFError):
        P0_TEST = None
else:
    NTERMLIST = [2, 3, 4, 5, 6]
    
def main():
    dfile = "example.data"
    data = avg_data(Dataset(dfile)) # compute avg and cov for data
    fitter = CorrFitter(models=build_models())
    p0 = P0_TEST if TEST else None
    sdata = copy.deepcopy(data)
    for nexp in NTERMLIST:
        print('========================== nexp =',  nexp)
        prior = build_prior(nexp)
        fit = fitter.lsqfit(data=data, prior=prior, p0=p0) # pfile)
        p0 = fit.pmean
        print_results(fit, prior, data)
        print('\n\n')
    if FASTFIT:
        def print_ffit(meson,  ffit,  fit=fit):
            osc = ffit.osc
            elabel = "Eo" if osc else "E"
            alabel = "ao" if osc else "a"
            E = ffit.E 
            if osc:
                E_fit = exp(fit.p["log(%s:dEo)" % meson][0])
                a_fit = exp(fit.p["log(%s:ao)" % meson][0])
            else:
                E_fit = exp(fit.p["log(%s:dE)" % meson][0])
                a_fit = exp(fit.p["log(%s:a)" % meson][0])
            a = (ffit.ampl)**0.5
            chi2_dof = ffit.chi2/ffit.dof
            Q = ffit.Q
            print("%8s: %2s = %s   %2s_fit = %s   chi2/dof = %.2f   Q = %.1f"
                % (meson,  elabel,  E.fmt(),  elabel,  E_fit.fmt(),  chi2_dof[0],  Q[0]))
            if osc:
                print("")
            else:
                print("%8s  %2s = %s   %2s_fit = %s   chi2/dof = %.2f   Q = %.1f\n"
                    % ("",  alabel,  a.fmt(),  alabel,  a_fit.fmt(),  chi2_dof[1],  Q[1]))
        ##
        print('fast_fit Analysis:')
        for k, model in zip(["etas", "Ds"], fitter.models[:2]):
            ffit = fastfit(data=data,  prior=prior,  model=model)
            print_ffit(k,  ffit)
            if k == 'Ds':
                ffit = fastfit(data=data,  prior=prior,  model=model,  osc='True')
                print_ffit(k,  ffit)
    if TEST == 'dump':
        with open(TEST_FILENAME, "w") as f:
            fit.pmean.dump(f,  use_json=True)
    if DISPLAYPLOTS:
        fitter.display_plots()
##
      
def print_results(fit,  prior,  data):
    """ print out additional results from the fit """
    ## etas parameters ##
    dEetas = exp(fit.p['log(etas:dE)'])
    print('etas:dE =', fmtlist(dEetas[:3]))
    print(' etas:E =', fmtlist([sum(dEetas[:i+1]) for i in range(3)]))
    aetas = exp(fit.p['log(etas:a)'])
    print(' etas:a =', fmtlist(aetas[:3]))
    print()
    ##
    ## Ds parameters ##
    dEDs = exp(fit.p['log(Ds:dE)'])
    print('Ds:dE =', fmtlist(dEDs[:3]))
    print(' Ds:E =', fmtlist([sum(dEDs[:i+1]) for i in range(3)]))
    aDs = exp(fit.p['log(Ds:a)'])
    print(' Ds:a =', fmtlist(aDs[:3]))
    print()
    ##
    ## vertices ##
    print('etas->V->Ds  =', fit.p['Vnn'][0, 0].fmt())
    print('etas->V->Dso =', fit.p['Vno'][0, 0].fmt())
    print()
    ##
    ## error budget ##
    outputs = BufferDict()
    outputs['Vnn'] = fit.p['Vnn'][0, 0]
    outputs['Vno'] = fit.p['Vno'][0, 0]
    outputs['Eetas'] = dEetas[0]
    outputs['EDs'] = dEDs[0]

    inputs = collections.OrderedDict()
    inputs['stat.'] = data          # statistical errors in data
    inputs['svd'] = fit.svdcorrection 
    inputs.update(prior)            # all entries in prior,  separately
    
    print(fmt_errorbudget(outputs, inputs, ndecimal=3))
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
    prior['log(etas:a)'] = [log(gvar(0.3, 0.3)) for i in range(nexp)]
    prior['log(etas:dE)'] = [log(gvar(0.5, 0.5)) for i in range(nexp)]
    prior['log(etas:dE)'][0] = log(gvar(0.4, 0.2))
    ##
    ## Ds ##
    prior['log(Ds:a)'] = [log(gvar(0.3, 0.3)) for i in range(nexp)]
    prior['log(Ds:dE)'] = [log(gvar(0.5, 0.5)) for i in range(nexp)]
    prior['log(Ds:dE)'][0] = log(gvar(1.2, 0.2))
    prior['log(Ds:ao)'] = [log(gvar(0.1, 0.1)) for i in range(nexp)]
    prior['log(Ds:dEo)'] = [log(gvar(0.5, 0.5)) for i in range(nexp)]
    prior['log(Ds:dEo)'][0] = log(exp(prior['log(Ds:dE)'][0])+gvar(0.3, 0.3))
    ##
    ## V ##
    prior['Vnn'] = [[gvar(0.01, 1.0) for i in range(nexp)] for j in range(nexp)]
    prior['Vno'] = [[gvar(0.01, 1.0) for i in range(nexp)] for j in range(nexp)]
    return prior
##

def build_models():
    """ build models """
    tmin = 5
    tdata = range(64)
    tfit = range(tmin, 64-tmin) # all ts
    
    tp = 64 # periodic
    models = [
        Corr2(datatag='etas', tp=tp,  tdata=tdata,  tfit=tfit,  
            a='etas:a',  b='etas:a',  dE='etas:dE'),  
            
        Corr2(datatag='Ds', tp=tp, tdata=tdata, tfit=tfit, 
            a=('Ds:a', 'Ds:ao'), b=('Ds:a', 'Ds:ao'), 
            dE=('Ds:dE', 'Ds:dEo'), s=(1., -1.)), 

        Corr3(datatag='3ptT15', tdata=range(16), T=15, tfit=range(tmin, 16-tmin), 
            a='etas:a', dEa='etas:dE', tpa=tp, 
            b=('Ds:a', 'Ds:ao'), dEb=('Ds:dE', 'Ds:dEo'), tpb=tp, sb=(1, -1.), 
            Vnn='Vnn', Vno='Vno', transpose_V=False), 
            
        Corr3(datatag='3ptT16', tdata=range(17), T=16, tfit=range(tmin, 17-tmin), 
            a='etas:a', dEa='etas:dE', tpa=tp, 
            b=('Ds:a', 'Ds:ao'), dEb=('Ds:dE', 'Ds:dEo'), tpb=tp, sb=(1, -1.), 
            Vnn='Vnn', Vno='Vno', transpose_V=False)
    ]
    return models
##    
    
if __name__ == '__main__':
    if True:
        main()
    else:
        import hotshot,  hotshot.stats
        prof = hotshot.Profile("lsqfit-test1.prof")
        prof.runcall(main)
        prof.close()
        stats = hotshot.stats.load("lsqfit-test1.prof")
        stats.strip_dirs()
        stats.sort_stats('time',  'calls')
        stats.print_stats(140)

