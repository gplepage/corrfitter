
#!/usr/bin/env python
# encoding: utf-8
"""
marginalization-chd.py -- tests marginalization w. chained fits

Created by Peter Lepage on 2010-11-26.
Copyright (c) 2010-2013 G. Peter Lepage.
"""

from __future__ import print_function   # makes this work for python2 and 3

import os
import collections
from corrfitter import Corr2,Corr3,CorrFitter
from gvar import gvar,log,exp,BufferDict,fmt_errorbudget
from gvar.dataset import Dataset,avg_data
from numpy import array,arange

import lsqfit

DISPLAYPLOTS = False         # display plots at end of fitting

TEST = True        # run test case: True, False, or "dump"

if TEST:
    TEST_FILE = {True:"marginalization-chd.truep", False:"marginalization-chd.falsep"}
    P0_TEST = {}
    for ratio in [True,False]:
        try:
            with open(TEST_FILE[ratio], "r") as f:
                P0_TEST[ratio] = BufferDict.load(f, use_json=True)
        except (IOError, EOFError):
            P0_TEST[ratio] = None
else:
    P0_TEST = {True:None, False:None}

try: 
    import matplotlib
except ImportError:
    DISPLAYPLOTS = False

def main():
    dfile = "example.data"  # data file
    data = avg_data(Dataset(dfile))
    prior = build_prior(6)              # 6 terms in prior
    nexp = 1                            # nexp terms in fit
    for ratio in [True,False][:]:
    # ratio = True
    # for nexp in [1,2,3,4,5,6][:1]:
        fitter = CorrFitter(models=build_models(),nterm=(nexp,nexp),ratio=ratio)
        print('========================== nexp =',nexp,'  ratio =',ratio)
        fit = fitter.chained_lsqfit(data=data,prior=prior,p0=P0_TEST[ratio])
        print_results(fit,prior,data)
        print('\n\n')
        if TEST == "dump":
            with open(TEST_FILE[ratio], "w") as f:
                fit.pmean.dump(f, use_json=True)
    if DISPLAYPLOTS:
        fitter.display_plots()
##
      
def print_results(fit,prior,data):
    """ print out additional results from the fit """
    dEetas = exp(fit.p['log(etas:dE)'])
    print('etas:dE =',fmtlist(dEetas[:3]))
    print(' etas:E =',fmtlist([sum(dEetas[:i+1]) for i in range(3)]))
    print(' etas:a =',fmtlist(exp(fit.p['log(etas:a)'])))
    dEDs = exp(fit.p['log(Ds:dE)'])
    print()
    print('Ds:dE =',fmtlist(dEDs[:3]))
    print(' Ds:E =',fmtlist([sum(dEDs[:i+1]) for i in range(3)]))
    print(' Ds:a =',fmtlist(exp(fit.p['log(Ds:a)'])))
    print()
    print('etas->V->Ds =',fit.p['Vnn'][0,0].fmt())
    print('etas->V->Dso =',fit.p['Vno'][0,0].fmt())
    print()
    ## error budget ##
    outputs = BufferDict()
    outputs['Vnn'] = fit.p['Vnn'][0,0]
    outputs['Vno'] = fit.p['Vno'][0,0]
    outputs['Eetas'] = dEetas[0]
    outputs['EDs'] = dEDs[0]

    inputs = collections.OrderedDict()
    inputs['stat.'] = data          # statistical errors in data
    # inputs['svd'] = fit.svdcorrection 
    inputs.update(prior)            # all entries in prior, separately
    
    # outputs = dict(Vnn=fit.p['Vnn'][0,0],Vno=fit.p['Vno'][0,0],
    #              Eetas=dEetas[0],EDs=dEDs[0])
    # inputs = {'stat.':[data[k] for k in data]} # statistical errors in data
    # inputs.update(prior)                       # all entries in prior
    print(fmt_errorbudget(outputs,inputs,ndecimal=3))
    ##
##

def fmtlist(x):
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
    prior['Vnn'] = [[gvar(0.1,1.0) for i in range(nexp)] for j in range(nexp)]
    prior['Vno'] = [[gvar(0.1,1.0) for i in range(nexp)] for j in range(nexp)]
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
            b=('log(Ds:a)','log(Ds:ao)'),dEb=('log(Ds:dE)','log(Ds:dEo)'),
            tpb=tp,sb=(1,-1.),
            Vnn='Vnn',Vno='Vno',transpose_V=False),
            
        Corr3(datatag='3ptT16',tdata=range(17),T=16,tfit=range(tmin,17-tmin),
            a='log(etas:a)',dEa='log(etas:dE)',tpa=tp,
            b=('log(Ds:a)','log(Ds:ao)'),dEb=('log(Ds:dE)','log(Ds:dEo)'),
            tpb=tp,sb=(1,-1.),
            Vnn='Vnn',Vno='Vno',transpose_V=False)
    ]
    return models
##    

if __name__ == '__main__':
    main()

