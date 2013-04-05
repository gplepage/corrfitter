#!/usr/bin/env python
# encoding: utf-8
"""
transposedV.py -- tests transposed V

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
try: 
    import matplotlib
except ImportError:
    DISPLAYPLOTS = False

TEST = True       # testing mode? (True, False, or "dump")

if TEST:
    NEXP_LIST = [3]
    TEST_FILENAME = 'transposedV.testp'
    with open(TEST_FILENAME, "r") as f:
        P0_TEST = BufferDict.load(f, use_json=True)
else:
    NEXP_LIST = [2,3,4,5,6,7,8]

def main():
    dfile = "symmetricV.data"
    data = avg_data(Dataset(dfile))
    pfile = "transposedV.p" # last fit
    fitter = CorrFitter(models=build_models())
    p0 = P0_TEST if TEST else pfile
    for nexp in NEXP_LIST:
        print('========================== nexp =',nexp)
        prior = build_prior(nexp)
        fit = fitter.lsqfit(data=data,prior=prior,p0=p0,svdcut=1e-5,maxit=2000, tol=1e-4)
        print_results(fit,prior,data)
        print('\n\n')
        if TEST == "dump":
            with open(TEST_FILENAME, "w") as f:
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
    aetas = fit.p['etas:a']
    print(' etas:a =',fmtlist(aetas[:3]))
    print()
    ##
    ## Ds parameters ##
    dEDs = exp(fit.p['log(Ds:dE)'])
    print('Ds:dE =',fmtlist(dEDs[:3]))
    print(' Ds:E =',fmtlist([sum(dEDs[:i+1]) for i in range(3)]))
    aDs = fit.p['Ds:a']
    print(' Ds:a =',fmtlist(aDs[:3]))
    print()
    ##
    ## vertices ##
    print('etas->V->Ds  =',fit.p['Vnn'][0,0].fmt())
    print('etas->V->Dso =',fit.p['Vno'][0,0].fmt())
    print()
    ##
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
    ## Ds ##
    prior['Ds:a'] = [gvar(0.01,1.0) for i in range(nexp)]
    prior['log(Ds:dE)'] = [log(gvar(0.5,0.5)) for i in range(nexp)]
    prior['log(Ds:dE)'][0] = log(gvar(1.30,0.2))
    prior['Ds:ao'] = [gvar(0.01,1.0) for i in range(nexp)]
    prior['log(Ds:dEo)'] = [log(gvar(0.5,0.5)) for i in range(nexp)]
    prior['log(Ds:dEo)'][0] = log(gvar(1.40,0.2))
    ##
    ## etas ##
    prior['etas:a'] = [gvar(0.01,1.0) for i in range(nexp)]
    prior['log(etas:dE)'] = [log(gvar(0.5,0.5)) for i in range(nexp)]
    prior['log(etas:dE)'][0] = log(gvar(0.67,0.2))
    prior['etas:ao'] = [gvar(0.01,1.0) for i in range(nexp)]
    prior['log(etas:dEo)'] = [log(gvar(0.5,0.5)) for i in range(nexp)]
    prior['log(etas:dEo)'][0] = log(gvar(0.88,0.2))
    ## V ##
    prior['Vno'] = [[gvar(0.01,1.0) for i in range(nexp)] for j in range(nexp)]
    prior['Vnn'] = [[gvar(0.01,1.0) for i in range(nexp)] for j in range(nexp)]
    prior['Voo'] = [[gvar(0.01,1.0) for i in range(nexp)] for j in range(nexp)]
    prior['Von'] = [[gvar(0.01,1.0) for i in range(nexp)] for j in range(nexp)]
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
            a=('Ds:a','Ds:ao'),b=('Ds:a','Ds:ao'),
            dE=('log(Ds:dE)','log(Ds:dEo)'),s=(1.,-1.)),
        #
        Corr2(datatag='etas',tp=tp,tdata=tdata,tfit=tfit,
            a=('etas:a','etas:ao'),b=('etas:a','etas:ao'),
            dE=('log(etas:dE)','log(etas:dEo)'),s=(1.,-1.)),
        #    
        Corr3(datatag='DsetasT18',tdata=range(19),T=18,tfit=range(tmin,18-tmin),
            b=('Ds:a','Ds:ao'),dEb=('log(Ds:dE)','log(Ds:dEo)'),tpa=tp,sa=(1.,-1),
            a=('etas:a','etas:ao'),dEa=('log(etas:dE)','log(etas:dEo)'),tpb=tp,sb=(1.,-1.),
            Vnn='Vnn',Voo='Voo',Vno='Vno',Von='Vno'),
        #
        Corr3(datatag='DsetasT15',tdata=range(16),T=15,tfit=range(tmin,15-tmin),
            b=('Ds:a','Ds:ao'),dEb=('log(Ds:dE)','log(Ds:dEo)'),tpa=tp,sa=(1.,-1),
            a=('etas:a','etas:ao'),dEa=('log(etas:dE)','log(etas:dEo)'),tpb=tp,sb=(1.,-1.),
            Vnn='Vnn',Voo='Voo',Vno='Vno',Von='Von'),
       #
       Corr3(datatag='etasDsT18',tdata=range(19),T=18,tfit=range(tmin,18-tmin),
           b=('etas:a','etas:ao'),dEb=('log(etas:dE)','log(etas:dEo)'),tpa=tp,sa=(1.,-1),
           a=('Ds:a','Ds:ao'),dEa=('log(Ds:dE)','log(Ds:dEo)'),tpb=tp,sb=(1.,-1.),
           Vnn='Vnn',Voo='Voo',Vno='Vno',Von='Von',transpose_V=True),
       #
       Corr3(datatag='etasDsT15',tdata=range(16),T=15,tfit=range(tmin,15-tmin),
           b=('etas:a','etas:ao'),dEb=('log(etas:dE)','log(etas:dEo)'),tpa=tp,sa=(1.,-1),
           a=('Ds:a','Ds:ao'),dEa=('log(Ds:dE)','log(Ds:dEo)'),tpb=tp,sb=(1.,-1.),
           Vnn='Vnn',Voo='Voo',Vno='Vno',Von='Von',transpose_V=True)
                    
    ][:]
    return models


if __name__ == '__main__':
    main()

