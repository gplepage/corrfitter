#!/usr/bin/env python
# encoding: utf-8
"""
antiperiodic.py   --- tests antiperiodic fits

Created by Peter Lepage on 2012-12-02.
Copyright (c) 2010-2013 G. Peter Lepage.
"""
from __future__ import print_function   # makes this work for python2 and 3

import os
import pickle
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

TEST = True        # testing mode? (True, False, or "dump")

TP = dict(pp=64, pa=-64)        # pa is anti-periodic, pp periodic
TMIN = 5
SVDCUT = 1e-3

if TEST:
    NEXP_LIST = [2]
    TEST_FILENAME = 'antiperiodic.testp'
    with open(TEST_FILENAME,"r") as f:
        P0_TEST = BufferDict.load(f, use_json=True)
else:
    NEXP_LIST = [2]
    
def main():
    dfile = "antiperiodic.data"   # data file
    data = avg_data(Dataset(dfile).trim()) # compute avg and cov for data
    pfile = "antiperiodic.p"  # last fit stored here if not TEST
    fitter = CorrFitter(models=build_models())
    p0 = P0_TEST if TEST else pfile
    for nexp in NEXP_LIST:
        print('========================== nexp =',nexp)
        prior = build_prior(nexp)
        fit = fitter.lsqfit(data=data, prior=prior, p0=p0, svdcut=SVDCUT) # pfile)
        print_results(fit,prior,data)
        print('\n\n')
        if TEST == 'dump':
            with open(TEST_FILENAME,"w") as f:
                fit.pmean.dump(f, use_json=True)
    if DISPLAYPLOTS:
        fitter.display_plots()
##

def print_results(fit, data, prior):
    """docstring for print_results"""
    print('E0 =',exp(fit.p['log(dE)'][0]).fmt(),'\n')
    
##

def build_models():
    def ampl(a):
        if a == 'p':
            return ('log(p)','po')
        else:
            return (a,a+'o')
    ##
    logdE = ('log(dE)','log(dEo)')
    models = []
    tdata = range(64)
    tfit = range(TMIN, 64-TMIN)
    for k in ['pp','pa']:
        models.append(Corr2(datatag='etas.'+k, tp=TP[k], tdata=tdata, 
            tfit=tfit, a=ampl(k[0]), b=ampl(k[1]), dE=logdE, s=(1,-1)))
    return models
##

def build_prior(nexp):
   prior = BufferDict()
   prior['a'] = gvar(nexp*["0(1)"]) 
   prior['log(p)'] = log(gvar(nexp*["0.5(5)"]))
   dE = gvar(nexp*["0.70(35)"])
   dE[0] = gvar("0.44(10)")
   prior['log(dE)'] = log(dE)
   prior['ao'] = gvar(nexp*["0(1)"]) 
   prior['po'] = gvar(nexp*["0(1)"])
   dEo = gvar(nexp*["0.5(5)"])
   prior['log(dEo)'] = log(dEo)
   return prior
##


if __name__ == '__main__':
	main()

