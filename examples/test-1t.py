#!/usr/bin/env python
# encoding: utf-8
"""
test-1t.py

Created by Peter Lepage on 2010-11-26.
Copyright (c) 2010-2012 Cornell University. All rights reserved.
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
    TEST_FILENAME = 'test-1t.testp'
    with open(TEST_FILENAME, "r") as f:
        P0_TEST = BufferDict.load(f, use_json=True)
else:
    NEXP_LIST = [2,3,4,5,6,7,8]

def main():
    # dfile = "coarse_3pt_etas_etas_kinCth1.10.576.bint8"  # data file
    dfile = "coarse_3pt_1link_Ds_and_etas_kinCth1.10.576.bint8"
    data = avg_data(Dataset(dfile))
    pfile = "test-1t.p" # last fit
    fitter = CorrFitter(models=build_models())
    # fitdata = fitter.compute_fitdata(dfile)
    p0 = P0_TEST if TEST else pfile
    for nexp in NEXP_LIST:
        print('========================== nexp =',nexp)
        prior = build_prior(nexp)
        fit = fitter.lsqfit(data=data,prior=prior,p0=p0,svdcut=1e-5,maxit=2000, tol=1e-4)
        print_results(fit,prior,data)
        print('\n\n')
        if TEST == "dump":
            fit.dump_pmean(TEST_FILENAME)
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
##    

P0_TEST = {'etas:ao': [0.06832113722650532, -0.04856424525422174,
0.13244215104455168], 'Vno': [[-0.15780975926841767, -0.10138450970109109,
0.03754390993676524], [-0.03028043405104003, -0.0014422366769212467,
-0.07590360369715335], [-0.0844583272134997, -0.03883915793203409,
0.01072279950831445]], 'Vnn': [[0.15865817872322574, -0.09958219668240623,
-0.003966513960031803], [0.009499363982511983, -0.381942448892726,
0.02537795274325082], [-0.0008979136255186132, 0.014628696247306049,
0.010271401780123526]], 'Von': [[-0.16124675838768052, -0.09054606786355532,
0.03453282386932355], [0.06621806750891938, 0.2560712079116751,
-0.044590086144259654], [-0.038871791527634744, -0.2301318831852351,
0.007420921215960247]], 'Voo': [[-0.018833837124288957, -0.30304599085542405,
0.3095056509997733], [0.5786174641162435, 0.13425658147144381,
0.007819484719777179], [0.4467342819268318, 0.054768171268679826,
0.009252848304865984]], 'Ds:a': [0.21075652444233753, 0.22193845186466948,
0.512739124744146], 'log(Ds:dEo)': [0.3188343293274145, -1.8021978670558103,
-1.083380431182324], 'log(etas:dE)': [-0.3973277828215011,
-0.4963863230507786, -0.06851152665600153], 'Ds:ao': [0.04467292682632941,
0.08747341988514248, 0.11646484365384113], 'log(Ds:dE)': [0.271004632904409,
-0.8663503670559614, -0.7773824952625505], 'log(etas:dEo)':
[-0.15639760080394724, -0.574443182323812, -0.4014766435400238], 'etas:a':
[0.17332580380244753, 0.2747232632782871, 0.6287293298724147]}

for k in P0_TEST:
    P0_TEST[k] = array(P0_TEST[k])


if __name__ == '__main__':
    main()

