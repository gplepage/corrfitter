#!/usr/bin/env python
# encoding: utf-8
"""
test-1.py

Created by Peter Lepage on 2010-11-26.
Copyright (c) 2010/2011 Cornell University. All rights reserved.
"""

import os
from corrfitter import Corr2,Corr3,CorrFitter
from gvar import gvar,log,exp,evalcov,Dataset,BufferDict,avg_data
from numpy import array,arange,dot

import lsqfit

lsqfit.nonlinear_fit.fmt_parameter = '%7.3f +- %7.3f'

DISPLAYPLOTS = False         # display plots at end of fitting
try: 
    import matplotlib
except ImportError:
    DISPLAYPLOTS = False

TEST = True     # testing mode?

def main():
    # dfile = "coarse_3pt_etas_etas_kinCth1.10.576.bint8"  # data file
    dfile = "coarse_3pt_1link_Ds_and_etas_kinCth1.10.576.bint8"
    data = avg_data(Dataset(dfile,keys=['DsDsT18','DsDsT15','Ds']))
    fitter = CorrFitter(models=build_models())
    pfile = "test-1s.p" # last fit
    nexp_list = NEXP_TEST if TEST else [2,3,4,5,6,7,8]
    p0 = P0_TEST if TEST else pfile
    for nexp in nexp_list:
        print '========================== nexp =',nexp
        prior = build_prior(nexp)
        fit = fitter.lsqfit(data=data,prior=prior,p0=p0)
        print_results(fit,prior,data)
        print '\n\n'
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
    inputs = {'stat':[data[k] for k in data]} # all data statistical errors
    inputs.update(prior)                      # all parts of the prior
    print fit.fmt_partialsdev(outputs,inputs)
##

def fmtlist(x):
    return '  '.join([xi.fmt() for xi in x])
##
    
def build_prior(nexp):
    """ build prior """
    prior = BufferDict()
    ## etas ##
    prior.add('a',[gvar(0.01,1.0) for i in range(nexp)])
    prior.add('logdE',[log(gvar(0.5,0.5)) for i in range(nexp)])
    prior['logdE'][0] = log(gvar(1.40,0.2))
    prior.add('ao',[gvar(0.01,1.0) for i in range(nexp)])
    prior.add('logdEo',[log(gvar(0.5,0.5)) for i in range(nexp)])
    prior['logdEo'][0] = log(gvar(1.59,0.2))
    ##
    ## V ##
    nV = (nexp*(nexp+1))/2
    prior.add('Vnn',[gvar(0.01,1.0) for i in range(nV)])
    prior.add('Voo',[gvar(0.01,1.0) for i in range(nV)])
    prior.add('Vno',[[gvar(0.01,1.0) for i in range(nexp)] for j in range(nexp)])
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

NEXP_TEST = [3]
P0_TEST = {'a': [0.21003188674777384, 0.20199445922587544,
0.5183807812003373], 'logdE': [0.27069323026592135, -0.9437824237909911,
-0.748844103009219], 'Vno': [[-0.21021700308173005, -0.008004198284218224,
0.0033815712672157078], [0.038037979385456035, -0.11929504858159,
0.011472143959295505], [0.002271934466225292, 0.1040458783723354,
0.009066523607814674]], 'Vnn': [0.10581237371666792, 0.0032368011577309395,
-0.027413371793672342, 0.015869286542962588, 0.027265956013929192,
0.010042601420795485], 'Voo': [-0.10552227151792866, 0.05790719018261943,
-0.08141353578942945, 0.003930130259594829, 0.01006758531249384,
0.009981448331701769], 'ao': [0.06116739717004869, 0.10794830765204166,
0.08915150248655652], 'logdEo': [0.33845814370750404, -1.366437774186797,
-1.0541868450784229]}

for k in P0_TEST:
    P0_TEST[k] = array(P0_TEST[k])

if __name__ == '__main__':
    main()

