
#!/usr/bin/env python
# encoding: utf-8
"""
bootstrap.py -- tests bootstrap and marginalization

Created by Peter Lepage on 2010-11-26.
Copyright (c) 2010-2013 G. Peter Lepage.
"""

from __future__ import print_function   # makes this work for python2 and 3

import os
from corrfitter import Corr2,Corr3,CorrFitter
from gvar import gvar,log,exp,ranseed,BufferDict
from gvar.dataset import Dataset,avg_data,bootstrap_iter
from numpy import array,arange,sum

import lsqfit


NBOOTSTRAP = 4              # number of bootstraps: 4 okay; 50 really good

DISPLAYPLOTS = False        # display plots at end of fitting

try: 
    import matplotlib
except ImportError:
    DISPLAYPLOTS = False
#

TEST = True         # testing mode? (True, False, or "dump")

P0 = {}

if TEST:
    TEST_FILENAME = 'bootstrap.testp'
    try:
        with open(TEST_FILENAME,"r") as f:
            P0[True] = BufferDict.load(f, use_json=True)
    except (IOError, EOFError):
        P0[True] = None
    P0[False] = P0[True]
else:
    P0[True] = dict([('log(etas:a)', array([-1.52132077])), ('log(etas:dE)', array([-0.87652893])), ('log(Ds:a)', array([-1.5379696])), ('log(Ds:dE)', array([ 0.18376377])), ('log(Ds:ao)', array([-2.59384428])), ('log(Ds:dEo)', array([ 0.37794027])), ('Vnn', array([[ 0.76905768]])), ('Vno', array([[-0.7629936]]))])
    P0[False] = dict([('log(etas:a)', array([-1.52128644])), ('log(etas:dE)', array([-0.8765244])), ('log(Ds:a)', array([-1.53828898])), ('log(Ds:dE)', array([ 0.18373669])), ('log(Ds:ao)', array([-2.63514295])), ('log(Ds:dEo)', array([ 0.3728692])), ('Vnn', array([[ 0.76630814]])), ('Vno', array([[-0.71139104]]))])
    

def main():
    dfile = "example.data"  # data file
    dset = Dataset(dfile)
    data = avg_data(dset)
    prior = build_prior(6)
    nexp = 1
    ratio = True
    ## original fit ##
    fitter = CorrFitter(models=build_models(),nterm=(nexp,nexp),ratio=ratio)
    print('========================== nexp =',nexp,'  ratio =',ratio)
    fit = fitter.lsqfit(data=data,prior=prior,p0=P0[ratio])
    print_results(fit)
    print('\n')
    if TEST == "dump":
        with open(TEST_FILENAME, "w") as f:
            fit.pmean.dump(f, use_json=True)
    bootstrap_last_fit(fitter,dset,NBOOTSTRAP)
    ##
    ## bootstrap of fit ##
    if DISPLAYPLOTS:
        fitter.display_plots()
##

def bootstrap_last_fit(fitter,dset,nbootstrap):
    fit = fitter.fit        # last fit done by fitter
    ranseed((1950,2000))
    print('\nNumber of bootstrap iterations =',nbootstrap,'\n')
    bs_output = Dataset()    # results accumulated here
    bs_datalist = (avg_data(d) for d in bootstrap_iter(dset,NBOOTSTRAP))
    for bs_fit in fitter.bootstrap_iter(bs_datalist):
        p = bs_fit.pmean
        bs_output.append('etas:dE',exp(p['log(etas:dE)']))
        bs_output.append('etas:a',exp(p['log(etas:a)']))
        bs_output.append('Ds:dE',exp(p['log(Ds:dE)']))
        bs_output.append('Ds:a',exp(p['log(Ds:a)']))
        bs_output.append('Vnn',p['Vnn'])
        bs_output.append('Vno',p['Vno'])
    bs_med = avg_data(bs_output,median=True,spread=True)
    bs_gvar = avg_data(bs_output,spread=True)
    print("%10s  %20s  %20s  %20s"%('label','bs-median','bs-gvar','fit'))
    print(76*'-')
    for lbl in ['etas:dE','etas:a','Ds:dE','Ds:a']:
        loglbl = 'log('+lbl+')'
        print("%10s  %20s  %20s  %20s" % (lbl,fmtlist(bs_med[lbl]),
                                    fmtlist(bs_gvar[lbl]),
                                    fmtlist(exp(fit.p[loglbl]))))
    lbl = 'Vnn'
    print("%10s  %20s  %20s  %20s" % (lbl,fmtlist(bs_med[lbl][0]),
                            fmtlist(bs_gvar[lbl][0]),fmtlist(fit.p[lbl][0])))
    lbl = 'Vno'
    print("%10s  %20s  %20s  %20s" % (lbl,fmtlist(bs_med[lbl][0]),
                            fmtlist(bs_gvar[lbl][0]),fmtlist(fit.p[lbl][0])))
##
      
def print_results(fit):
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
    main()

