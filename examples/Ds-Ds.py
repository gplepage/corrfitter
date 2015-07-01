from __future__ import print_function # makes this work for python2 and python3

import collections
import gvar as gv 
import numpy as np 
import corrfitter as cf 

DISPLAYPLOTS = False        # display plots at end?
try:
    import matplotlib
except ImportError:
    DISPLAYPLOTS = False

def main():
    data = make_data('Ds-Ds.data') 
    fitter = cf.CorrFitter(models=make_models())
    p0 = None
    for N in [1, 2, 3, 4]:  
        print(30 * '=', 'nterm =', N)
        prior = make_prior(N)               
        fit = fitter.lsqfit(data=data, prior=prior, p0=p0) 
        p0 = fit.pmean
    print_results(fit, prior, data) 
    if DISPLAYPLOTS:
        fitter.display_plots()

def make_data(datafile):
    """ Read data from datafile and average it. """
    return gv.dataset.avg_data(cf.read_dataset(datafile))

def make_models():
    """ Create models to fit data. """
    tmin = 3
    tp = 64
    models = [
        cf.Corr2(
            datatag='Ds',
            tp=tp, tdata=range(tp), tfit=range(tmin, tp-tmin+1),
            a=('a', 'ao'), b=('a', 'ao'), dE=('dE', 'dEo'), s=(1., -1.),
            ),
        cf.Corr3(
            datatag='DsDsT18', 
            tdata=range(19), T=18, tfit=range(tmin, 18 - tmin + 1),
            a=('a', 'ao'), dEa=('dE', 'dEo'), tpa=tp, sa=(1., -1),
            b=('a', 'ao'), dEb=('dE', 'dEo'), tpb=tp, sb=(1., -1.),
            Vnn='Vnn', Voo='Voo', Vno='Vno', Von='Vno', symmetric_V=True,
            ),
        cf.Corr3(
            datatag='DsDsT15',
            tdata=range(16), T=15, tfit=range(tmin, 15-tmin +1),
            a=('a', 'ao'), dEa=('dE', 'dEo'), tpa=tp, sa=(1., -1),
            b=('a', 'ao'), dEb=('dE', 'dEo'), tpb=tp, sb=(1., -1.),
            Vnn='Vnn', Voo='Voo', Vno='Vno', Von='Vno', symmetric_V=True,
            )         
        ]
    return models

def make_prior(N):
    """ Create priors for fit parameters. """
    prior = gv.BufferDict()
    # Ds 
    mDs = gv.gvar('1.2(2)')
    prior['loga'] = gv.log(gv.gvar(N * ['0.3(3)']))
    prior['logdE'] =  gv.log(gv.gvar(N * ['0.5(5)']))
    prior['logdE'][0] = gv.log(mDs)

    # Ds -- oscillating part
    prior['logao'] = gv.log(gv.gvar(N * ['0.1(1)']))
    prior['logdEo'] = gv.log(gv.gvar(N * ['0.5(5)']))
    prior['logdEo'][0] = gv.log(mDs + gv.gvar('0.3(3)'))

    # V
    nV = int((N * (N + 1)) / 2)
    prior['Vnn'] = gv.gvar(nV * ['0(1)'])
    prior['Voo'] = gv.gvar(nV * ['0(1)'])
    prior['Vno'] = gv.gvar(N * [N * ['0(1)']])
    return prior

def print_results(fit, prior, data):
    """ Print results of fit. """
    outputs = collections.OrderedDict()
    outputs['mDs'] = fit.p['dE'][0]
    outputs['Vnn'] = fit.p['Vnn'][0]

    inputs = collections.OrderedDict()      
    inputs['statistics'] = data             # statistical errors in data
    inputs.update(prior)                    # errors from priors
    inputs['svd'] = fit.svdcorrection       # errors from svd cut (if present)

    print('\n' + gv.fmt_values(outputs))
    print(gv.fmt_errorbudget(outputs, inputs))

if __name__ == '__main__':
    main()