from __future__ import print_function   # makes this work for python2 and 3

import gvar as gv
import numpy as np
import collections
from corrfitter import CorrFitter, Corr2, Corr3

DISPLAYPLOTS = False         # display plots at end of fitting
try: 
    import matplotlib
except ImportError:
    DISPLAYPLOTS = False

POLISH = True               # polish solution at end

STRUCTURE_MODELS = True

def main():
    data = make_data('example.data') 
    models = make_models()   
    prior = make_prior(8)
    fitter = CorrFitter(models=make_models())
    p0 = None
    for N in [1, 2, 3, 4]:  
        print(30 * '=', 'nterm =', N)
        prior = make_prior(N)               
        fit = fitter.chained_lsqfit(data=data, prior=prior, p0=p0) 
        p0 = fit.pmean
    if POLISH:
        if DISPLAYPLOTS:
                fitter.display_plots()
        fit = fitter.lsqfit(data=data, prior=fit.p, svdcut=1e-4)
    print_results(fit, prior, data)  
    if DISPLAYPLOTS:
            fitter.display_plots()

def make_data(datafile):
    """ Read data from datafile and average it. """
    return gv.dataset.avg_data(gv.dataset.Dataset(datafile))

def make_models():
    """ Create models to fit data. """
    tmin = 5
    tp = 64
    models = [
        Corr2(
            datatag='etas', 
            tp=tp,  tdata=range(tp),  tfit=range(tmin, tp-tmin),  
            a='etas:a',  b='etas:a',  dE='etas:dE'
            ),  
            
        Corr2(
            datatag='Ds',
            tp=tp,  tdata=range(tp),  tfit=range(tmin, tp-tmin),  
            a=('Ds:a', 'Ds:ao'), b=('Ds:a', 'Ds:ao'), 
            dE=('Ds:dE', 'Ds:dEo'), s=(1., -1.)
            ),

        Corr3(
            datatag='3ptT15', tdata=range(16), T=15, tfit=range(tmin, 16-tmin), 
            a='etas:a', dEa='etas:dE', tpa=tp, 
            b=('Ds:a', 'Ds:ao'), dEb=('Ds:dE', 'Ds:dEo'), tpb=tp, sb=(1, -1.), 
            Vnn='Vnn', Vno='Vno'
            ), 
            
        Corr3(
            datatag='3ptT16', tdata=range(17), T=16, tfit=range(tmin, 17-tmin), 
            a='etas:a', dEa='etas:dE', tpa=tp, 
            b=('Ds:a', 'Ds:ao'), dEb=('Ds:dE', 'Ds:dEo'), tpb=tp, sb=(1, -1.), 
            Vnn='Vnn', Vno='Vno'
            )
        ]
    if STRUCTURE_MODELS:
        return [models[:2]] + models[2:]
    else:
        return models

def make_prior(N):
    """ Create priors for fit parameters. """
    prior = gv.BufferDict()
    # etas
    metas = gv.gvar('0.4(2)')
    prior['log(etas:a)'] = gv.log(gv.gvar(N * ['0.3(3)']))   
    prior['log(etas:dE)'] = gv.log(gv.gvar(N * ['0.5(5)']))
    prior['log(etas:dE)'][0] = gv.log(metas)

    # Ds
    mDs = gv.gvar('1.2(2)')
    prior['log(Ds:a)'] = gv.log(gv.gvar(N * ['0.3(3)']))
    prior['log(Ds:dE)'] = gv.log(gv.gvar(N * ['0.5(5)']))
    prior['log(Ds:dE)'][0] = gv.log(mDs)

    # Ds -- oscillating part
    prior['log(Ds:ao)'] = gv.log(gv.gvar(N * ['0.1(1)']))
    prior['log(Ds:dEo)'] = gv.log(gv.gvar(N * ['0.5(5)']))
    prior['log(Ds:dEo)'][0] = gv.log(mDs + gv.gvar('0.3(3)'))

    # V
    prior['Vnn'] = gv.gvar(N * [N * ["0(1)"]])
    prior['Vno'] = gv.gvar(N * [N * ["0(1)"]])

    return prior

def print_results(fit, prior, data):
    """ Report best-fit results. """
    print('Fit results:')
    # etas
    E_etas = np.cumsum(gv.exp(fit.p['log(etas:dE)']))
    a_etas = gv.exp(fit.p['log(etas:a)'])
    print('  Eetas:', E_etas[:3])
    print('  aetas:', a_etas[:3])

    # Ds
    E_Ds = np.cumsum(gv.exp(fit.p['log(Ds:dE)']))
    a_Ds = gv.exp(fit.p['log(Ds:a)'])
    print('\n  EDs:', E_Ds[:3])
    print(  '  aDs:', a_Ds[:3])

    # Dso -- oscillating piece
    E_Dso = np.cumsum(gv.exp(fit.p['log(Ds:dEo)']))
    a_Dso = gv.exp(fit.p['log(Ds:ao)'])
    print('\n  EDso:', E_Dso[:3])
    print(  '  aDso:', a_Dso[:3])

    # V
    Vnn = fit.p['Vnn']
    Vno = fit.p['Vno']
    print('\n  etas->V->Ds  =', Vnn[0, 0].fmt())
    print('  etas->V->Dso =', Vno[0, 0].fmt())

    # error budget
    outputs = gv.BufferDict()
    outputs['metas'] = E_etas[0]
    outputs['mDs'] = E_Ds[0]
    outputs['mDso-mDs'] = E_Dso[0] - E_Ds[0]
    outputs['Vnn'] = Vnn[0, 0]
    outputs['Vno'] = Vno[0, 0]

    inputs = collections.OrderedDict()
    inputs['statistics'] = data             # statistical errors in data
    inputs.update(prior)                    # all entries in prior
    inputs['svd'] = fit.svdcorrection        # svd cut (if present)

    print('\n' + gv.fmt_values(outputs))
    print(gv.fmt_errorbudget(outputs, inputs))

    print('\n')


if __name__ == '__main__':
    main()