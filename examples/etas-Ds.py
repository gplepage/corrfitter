from __future__ import print_function   # makes this work for python2 and 3

import collections
import gvar as gv
import numpy as np
from corrfitter import CorrFitter, Corr2, Corr3, read_dataset

DISPLAYPLOTS = True         # display plots at end of fitting?
try: 
    import matplotlib
except ImportError:
    DISPLAYPLOTS = False

def main():
    data = make_data('etas-Ds.data') 
    fitter = CorrFitter(models=make_models())
    p0 = None
    for N in [1, 2, 3, 4]:  
        print(30 * '=', 'nterm =', N)
        prior = make_prior(N)               
        fit = fitter.lsqfit(data=data, prior=prior, p0=p0) 
        p0 = fit.pmean
    print_results(fit, prior, data) 
    if DISPLAYPLOTS:
        fitter.display_plots()
    test_fit(fitter, 'etas-Ds.data') 

def test_fit(fitter, datafile):
    """ Test the fit with simulated data """
    gv.ranseed((1487942813, 775399747, 906327435))
    print('\nRandom seed:', gv.ranseed.seed)
    dataset = read_dataset(datafile)
    pexact = fitter.fit.pmean
    prior = fitter.fit.prior
    for sdata in fitter.simulated_data_iter(n=2, dataset=dataset, pexact=pexact):
        print('\n============================== simulation')
        sfit = fitter.lsqfit(data=sdata, prior=prior, p0=pexact)
        diff = []
        # check chi**2 for leading parameters
        for k in sfit.p: 
            diff.append(sfit.p[k].flat[0] - pexact[k].flat[0])
        chi2diff = gv.chi2(diff)
        print(
            'Leading parameter chi2/dof [dof] = %.2f' % 
            (chi2diff / chi2diff.dof),
            '[%d]' % chi2diff.dof, 
            '  Q = %.1f' % chi2diff.Q
            ) 

def make_data(datafile):
    """ Read data from datafile and average it. """
    return gv.dataset.avg_data(read_dataset(datafile))

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
            a=('Ds:a', 'Dso:a'), b=('Ds:a', 'Dso:a'), 
            dE=('Ds:dE', 'Dso:dE'), s=(1., -1.)
            ),

        Corr3(
            datatag='3ptT15', tdata=range(16), T=15, tfit=range(tmin, 16-tmin), 
            a='etas:a', dEa='etas:dE', tpa=tp, 
            b=('Ds:a', 'Dso:a'), dEb=('Ds:dE', 'Dso:dE'), tpb=tp, sb=(1, -1.), 
            Vnn='Vnn', Vno='Vno'
            ), 
            
        Corr3(
            datatag='3ptT16', tdata=range(17), T=16, tfit=range(tmin, 17-tmin), 
            a='etas:a', dEa='etas:dE', tpa=tp, 
            b=('Ds:a', 'Dso:a'), dEb=('Ds:dE', 'Dso:dE'), tpb=tp, sb=(1, -1.), 
            Vnn='Vnn', Vno='Vno'
            )
        ]
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
    prior['log(Dso:a)'] = gv.log(gv.gvar(N * ['0.1(1)']))
    prior['log(Dso:dE)'] = gv.log(gv.gvar(N * ['0.5(5)']))
    prior['log(Dso:dE)'][0] = gv.log(mDs + gv.gvar('0.3(3)'))

    # V
    prior['Vnn'] = gv.gvar(N * [N * ["0(1)"]])
    prior['Vno'] = gv.gvar(N * [N * ["0(1)"]])
    return prior

def print_results(fit, prior, data):
    """ Report best-fit results. """
    print('Fit results:')
    p = fit.transformed_p                       # best-fit parameters

    # etas
    E_etas = np.cumsum(p['etas:dE'])
    a_etas = p['etas:a']
    print('  Eetas:', E_etas[:3])
    print('  aetas:', a_etas[:3])

    # Ds
    E_Ds = np.cumsum(p['Ds:dE'])
    a_Ds = p['Ds:a']
    print('\n  EDs:', E_Ds[:3])
    print(  '  aDs:', a_Ds[:3])

    # Dso -- oscillating piece
    E_Dso = np.cumsum(p['Dso:dE'])
    a_Dso = p['Dso:a']
    print('\n  EDso:', E_Dso[:3])
    print(  '  aDso:', a_Dso[:3])

    # V
    Vnn = p['Vnn']
    Vno = p['Vno']
    print('\n  etas->V->Ds  =', Vnn[0, 0].fmt())
    print('  etas->V->Dso =', Vno[0, 0].fmt())

    # error budget
    outputs = collections.OrderedDict()
    outputs['metas'] = E_etas[0]
    outputs['mDs'] = E_Ds[0]
    outputs['mDso-mDs'] = E_Dso[0] - E_Ds[0]
    outputs['Vnn'] = Vnn[0, 0]
    outputs['Vno'] = Vno[0, 0]

    inputs = collections.OrderedDict()
    inputs['statistics'] = data                 # statistical errors in data
    inputs.update(prior)                        # all entries in prior
    inputs['svd'] = fit.svdcorrection           # svd cut (if present)

    print('\n' + gv.fmt_values(outputs))
    print(gv.fmt_errorbudget(outputs, inputs))
    print('\n')


if __name__ == '__main__':
    main()