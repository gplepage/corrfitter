from __future__ import print_function   # makes this work for python2 and 3

import collections
import gvar as gv
import numpy as np
import corrfitter as cf

SHOWPLOTS = True
SVDCUT = 8e-5

def main():
    data = make_data('etas-Ds.h5')
    fitter = cf.CorrFitter(models=make_models())
    p0 = None
    for N in [1, 2, 3, 4]:
        print(30 * '=', 'nterm =', N)
        prior = make_prior(N)
        fit = fitter.lsqfit(data=data, prior=prior, p0=p0, svdcut=SVDCUT)
        print(fit.format(pstyle=None if N < 4 else 'v'))
        p0 = fit.pmean
    print_results(fit, prior, data)
    if SHOWPLOTS:
        fit.show_plots()

    # check fit quality by adding noise
    print('\n==================== add svd, prior noise')
    noisy_fit = fitter.lsqfit(
        data=data, prior=prior, p0=fit.pmean, svdcut=SVDCUT,
        noise=True,
        )
    print(noisy_fit.format(pstyle=None))
    p = key_parameters(fit.p)
    noisy_p = key_parameters(noisy_fit.p)
    print('      fit:', p)
    print('noisy fit:', noisy_p)
    print('          ', gv.fmt_chi2(gv.chi2(p - noisy_p)))

    # simulated fit
    for sim_pdata in fitter.simulated_pdata_iter(
        n=2, dataset=cf.read_dataset('etas-Ds.h5'), p_exact=fit.pmean
        ):
        print('\n==================== simulation')
        sim_fit = fitter.lsqfit(
            pdata=sim_pdata, prior=prior, p0=fit.pmean, svdcut=SVDCUT,
            )
        print(sim_fit.format(pstyle=None))
        p = key_parameters(fit.pmean)
        sim_p = key_parameters(sim_fit.p)
        print('simulated - exact:', sim_p - p)
        print('          ', gv.fmt_chi2(gv.chi2(p - sim_p)))

def key_parameters(p):
    """ collect key fit parameters in dictionary """
    ans = gv.BufferDict()
    for k in ['etas:a', 'etas:dE', 'Ds:a', 'Ds:dE']:
        ans[k] = p[k][0]
    ans['Vnn'] = p['Vnn'][0, 0]
    return ans

def make_data(datafile):
    """ Read data from datafile and average it. """
    dset = cf.read_dataset(datafile)
    return gv.dataset.avg_data(dset)

def make_models():
    """ Create models to fit data. """
    tmin = 5
    tp = 64
    models = [
        cf.Corr2(
            datatag='etas', tp=tp,  tmin=tmin,
            a='etas:a',  b='etas:a',  dE='etas:dE',
            ),

        cf.Corr2(
            datatag='Ds', tp=tp,  tmin=tmin,
            a=('Ds:a', 'Dso:a'), b=('Ds:a', 'Dso:a'), dE=('Ds:dE', 'Dso:dE'),
            ),

        cf.Corr3(
            datatag='3ptT15', T=15, tmin=tmin, a='etas:a', dEa='etas:dE',
            b=('Ds:a', 'Dso:a'), dEb=('Ds:dE', 'Dso:dE'),
            Vnn='Vnn', Vno='Vno',
            ),

        cf.Corr3(
            datatag='3ptT16', T=16, tmin=tmin, a='etas:a', dEa='etas:dE',
            b=('Ds:a', 'Dso:a'), dEb=('Ds:dE', 'Dso:dE'),
            Vnn='Vnn', Vno='Vno',
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
    prior['Vnn'] = gv.gvar(N * [N * ['0(1)']])
    prior['Vno'] = gv.gvar(N * [N * ['0(1)']])
    return prior

def print_results(fit, prior, data):
    """ Report best-fit results. """
    print('Fit results:')
    p = fit.p                       # best-fit parameters

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
    print('\n  etas->V->Ds  =', Vnn[0, 0])

    # error budget
    outputs = collections.OrderedDict()
    outputs['metas'] = E_etas[0]
    outputs['mDs'] = E_Ds[0]
    outputs['Vnn'] = Vnn[0, 0]

    inputs = collections.OrderedDict()
    inputs['statistics'] = data         # statistical errors in data
    inputs['svd'] = fit.correction
    inputs.update(prior)                # all entries in prior

    print('\n' + gv.fmt_values(outputs))
    print(gv.fmt_errorbudget(outputs, inputs))
    print('\n')

import sys
if sys.argv[1:]:
    SHOWPLOTS = eval(sys.argv[1]) # show plots at end of fitting?

if SHOWPLOTS:
    try:
        import matplotlib
    except ImportError:
        SHOWPLOTS = False

if __name__ == '__main__':
    gv.ranseed(123456)
    main()
    # if True:
    #     main()
    # else:
    #     import cProfile, pstats, StringIO
    #     pr = cProfile.Profile()
    #     pr.enable()
    #     main()
    #     pr.disable()
    #     s = StringIO.StringIO()
    #     sortby = 'tottime'
    #     ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    #     ps.print_stats()
    #     print (s.getvalue())
