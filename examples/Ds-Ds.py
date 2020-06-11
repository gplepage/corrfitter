from __future__ import print_function # makes this work for python2 and python3

import collections
import h5py
import gvar as gv
import corrfitter as cf

SHOWPLOTS = False        # display plots at end?
SVDCUT = 0.002

try:
    import matplotlib
except ImportError:
    SHOWPLOTS = False

def main():
    data = make_data('Ds-Ds.h5')
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
        fit.show_plots(save='Ds-Ds.{}.png', view='ratio')

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
        n=2, dataset=h5py.File('Ds-Ds.h5', 'r'), p_exact=fit.pmean
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
    for k in ['a', 'dE', 'Vnn']:
        ans[k] = p[k][0]
    return ans

def make_data(filename):
    """ Read data from file and average it. """
    return gv.dataset.avg_data(h5py.File(filename, 'r'))

def make_models():
    """ Create models to fit data. """
    tmin = 3
    tp = 64
    models = [
        cf.Corr2(
            datatag='Ds', tp=tp, tmin=tmin,
            a=('a', 'ao'), b=('a', 'ao'), dE=('dE', 'dEo'), s=(1., -1.),
            ),
        cf.Corr3(
            datatag='DsDsT18', T=18, tmin=tmin,
            a=('a', 'ao'), dEa=('dE', 'dEo'), sa=(1., -1),
            b=('a', 'ao'), dEb=('dE', 'dEo'), sb=(1., -1.),
            Vnn='Vnn', Voo='Voo', Vno='Vno', symmetric_V=True,
            ),
        cf.Corr3(
            datatag='DsDsT15', T=15, tmin=tmin,
            a=('a', 'ao'), dEa=('dE', 'dEo'), sa=(1., -1),
            b=('a', 'ao'), dEb=('dE', 'dEo'), sb=(1., -1.),
            Vnn='Vnn', Voo='Voo', Vno='Vno', symmetric_V=True,
            )
        ]
    return models

def make_prior(N):
    """ Create priors for fit parameters. """
    prior = gv.BufferDict()
    # Ds
    mDs = gv.gvar('1.2(2)')
    prior['log(a)'] = gv.log(gv.gvar(N * ['0.3(3)']))
    prior['log(dE)'] =  gv.log(gv.gvar(N * ['0.5(5)']))
    prior['log(dE)'][0] = gv.log(mDs)

    # Ds -- oscillating part
    prior['log(ao)'] = gv.log(gv.gvar(N * ['0.1(1)']))
    prior['log(dEo)'] = gv.log(gv.gvar(N * ['0.5(5)']))
    prior['log(dEo)'][0] = gv.log(mDs + gv.gvar('0.3(3)'))

    # V
    nV = int((N * (N + 1)) / 2)
    prior['Vnn'] = gv.gvar(nV * ['0.0(5)'])
    prior['Voo'] = gv.gvar(nV * ['0.0(5)'])
    prior['Vno'] = gv.gvar(N * [N * ['0.0(5)']])
    return prior

def print_results(fit, prior, data):
    """ Print results of fit. """
    outputs = collections.OrderedDict()
    outputs['mDs'] = fit.p['dE'][0]
    outputs['Vnn'] = fit.p['Vnn'][0]

    inputs = collections.OrderedDict()
    inputs['statistics'] = data             # statistical errors in data
    inputs['Ds priors'] = {
        k:prior[k] for k in ['log(a)', 'log(dE)', 'log(ao)', 'log(dEo)']
        }
    inputs['V priors'] = {
        k:prior[k] for k in ['Vnn', 'Vno', 'Voo']
        }

    print('\n' + gv.fmt_values(outputs))
    print(gv.fmt_errorbudget(outputs, inputs))


if __name__ == '__main__':
    gv.ranseed(1234)
    main()
