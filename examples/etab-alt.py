from __future__ import print_function   # makes this work for python2 and 3

import collections
import gvar as gv
import numpy as np
import corrfitter as cf

DISPLAYPLOTS = False         # display plots at end of fits?

SOURCES = ['l', 'g', 'd', 'e']
EIG_SOURCES = ['0', '1', '2', '3']                                      # 1a
KEYFMT = '1s0.{s1}{s2}'
TDATA = range(1, 24)
SVDCUT = 0.06                                                           # 2

try:
    import matplotlib
except ImportError:
    DISPLAYPLOTS = False

def main():
    data, basis = make_data('etab.h5')
    fitter = cf.CorrFitter(models=make_models())
    p0 = None
    for N in range(1, 8):
        print(30 * '=', 'nterm =', N)
        prior = make_prior(N, basis)
        fit = fitter.lsqfit(data=data, prior=prior, p0=p0, svdcut=SVDCUT)
        print(fit.format(pstyle=None if N < 7 else 'm'))
        p0 = fit.pmean
    print_results(fit, basis, prior, data)
    if DISPLAYPLOTS:
        fitter.display_plots()
    print('\n==================== add svd, prior noise')
    noisy_fit = fitter.lsqfit(
        data=data, prior=prior, p0=fit.pmean, svdcut=SVDCUT,
        noise=True,
        )
    print(noisy_fit.format(pstyle=None))
    dE = fit.p['etab.dE'][:3]
    noisy_dE = noisy_fit.p['etab.dE'][:3]
    print('      dE:', dE)
    print('noisy dE:', noisy_dE)
    print('          ', gv.fmt_chi2(gv.chi2(dE - noisy_dE)))

def make_data(filename):
    data = gv.dataset.avg_data(cf.read_dataset(filename))
    basis = cf.EigenBasis(
        data, keyfmt='1s0.{s1}{s2}', srcs=['l', 'g', 'd', 'e'],
        t=(1,2), tdata=range(1,24),
        )
    return basis.apply(data, keyfmt='1s0.{s1}{s2}'), basis              # 3

def make_models():
    models = []
    for i, s1 in enumerate(EIG_SOURCES):                                # 1b
        for s2 in EIG_SOURCES[i:]:                                      # 1c
            tfit = TDATA if s1 == s2 else TDATA[:12]
            otherdata = None if s1 == s2 else KEYFMT.format(s1=s2, s2=s1)
            models.append(
                cf.Corr2(
                    datatag=KEYFMT.format(s1=s1, s2=s2),
                    tdata=TDATA, tfit=tfit,
                    a='etab.' + s1, b='etab.' + s2, dE='etab.dE',
                    otherdata=otherdata,
                    )
                )
    return models

def make_prior(N, basis):
    return basis.make_prior(nterm=N, keyfmt='etab.{s1}', eig_srcs=True) # 4

def print_results(fit, basis, prior, data):
    print(30 * '=', 'Results\n')
    print(basis.tabulate(fit.p, keyfmt='etab.{s1}'))
    print(basis.tabulate(fit.p, keyfmt='etab.{s1}', eig_srcs=True))
    E = np.cumsum(fit.p['etab.dE'])
    outputs = collections.OrderedDict()
    outputs['a*E(2s-1s)'] = E[1] - E[0]
    outputs['a*E(3s-1s)'] = E[2] - E[0]
    outputs['E(3s-1s)/E(2s-1s)'] = (E[2] - E[0]) / (E[1] - E[0])
    inputs = collections.OrderedDict()
    inputs['prior'] = prior
    inputs['data'] = data
    inputs['svdcut'] = fit.correction
    print(gv.fmt_values(outputs))
    print(gv.fmt_errorbudget(outputs, inputs, colwidth=18))

if __name__ == '__main__':
    gv.ranseed(1)
    main()