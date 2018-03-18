from __future__ import print_function   # makes this work for python2 and 3

import collections
import gvar as gv
import numpy as np
import corrfitter as cf

SHOWPLOTS = False         # display plots at end of fits?
try:
    import matplotlib
except ImportError:
    SHOWPLOTS = False

def main():
    data, basis = make_data('etab.h5')
    fitter = cf.CorrFitter(models=make_models(basis))
    p0 = None
    for N in range(1, 8):
        print(30 * '=', 'nterm =', N)
        prior = make_prior(N, basis)
        fit = fitter.lsqfit(data=data, prior=prior, p0=p0, svdcut=0.0004)
        print(fit.format(pstyle=None if N < 7 else 'v'))
        p0 = fit.pmean
    print_results(fit, basis, prior, data)
    if SHOWPLOTS:
        fit.show_plots(save='etab.{}.png', view='ratio')

def make_data(filename):
    data = gv.dataset.avg_data(cf.read_dataset(filename, grep='1s0'))
    basis = cf.EigenBasis(
        data, keyfmt='1s0.{s1}{s2}', srcs=['l', 'g', 'd', 'e'],
        t=(1, 2), tdata=range(1, 24),
        )
    return data, basis

def make_models(basis):
    models = []
    for s1 in basis.srcs:
        for s2 in basis.srcs:
            tfit = basis.tdata if s1 == s2 else basis.tdata[:14]
            models.append(
                cf.Corr2(
                    datatag=basis.keyfmt.format(s1=s1, s2=s2),
                    tdata=basis.tdata, tfit=tfit,
                    a='etab.' + s1, b='etab.' + s2, dE='etab.dE',
                    )
                )
    return models

def make_prior(N, basis):
    return basis.make_prior(nterm=N, keyfmt='etab.{s1}')

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
    inputs['svdcut'] = fit.svdcorrection
    print(gv.fmt_values(outputs))
    print(gv.fmt_errorbudget(outputs, inputs, colwidth=18))
    print('Prior:\n')
    for k in ['etab.l', 'etab.g', 'etab.d', 'etab.e']:
        print(
            '{:13}{}'.format(k, list(prior[k]))
            )
    print()
    prior_eig = basis.apply(prior, keyfmt='etab.{s1}')
    for k in ['etab.0', 'etab.1', 'etab.2', 'etab.3']:
        print(
            '{:13}{}'.format(k, list(prior_eig[k]))
            )

if __name__ == '__main__':
    if True:
        main()
    else:
        import cProfile, pstats, StringIO
        pr = cProfile.Profile()
        pr.enable()
        main()
        pr.disable()
        s = StringIO.StringIO()
        sortby = 'tottime'
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        print (s.getvalue())
