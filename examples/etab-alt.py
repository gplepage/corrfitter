from __future__ import print_function   # makes this work for python2 and 3

import collections
import gvar as gv
import numpy as np
from corrfitter import CorrFitter, Corr2, EigenBasis, read_dataset

DISPLAYPLOTS = False         # display plots at end of fits?
try: 
    import matplotlib
except ImportError:
    DISPLAYPLOTS = False

def main():
    data, basis = make_data('etab.data') 
    fitter = CorrFitter(models=make_models(basis))
    p0 = None
    for N in range(1, 10):  
        print(30 * '=', 'nterm =', N)
        prior = make_prior(N, basis)          
        fit = fitter.lsqfit(data=data, prior=prior, p0=p0, svdcut=0.005 )    #1
        p0 = fit.pmean
    print_results(fit, basis, prior, data) 
    if DISPLAYPLOTS:
        fitter.display_plots()

def make_data(filename):
    data = gv.dataset.avg_data(read_dataset(filename))
    basis = EigenBasis(
        data, keyfmt='1s0.{s1}{s2}', srcs=['l', 'g', 'd', 'e'], 
        t=(1,2), tdata=range(1,24), 
        ) 
    return basis.apply(data, keyfmt='1s0.{s1}{s2}'), basis                  #2

def make_models(basis):
    models = []
    for s1 in basis.eig_srcs:                                               #3
        for s2 in basis.eig_srcs:                                           #4
            tfit = basis.tdata if s1 == s2 else basis.tdata[:14]
            models.append(
                Corr2(
                    datatag=basis.keyfmt.format(s1=s1, s2=s2),
                    tdata=basis.tdata, tfit=tfit,
                    a='etab.' + s1, b='etab.' + s2, dE='etab.dE',
                    )
                )
    return models

def make_prior(N, basis):
    return basis.make_prior(nterm=N, keyfmt='etab.{s1}', eig_srcs=True)     #5

def print_results(fit, basis, prior, data):
    print(30 * '=', 'Results\n')
    print(basis.tabulate(fit.transformed_p, keyfmt='etab.{s1}'))
    print(basis.tabulate(fit.transformed_p, keyfmt='etab.{s1}', eig_srcs=True))
    E = np.cumsum(fit.transformed_p['etab.dE']) 
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

if __name__ == '__main__':
    main()