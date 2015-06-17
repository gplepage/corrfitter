from __future__ import print_function   # makes this work for python2 and 3

import gvar as gv
import numpy as np
import collections
from corrfitter import CorrFitter, Corr2, Corr3, read_dataset

DISPLAYPLOTS = False         # display plots at end of fitting?
try: 
    import matplotlib
except ImportError:
    DISPLAYPLOTS = False

def main():
    data = make_data('etas-Ds.data') 
    models = make_models()   
    prior = make_prior(8)
    fitter = CorrFitter(models=make_models(), ratio=False)                  # 1
    p0 = None
    for N in [1, 2]:                                                        # 2
        print(30 * '=', 'nterm =', N)
        prior = make_prior(8)                                               # 3          
        fit = fitter.lsqfit(data=data, prior=prior, p0=p0, nterm=(N, N))    # 4 
        p0 = fit.pmean
    print_results(fit, prior, data)  
    if DISPLAYPLOTS:
        fitter.display_plots()
    test_fit(fitter, 'etas-Ds.data')

def test_fit(fitter, datafile):
    """ Test the fit with simulated data """ 
    gv.ranseed((623738625, 435880512, 1745400596))
    print('\nRandom seed:', gv.ranseed.seed)
    dataset = read_dataset(datafile)
    pexact = fitter.fit.pmean
    prior = fitter.fit.prior
    for sdata in fitter.simulated_data_iter(n=2, dataset=dataset, pexact=pexact):
        print('\n============================== simulation')
        sfit = fitter.lsqfit(data=sdata, prior=prior, p0=pexact, nterm=(2, 2))
        diff = []
        # check chi**2 for leading parameters
        for k in sfit.p: 
            diff.append(sfit.p[k].flat[0] - pexact[k].flat[0])
        chi2_diff = gv.chi2(diff)
        print(
            'Leading parameter chi2/dof [dof] = %.2f' % 
            (chi2_diff / chi2_diff.dof),
            '[%d]' % chi2_diff.dof, 
            '  Q = %.1f' % chi2_diff.Q
            )

# reuse code from etas-Ds.py
import importlib
etas_Ds = importlib.import_module('etas-Ds')
make_data = etas_Ds.make_data
make_models = etas_Ds.make_models
make_prior = etas_Ds.make_prior
print_results = etas_Ds.print_results

if __name__ == '__main__':
    main()