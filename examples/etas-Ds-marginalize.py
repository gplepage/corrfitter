from __future__ import print_function   # makes this work for python2 and 3

import gvar as gv
import h5py
import numpy as np
import collections
from corrfitter import CorrFitter, Corr2, Corr3, read_dataset

DISPLAYPLOTS = False         # display plots at end of fitting?
try:
    import matplotlib
except ImportError:
    DISPLAYPLOTS = False

def main():
    data = make_data('etas-Ds.h5')
    fitter = CorrFitter(models=make_models())
    p0 = None
    for N in [1, 2]:                                                        # 1
        print(30 * '=', 'nterm =', N)
        prior = make_prior(8)                                               # 2
        fit = fitter.lsqfit(data=data, prior=prior, p0=p0, nterm=(N, N))    # 3
        print(fit)                                                          # 4
        p0 = fit.pmean
    print_results(fit, prior, data)
    if DISPLAYPLOTS:
        fit.show_plots()
    test_fit(
        fitter=fitter, p_exact=fit.pmean, prior=prior, datafile='etas-Ds.h5'
        )

def test_fit(fitter, p_exact, prior, datafile):
    """ Test the fit with simulated data """
    gv.ranseed(9876)
    print('\nRandom seed:', gv.ranseed.seed)
    dataset = h5py.File(datafile)
    for spdata in fitter.simulated_pdata_iter(
        n=2, dataset=dataset, p_exact=p_exact
        ):
        print('\n============================== simulation')
        sfit = fitter.lsqfit(pdata=spdata, prior=prior, p0=p_exact, nterm=(2,2))
        print(sfit.format(pstyle=None))
        # check chi**2 for key parameters
        diff = {}
        for k in ['etas:a', 'etas:dE', 'Ds:a', 'Ds:dE', 'Vnn']:
            p_k = sfit.p[k].flat[0]
            pex_k = p_exact[k].flat[0]
            print(
                '{:>10}:  fit = {:<11}    exact = {:<9.5}    diff = {}'
                    .format(k, p_k, pex_k, p_k - pex_k)
                )

            diff[k] = p_k - pex_k
        print('\nAccuracy of key parameters: ' + gv.fmt_chi2(gv.chi2(diff)))

# reuse code from etas-Ds.py
import importlib
etas_Ds = importlib.import_module('etas-Ds')
make_data = etas_Ds.make_data
make_models = etas_Ds.make_models
make_prior = etas_Ds.make_prior
print_results = etas_Ds.print_results

if __name__ == '__main__':
    main()